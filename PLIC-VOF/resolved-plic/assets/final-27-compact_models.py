"""Compact pure-neural monotone flux operators for validation-selected Pareto tests.

Inputs: raw VOF stencils[B,k,k] (top row is high y), courants[B].
Output: right-face swept-water / center-water ratio[B]. No geometry is inferred
outside the learned MLP. CDF and capacity arithmetic match the original head.
"""
from typing import Any, Dict, List, Mapping

import torch
from torch import Tensor, nn
from torch.nn import functional as F

ARCHITECTURE = "compact_monotone_flux_v1"


class CompactFlux(nn.Module):
    __constants__ = ["k", "stencil_size", "width", "num_bins", "minimum_mass"]

    def __init__(self, config: Mapping[str, Any]):
        super().__init__()
        if config.get("operator", "flux") != "flux":
            raise ValueError("CompactFlux supports operator='flux' only.")
        if config.get("architecture", ARCHITECTURE) != ARCHITECTURE:
            raise ValueError("Unsupported compact flux architecture.")
        self.k = int(config["stencil_size"])
        self.stencil_size = self.k
        self.width = int(config["width"])
        self.num_bins = int(config.get("num_bins", 16))
        self.minimum_mass = 1.0e-6
        if self.k < 1 or self.k % 2 != 1:
            raise ValueError("stencil_size must be a positive odd integer.")
        if self.width < 1 or self.num_bins < 2:
            raise ValueError("width must be positive and num_bins must be >=2.")
        self.net = nn.Sequential(
            nn.Linear(self.k * self.k, self.width),
            nn.GELU(),
            nn.Linear(self.width, self.width),
            nn.GELU(),
            nn.Linear(self.width, self.num_bins),
        )
        self.register_buffer(
            "bin_starts", torch.arange(self.num_bins, dtype=torch.float32)[None, :]
        )

    def forward(self, stencils: Tensor, courants: Tensor) -> Tensor:
        if stencils.dim() != 3:
            raise ValueError("Expected stencils[B,k,k].")
        if stencils.size(1) != self.k or stencils.size(2) != self.k:
            raise ValueError("Stencil shape does not match model configuration.")
        batch = stencils.size(0)
        if courants.numel() != batch:
            raise ValueError("Expected one Courant value per stencil.")
        reflected = torch.flip(stencils, [1])
        features = torch.cat((stencils, reflected), dim=0).reshape(
            2 * batch, self.k * self.k
        )
        logits = self.net(features).reshape(2, batch, self.num_bins)
        # Average positive increments, exactly as in the original row-mean head.
        masses = (F.softplus(logits) + self.minimum_mass).mean(dim=0)
        c = courants.reshape(-1).clamp(0.0, 1.0)
        coverage = (float(self.num_bins) * c[:, None] - self.bin_starts).clamp(
            0.0, 1.0
        )
        ratio = (coverage * masses).sum(dim=1) / masses.sum(dim=1)
        f = stencils[:, self.k // 2, self.k // 2].clamp(0.0, 1.0)
        safe = f.clamp_min(1.0e-12)
        lower = ((f + c - 1.0) / safe).clamp(0.0, 1.0)
        upper = (c / safe).clamp(0.0, 1.0)
        ratio = torch.minimum(torch.maximum(ratio, lower), upper)
        return torch.where(f > 0.0, ratio, torch.zeros_like(ratio))

    @torch.jit.unused
    def get_config(self) -> Dict[str, Any]:
        return dict(
            operator="flux",
            stencil_size=self.k,
            width=self.width,
            num_bins=self.num_bins,
            architecture=ARCHITECTURE,
        )


def build_model(config: Mapping[str, Any]) -> CompactFlux:
    """Construct on CPU; caller controls seed, dtype, device and training."""
    return CompactFlux(config)


def candidate_configs() -> List[Dict[str, Any]]:
    return [
        dict(operator="flux", stencil_size=k, width=w, num_bins=16,
             architecture=ARCHITECTURE)
        for k in (3, 5) for w in (32, 64)
    ]


def parameter_count(config: Mapping[str, Any]) -> int:
    k = int(config["stencil_size"])
    w = int(config["width"])
    b = int(config.get("num_bins", 16))
    return (k * k + 1) * w + (w + 1) * w + (w + 1) * b
