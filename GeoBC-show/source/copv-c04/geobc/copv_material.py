"""Small-strain 3-D orthotropic elastic utilities, with no production constants.

Voigt order: (11,22,33,23,13,12). Stress uses physical shear stress;
strain uses engineering shear (e11,e22,e33,2e23,2e13,2e12).
Thus sigma_voigt @ epsilon_voigt == sigma_tensor : epsilon_tensor.
E/G/stress units are MPa; compliance is MPa^-1; strain dimensionless.
nu_ij means -epsilon_jj / epsilon_ii under uniaxial sigma_ii.
Axes Q has local unit axes as its COLUMNS in global coordinates:
T_global = Q @ T_local @ Q.T. Only proper orthogonal frames are accepted.
"""
from __future__ import annotations
import numpy as np

VOIGT_ORDER = ('11', '22', '33', '23', '13', '12')
_PAIRS = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))


def engineering_compliance(*, E1, E2, E3, nu12, nu13, nu23, G12, G13, G23):
    """Return symmetric positive-definite S with reciprocal Poisson ratios.

    Only the three supplied major ratios are independent. In particular,
    nu21=nu12*E2/E1, nu31=nu13*E3/E1, nu32=nu23*E3/E2.
    Nonphysical constants (including loss of positive definiteness) raise.
    """
    constants = np.asarray([E1, E2, E3, nu12, nu13, nu23, G12, G13, G23], float)
    if not np.isfinite(constants).all() or np.any(constants[[0, 1, 2, 6, 7, 8]] <= 0):
        raise ValueError('All constants must be finite; E and G must be positive')
    S = np.diag([1/E1, 1/E2, 1/E3, 1/G23, 1/G13, 1/G12])
    S[0, 1] = S[1, 0] = -nu12/E1
    S[0, 2] = S[2, 0] = -nu13/E1
    S[1, 2] = S[2, 1] = -nu23/E2
    require_spd(S)
    return S


def require_spd(matrix):
    """Check symmetry and numerical positive definiteness of 6x6 matrices."""
    matrix = np.asarray(matrix, float)
    if matrix.shape[-2:] != (6, 6) or not np.isfinite(matrix).all():
        raise ValueError('Expected finite (...,6,6) matrix')
    scale = np.max(np.abs(matrix), axis=(-2, -1), keepdims=True)
    if np.any(np.abs(matrix - np.swapaxes(matrix, -1, -2)) > 1e-12 * scale):
        raise ValueError('Constitutive matrix is not symmetric')
    eigenvalues = np.linalg.eigvalsh(matrix)
    if np.any(eigenvalues[..., 0] <= 1e-13 * eigenvalues[..., -1]):
        raise ValueError('Constitutive matrix is not numerically positive definite')
    return eigenvalues


def engineering_stiffness(**constants):
    C = np.linalg.inv(engineering_compliance(**constants))
    require_spd(C)
    return C


def _to_tensor(vector, engineering_strain):
    vector = np.asarray(vector, float)
    if vector.shape[-1:] != (6,):
        raise ValueError('Expected (...,6) Voigt vector')
    out = np.zeros(vector.shape[:-1] + (3, 3))
    for k, (i, j) in enumerate(_PAIRS):
        out[..., i, j] = out[..., j, i] = vector[..., k] / (2 if engineering_strain and i != j else 1)
    return out


def _to_voigt(tensor, engineering_strain):
    tensor = np.asarray(tensor, float)
    if tensor.shape[-2:] != (3, 3) or not np.isfinite(tensor).all():
        raise ValueError('Expected symmetric (...,3,3) tensor')
    scale = np.maximum(np.max(np.abs(tensor), axis=(-2, -1), keepdims=True), 1e-30)
    if np.any(np.abs(tensor - np.swapaxes(tensor, -1, -2)) > 1e-12 * scale):
        raise ValueError('Expected symmetric (...,3,3) tensor')
    return np.stack([tensor[..., i, j] * (2 if engineering_strain and i != j else 1) for i, j in _PAIRS], axis=-1)


def stress_to_tensor(stress):
    return _to_tensor(stress, False)


def strain_to_tensor(strain):
    return _to_tensor(strain, True)


def stress_to_voigt(stress):
    return _to_voigt(stress, False)


def strain_to_voigt(strain):
    return _to_voigt(strain, True)


def validate_axes(axes):
    Q = np.asarray(axes, float)
    if Q.shape[-2:] != (3, 3) or not np.isfinite(Q).all():
        raise ValueError('Expected finite (...,3,3) axes')
    if not np.allclose(np.swapaxes(Q, -1, -2) @ Q, np.eye(3), atol=1e-10, rtol=0) or not np.allclose(np.linalg.det(Q), 1., atol=1e-10, rtol=0):
        raise ValueError('Axes must be orthonormal and right handed')
    return Q


def rotate_stress(stress, axes, *, to_global=True):
    Q = validate_axes(axes)
    if not to_global:
        Q = np.swapaxes(Q, -1, -2)
    return stress_to_voigt(Q @ stress_to_tensor(stress) @ np.swapaxes(Q, -1, -2))


def rotate_strain(strain, axes, *, to_global=True):
    Q = validate_axes(axes)
    if not to_global:
        Q = np.swapaxes(Q, -1, -2)
    return strain_to_voigt(Q @ strain_to_tensor(strain) @ np.swapaxes(Q, -1, -2))


def rotate_stiffness(stiffness_local, axes):
    """Return global engineering stiffness by six unit-strain responses.

    Separate strain and stress transformations preserve engineering shear;
    applying a single stress transform on both sides would be incorrect.
    Supports broadcasting of stiffness/axes batch dimensions.
    """
    C = np.asarray(stiffness_local, float)
    require_spd(C)
    Q = validate_axes(axes)
    columns = []
    for unit_strain in np.eye(6):
        strain_local = rotate_strain(unit_strain, Q, to_global=False)
        stress_local = np.einsum('...ij,...j->...i', C, strain_local)
        columns.append(rotate_stress(stress_local, Q))
    rotated = np.stack(columns, axis=-1)
    require_spd(rotated)
    return rotated


def fiber_axes(meridional, circumferential, outward_normal, alpha_radians):
    """Material 1=fiber, 2=in-plane transverse, 3=outward normal.

    Base (meridional,circumferential,normal) must already be right handed.
    Alpha is measured from meridian toward circumference; no winding-path
    or manufacturability assumption is implied by this geometric rotation.
    """
    base = np.stack(np.broadcast_arrays(meridional, circumferential, outward_normal), axis=-1)
    validate_axes(base)
    alpha = np.asarray(alpha_radians, float)
    rotation = np.zeros(alpha.shape + (3, 3))
    rotation[..., 0, 0] = rotation[..., 1, 1] = np.cos(alpha)
    rotation[..., 1, 0] = np.sin(alpha)
    rotation[..., 0, 1] = -np.sin(alpha)
    rotation[..., 2, 2] = 1.
    return validate_axes(base @ rotation)


def strain_energy_density(strain, stiffness):
    """Return 1/2 eps_eng.T C eps_eng, MPa = N/mm^2."""
    return .5 * np.einsum('...i,...ij,...j->...', strain, stiffness, strain)
