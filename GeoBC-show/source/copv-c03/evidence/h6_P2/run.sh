#!/bin/bash
set -euo pipefail
export PETSC_OPTIONS="-use_gpu_aware_mpi 0 -no_signal_handler"
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONFAULTHANDLER=1
export TMPDIR='/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/.runtime/fenicsx_cuda/tmp'
export XDG_CACHE_HOME='/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/.runtime/fenicsx_cuda/cache'
export PATH='/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/.runtime/fenicsx_cuda/env/bin':"$PATH"
cd '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/publication/GeoBC-Neural-Operator'
exec '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/.runtime/fenicsx_cuda/env/bin/python' -B -m tools.copv_fenicsx_c03 --mesh '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/data/copv/c02_structural_v1/mesh_h6_netgen/mesh.npz' --review '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/data/copv/c02_structural_v1/mesh_h6_netgen/mesh_review.json' --contract '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/publication/GeoBC-Neural-Operator/configs/copv_c03_mesh_v1.json' --output '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/data/copv/c03_mesh_v1/h6_P2' --degree 2 --cap-convention circular --preconditioner p1-gamg
