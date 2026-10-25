#!/bin/bash
set -e
export PETSC_OPTIONS="-use_gpu_aware_mpi 0 -no_signal_handler"
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
export TMPDIR='/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/.runtime/fenicsx_cuda/tmp'
export XDG_CACHE_HOME='/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/.runtime/fenicsx_cuda/cache'
export PATH='/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/.runtime/fenicsx_cuda/env/bin':"$PATH"
export PYTHONPATH='/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/publication/GeoBC-Neural-Operator'
cd '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/publication/GeoBC-Neural-Operator'
exec '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/.runtime/fenicsx_cuda/env/bin/python' '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/data/copv/c04_pilot_v1/solver_schema_validation_P2_attempt002/solver_executed_snapshot.py' --mesh '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/data/copv/c02_structural_v1/mesh_h6_netgen/mesh.npz' --review '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/data/copv/c02_structural_v1/mesh_h6_netgen/mesh_review.json' --surface '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/data/copv/c04_pilot_v1/viewer_packs/validation/surface/surface.npz' --contract '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/data/copv/c04_pilot_v1/solver_schema_validation_P2_attempt002/executed_contract.json' --output '/mnt/c/Users/ASUS/Desktop/GeoBC Neural Operator/data/copv/c04_pilot_v1/solver_schema_validation_P2_attempt002' --degree 2 --preconditioner p1-gamg --cap-convention discrete --research-validation-only --validate-native-schema --save-gp-local
