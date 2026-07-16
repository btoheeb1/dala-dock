#!/bin/bash
# Jupyter kernel launcher for the SHARED dala-dock env on NERSC Perlmutter.
#
# This runs the shared env's python directly (no conda activate needed) and sets up
# the GPU engine, so every student gets the same working kernel with zero install.
#
# TO DEPLOY (once, on Perlmutter): set PROJECT_ROOT below to your CFS project dir,
# and set CUDA_LIB to the cudatoolkit lib path used to BUILD AutoDock-GPU
# (find it after `module load cudatoolkit` via:  echo $CUDATOOLKIT_HOME/lib64 ,
#  or `dirname $(dirname $(which nvcc)))/lib64`). Then register this kernel:
#   mkdir -p $HOME/.local/share/jupyter/kernels/dala-dock
#   cp kernel.json kernel-helper.sh $HOME/.local/share/jupyter/kernels/dala-dock/
#   chmod +x $HOME/.local/share/jupyter/kernels/dala-dock/kernel-helper.sh
# (For the class-wide shared kernel, the coordinator installs it in the project's
#  shared kernels location instead of $HOME.)

PROJECT_ROOT="${DALADOCK_ROOT:-$CFS/m4833/projects/project<N>}"   # <-- EDIT project number
ENV="$PROJECT_ROOT/env"
CUDA_LIB="${DALADOCK_CUDA_LIB:-}"                                 # <-- set after building AutoDock-GPU

export ADGPU="$PROJECT_ROOT/AutoDock-GPU/bin/autodock_gpu_128wi"  # GPU engine (optional)
unset PYTHONPATH                                                  # avoid leaked-env clashes
export LD_LIBRARY_PATH="$ENV/lib${CUDA_LIB:+:$CUDA_LIB}:$LD_LIBRARY_PATH"

exec "$ENV/bin/python" -m ipykernel_launcher "$@"
