#!/bin/bash
# Build the SHARED dala-dock environment + AutoDock-GPU on NERSC Perlmutter.
# Run ONCE (as the project owner). Everything lands in the project's CFS dir so
# all students share it. Adjust PROJECT_ROOT, then run from the repo root:
#   bash scripts/build_on_perlmutter.sh
#
# NOTE: module names/versions below are the expected Perlmutter ones; confirm with
# `module avail` and tweak if needed. AutoDock-GPU needs gcc in the 9-12 range and
# CUDA >= 11 (Perlmutter A100 = compute capability 8.0 -> TARGETS=80).
set -e

PROJECT_ROOT="${DALADOCK_ROOT:-$CFS/m4833/projects/project<N>}"   # <-- EDIT project number
ENV="$PROJECT_ROOT/env"
REPO="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== [1/3] create shared conda env at $ENV ==="
module load python
if [ ! -x "$ENV/bin/python" ]; then
  conda env create -f "$REPO/environment_nersc.yml" --prefix "$ENV"
else
  echo "  env already exists; skipping"
fi

echo "=== [2/3] build AutoDock-GPU for A100 (sm_80) ==="
module load gpu cudatoolkit gcc
cd "$PROJECT_ROOT"
[ -d AutoDock-GPU ] || git clone --depth 1 https://github.com/ccsb-scripps/AutoDock-GPU.git
cd AutoDock-GPU
CUDA_ROOT="$(dirname "$(dirname "$(command -v nvcc)")")"
export GPU_INCLUDE_PATH="$CUDA_ROOT/include"
export GPU_LIBRARY_PATH="$CUDA_ROOT/lib64"
# NVIDIA HPC SDK (Perlmutter) keeps curand/cublas in a separate math_libs dir.
# Add it so the linker (-lcurand) and runtime can find them.
MATH_LIBS=""
if [ -d "${CUDA_ROOT%/cuda/*}/math_libs/${CUDA_ROOT##*/}/lib64" ]; then
  MATH_LIBS="${CUDA_ROOT%/cuda/*}/math_libs/${CUDA_ROOT##*/}/lib64"
elif [ ! -e "$CUDA_ROOT/lib64/libcurand.so" ]; then
  MATH_LIBS="$(dirname "$(find "${CUDA_ROOT%/cuda/*}" -name 'libcurand.so' 2>/dev/null | head -1)")"
fi
if [ -n "$MATH_LIBS" ] && [ -d "$MATH_LIBS" ]; then
  export LIBRARY_PATH="$MATH_LIBS:$LIBRARY_PATH"
  export LD_LIBRARY_PATH="$MATH_LIBS:$LD_LIBRARY_PATH"
  echo "  math_libs (curand): $MATH_LIBS"
fi
echo "  CUDA_ROOT=$CUDA_ROOT ; host gcc=$(gcc -dumpfullversion)"
make DEVICE=CUDA TARGETS=80 NUMWI=128 OVERLAP=OFF
echo "  built: $PROJECT_ROOT/AutoDock-GPU/bin/autodock_gpu_128wi"

echo "=== [3/3] values for the Jupyter kernel helper ==="
echo "  PROJECT_ROOT = $PROJECT_ROOT"
echo "  DALADOCK_CUDA_LIB = $CUDA_ROOT/lib64${MATH_LIBS:+:$MATH_LIBS}"
echo
echo "Next: put these into kernels/dala-dock/kernel-helper.sh, then register the kernel:"
echo "  mkdir -p \$HOME/.local/share/jupyter/kernels/dala-dock"
echo "  cp $REPO/kernels/dala-dock/kernel.json $REPO/kernels/dala-dock/kernel-helper.sh \\"
echo "     \$HOME/.local/share/jupyter/kernels/dala-dock/"
echo "  chmod +x \$HOME/.local/share/jupyter/kernels/dala-dock/kernel-helper.sh"
echo
echo "Then open notebooks/docking_demo.ipynb on https://jupyter.nersc.gov (Perlmutter GPU),"
echo "select the 'dala-dock' kernel, and Run All."
