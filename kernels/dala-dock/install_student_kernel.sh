#!/bin/bash
# STUDENTS: run this ONCE on Perlmutter to install the "dala-dock" Jupyter kernel.
# It points at the shared Project 2 environment on CFS - you install nothing else.
# After it runs, refresh https://jupyter.nersc.gov and pick the "dala-dock" kernel.
#
#   bash kernels/dala-dock/install_student_kernel.sh
set -e

# --- shared Project 2 locations (built once by the project owner) ---
SHARED=/global/cfs/projectdirs/m4388/projects/project2/HPC_Bootcamp_2026
ENV="$SHARED/env"
ADGPU="$SHARED/AutoDock-GPU/bin/autodock_gpu_128wi"
# NVIDIA HPC SDK CUDA runtime + math libs (curand) used to build AutoDock-GPU:
CUDA_LIB="/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/cuda/12.9/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/math_libs/12.9/lib64"

[ -x "$ENV/bin/python" ] || { echo "ERROR: shared env not found at $ENV"; exit 1; }

KDIR="$HOME/.local/share/jupyter/kernels/dala-dock"
mkdir -p "$KDIR"

cat > "$KDIR/kernel.json" <<'JSON'
{
  "argv": ["{resource_dir}/kernel-helper.sh", "-f", "{connection_file}"],
  "display_name": "dala-dock",
  "language": "python"
}
JSON

cat > "$KDIR/kernel-helper.sh" <<EOF
#!/bin/bash
export ADGPU="$ADGPU"
unset PYTHONPATH
export LD_LIBRARY_PATH="$ENV/lib:$CUDA_LIB:\$LD_LIBRARY_PATH"
exec "$ENV/bin/python" -m ipykernel_launcher "\$@"
EOF
chmod +x "$KDIR/kernel-helper.sh"

echo "Installed the 'dala-dock' kernel -> $KDIR"
echo "Now refresh https://jupyter.nersc.gov (Perlmutter GPU node) and select 'dala-dock'."
