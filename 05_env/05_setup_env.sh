#!/bin/bash
# Step 05: build a self-contained environment for PURGE and install Puffin / Puffin-D.
#
# Everything heavy (conda, caches, model weights, genome) goes to scratch, because it
# can be rebuilt from this script; only the script and the exported environment
# specification live in project_data, which is backed up.
#
# Run on the LOGIN node: compute nodes have no internet access.
# Downloads about 4 GB and takes roughly 20-40 minutes on a first run.
# The script is idempotent: steps that already completed are skipped.
set -euo pipefail

SCRATCH=${PURGE_SCRATCH:?}
PROJ=${PURGE_PROJECT:?}
MF=$SCRATCH/miniforge3
ENV=$SCRATCH/envs/purge
STAMP=$(date +%Y%m%d)

# The home directory sits at its 1 GB quota, so every cache is redirected to scratch.
export CONDA_PKGS_DIRS=$SCRATCH/.conda_pkgs
export PIP_CACHE_DIR=$SCRATCH/.pip_cache
export XDG_CACHE_HOME=$SCRATCH/.cache
export MPLCONFIGDIR=$SCRATCH/.cache/matplotlib
export TMPDIR=$SCRATCH/.tmp
mkdir -p "$TMPDIR"
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$MPLCONFIGDIR" "$SCRATCH/tools" "$SCRATCH/ref" "$PROJ/05_env"

echo "== 1. Miniforge =="
if [ ! -x "$MF/bin/conda" ]; then
    curl -fL -o "$SCRATCH/miniforge.sh" \
        https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash "$SCRATCH/miniforge.sh" -b -p "$MF"
    rm -f "$SCRATCH/miniforge.sh"
fi
source "$MF/etc/profile.d/conda.sh"
conda --version

echo "== 2. conda environment =="
if [ ! -d "$ENV" ]; then
    conda create -y -p "$ENV" -c conda-forge -c bioconda \
        python=3.11 numpy pandas scipy matplotlib pyfaidx pybigwig h5py tqdm
fi
conda activate "$ENV"
python --version

echo "== 3. PyTorch with CUDA 12.4 runtime =="
# Wheels bundle their own CUDA runtime, so only the NVIDIA driver matters; the cluster
# provides driver 550 (module cuda/12.4.1_550.54.15), which supports CUDA 12.4.
if ! python -c "import torch" 2>/dev/null; then
    python -m pip install torch --index-url https://download.pytorch.org/whl/cu124
fi
python -c "import torch; print('torch', torch.__version__, 'cuda build', torch.version.cuda)"

echo "== 4. Puffin repository =="
cd "$SCRATCH/tools"
if [ ! -d puffin ]; then
    git clone https://github.com/jzhoulab/puffin.git
fi
cd puffin
echo "commit: $(git log -1 --format='%H %cd')"
echo "--- requirements.txt as published by the authors ---"
cat requirements.txt
echo "---------------------------------------------------"

echo "== 5. Puffin-D weights =="
mkdir -p resources
if [ ! -s resources/puffin_D.pth ]; then
    curl -fL -o resources/puffin_D.pth https://zenodo.org/record/7941729/files/puffin_D.pth
fi
ls -la resources/puffin_D.pth

echo "== 6. hg38 genome =="
# Downloaded fresh rather than copied from another project, so PURGE stays self-contained.
cd "$SCRATCH/ref"
if [ ! -s hg38.fa ]; then
    curl -fL -o hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gzip -d hg38.fa.gz
fi
ln -sf "$SCRATCH/ref/hg38.fa" "$SCRATCH/tools/puffin/resources/hg38.fa"
ls -la hg38.fa "$SCRATCH/tools/puffin/resources/hg38.fa"

echo "== 7. record what was installed =="
conda env export -p "$ENV" > "$PROJ/05_env/environment_purge_$STAMP.yml"
python -m pip freeze > "$PROJ/05_env/pip_freeze_$STAMP.txt"
(cd "$SCRATCH/tools/puffin" && git log -1 --format='%H') > "$PROJ/05_env/puffin_commit_$STAMP.txt"
md5sum "$SCRATCH/tools/puffin/resources/puffin_D.pth" "$SCRATCH/ref/hg38.fa" \
    > "$PROJ/05_env/MANIFEST_downloads_$STAMP.md5"
cat "$PROJ/05_env/MANIFEST_downloads_$STAMP.md5"

echo
echo "DONE. Environment: $ENV"
echo "Activate with: source $MF/etc/profile.d/conda.sh && conda activate $ENV"
