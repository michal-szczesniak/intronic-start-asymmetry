#!/bin/bash
# Step 06: check that Puffin-D runs on a GPU here, and measure how long one prediction takes.
#
# Puffin-D takes a 100 kb sequence and returns ten tracks, in this order:
#   0 FANTOM_CAGE(+)  1 ENCODE_CAGE(+)  2 ENCODE_RAMPAGE(+)  3 GRO_CAP(+)  4 PRO_CAP(+)
#   5 PRO_CAP(-)      6 GRO_CAP(-)      7 ENCODE_RAMPAGE(-)  8 ENCODE_CAGE(-)  9 FANTOM_CAGE(-)
# PURGE needs tracks 4 and 5: PRO-cap on each strand.
#
# Run inside a GPU allocation, e.g.
#   srun --pty --partition=tesla --qos=tesla --gres=gpu:1 --mem=32G --time=04:00:00 /bin/bash
set -euo pipefail

SCRATCH=${PURGE_SCRATCH:?}
export TMPDIR=$SCRATCH/.tmp
export XDG_CACHE_HOME=$SCRATCH/.cache
export MPLCONFIGDIR=$SCRATCH/.cache/matplotlib
mkdir -p "$TMPDIR"

echo "== 1. allocation and GPU =="
echo "job: ${SLURM_JOB_ID:-none}  node: $(hostname)"
nvidia-smi || echo "nvidia-smi not available"

echo
echo "== 2. environment =="
source "$SCRATCH/miniforge3/etc/profile.d/conda.sh"
conda activate "$SCRATCH/envs/purge"
python - <<'PY'
import torch
print("torch", torch.__version__, "| cuda build", torch.version.cuda)
print("cuda available:", torch.cuda.is_available())
if torch.cuda.is_available():
    i = torch.cuda.current_device()
    print("device:", torch.cuda.get_device_name(i),
          "| capability", torch.cuda.get_device_capability(i),
          "| memory GB", round(torch.cuda.get_device_properties(i).total_memory / 1e9, 1))
PY

echo
echo "== 3. Puffin-D on one promoter (ACTB, the example from the authors' README) =="
cd "$SCRATCH/tools/puffin"
ls -la resources/
before=$(ls -1 | sort)
/usr/bin/time -v python puffin_D.py coord chr7:5530600- --use_cuda 2>&1 | tail -25 || true
echo "--- files created ---"
comm -13 <(echo "$before") <(ls -1 | sort) || true

echo
echo "== 4. timing: three predictions through the API =="
python - <<'PY'
import sys, time
sys.path.append(".")
import numpy as np
from puffin_D import *

t0 = time.time()
model = Puffin_D(use_cuda=True)
print(f"model load: {time.time() - t0:.1f} s")

from pyfaidx import Fasta
fa = Fasta("resources/hg38.fa")
# three 100 kb windows on chr7, including the ACTB locus
centres = [5530600, 5630600, 5730600]
for i, c in enumerate(centres):
    seq = str(fa["chr7"][c - 50000:c + 50000]).upper()
    t0 = time.time()
    pred = model.predict(seq)
    dt = time.time() - t0
    arr = np.asarray(pred)
    print(f"window {i + 1}: {dt:.2f} s | output shape {arr.shape} | "
          f"PRO-cap(+) max {arr[4].max():.3f} | PRO-cap(-) max {arr[5].max():.3f}")
PY

echo
echo "DONE"
