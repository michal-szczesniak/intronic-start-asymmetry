#!/bin/bash
# Step 27a (run on a second server, outside the cluster): sequence around the two starts of every depth-matched pair (PURGE step 23b) - base counts for
# logos and pre-specified core-promoter elements. Everything stays inside $PURGE27_DIR
# (the genome, about 1 GB gzipped, is downloaded there once). Run:  PURGE27_DIR=<folder> bash run_27a.sh
set -e -o pipefail
P27_DIR=${PURGE27_DIR:?}
cd $P27_DIR
mkdir -p ref results
echo "== 27a on $(hostname) $(date '+%F %T')"

# 1. bundle files unchanged
md5sum -c --quiet MD5SUMS || { echo "STOP: a bundle file differs from MD5SUMS"; exit 1; }
echo "bundle md5 OK"

# 2. python with numpy, pandas, scipy (installed into $P27_DIR/pylib only if missing)
P27_PY=$(command -v python3) || { echo "STOP: python3 not found"; exit 1; }
$P27_PY -c 'import sys; assert sys.version_info >= (3, 8), sys.version' || { echo "STOP: python3 >= 3.8 needed"; exit 1; }
if ! $P27_PY -c 'import numpy, pandas, scipy' 2>/dev/null; then
  echo "installing numpy, pandas, scipy into $P27_DIR/pylib"
  $P27_PY -m pip install --quiet --target $P27_DIR/pylib numpy pandas scipy || { echo "STOP: pip install failed"; exit 1; }
  export PYTHONPATH=$P27_DIR/pylib
fi
$P27_PY -c 'import numpy, pandas, scipy, sys; print("python", sys.version.split()[0], "numpy", numpy.__version__, "pandas", pandas.__version__, "scipy", scipy.__version__)'

# 3. hg38 from UCSC (once), md5 checked against UCSC's md5sum.txt
P27_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips
if [ ! -s ref/hg38.fa.gz.ok ]; then
  P27_FREE=$(df -Pk $P27_DIR | awk 'NR == 2 {print $4}')
  [ "$P27_FREE" -gt 2000000 ] || { echo "STOP: less than 2 GB free in $P27_DIR"; exit 1; }
  if command -v wget >/dev/null; then P27_GET="wget -q -O"; else P27_GET="curl -sSfL -o"; fi
  $P27_GET ref/md5sum.txt $P27_URL/md5sum.txt
  echo "downloading hg38.fa.gz (about 1 GB) ..."
  $P27_GET ref/hg38.fa.gz $P27_URL/hg38.fa.gz
  P27_WANT=$(awk '$2 == "hg38.fa.gz" {print $1}' ref/md5sum.txt)
  [ -n "$P27_WANT" ] || { echo "STOP: hg38.fa.gz not listed in md5sum.txt"; exit 1; }
  [ "$(md5sum < ref/hg38.fa.gz | cut -d' ' -f1)" = "$P27_WANT" ] || { echo "STOP: md5 of hg38.fa.gz differs from UCSC"; exit 1; }
  echo "$P27_WANT  hg38.fa.gz" > ref/hg38.fa.gz.ok
fi
echo "genome: $P27_DIR/ref/hg38.fa.gz, md5 $(cut -d' ' -f1 ref/hg38.fa.gz.ok) (UCSC)"

# 4. windows around the starts (one pass through the gzipped genome), then the analysis (10-30 min, low priority)
P27_ARMS="--arms K562=data/arms_K562.tsv --arms HCT116=data/arms_HCT116.tsv --matched data/23b_all_matched_pairs.tsv"
nice -n 19 $P27_PY 27a_start_motifs.py $P27_ARMS --genome ref/hg38.fa.gz --extract-windows results/27a_windows.fa.gz \
    2>&1 | tee results/27a_windows.log
nice -n 19 $P27_PY 27a_start_motifs.py $P27_ARMS --windows results/27a_windows.fa.gz --step10 10_shuffle_test.py \
    --out-prefix results/27a 2>&1 | tee results/27a.log
P27_OUT="results/27a_windows.log results/27a_windows.fa.gz results/27a.log results/27a_pfm.tsv results/27a_yr_offset.tsv results/27a_starts.tsv.gz results/27a_tests.tsv"
md5sum 27a_start_motifs.py ref/hg38.fa.gz.ok $P27_OUT > results/27a_outputs.md5
tar -czf results/27a_results.tar.gz $P27_OUT results/27a_outputs.md5
cat results/27a_outputs.md5
echo "== 27a done $(date '+%F %T'); please send: $P27_DIR/results/27a_results.tar.gz"
