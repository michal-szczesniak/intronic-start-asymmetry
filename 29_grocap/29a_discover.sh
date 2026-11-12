#!/bin/bash
# Step 29a (version b: only file links named GSM..._ are kept; version a also took a link to an HHS web page)
# list and download the K562 GRO-cap files of GEO GSM1480321 (with TAP) and GSM1480322 (without TAP),
# and the UCSC hg19ToHg38 liftover chain; then check their format with 29a_inspect.py. No strand counts inside
# genes are computed here. Takes a few minutes; run it in the current srun session (needs internet access):
#   bash ${PURGE_SCRATCH:?}/purge29a/29a_discover.sh
set -e -o pipefail
P29A_HERE=$(cd "$(dirname "$0")" && pwd)
P29A_PROJ=${PURGE_PROJECT:?}
P29A_SCRATCH=${PURGE_SCRATCH:?}
P29A_DATA=$P29A_SCRATCH/29_grocap/data
P29A_DIR=$P29A_PROJ/29_grocap
export TMPDIR=$P29A_SCRATCH/.tmp XDG_CACHE_HOME=$P29A_SCRATCH/.cache MPLCONFIGDIR=$P29A_SCRATCH/.cache/matplotlib
mkdir -p "$TMPDIR" "$P29A_DATA" "$P29A_DIR/results" "$P29A_DIR/tests"
( cd "$P29A_HERE" && md5sum -c --quiet MD5SUMS ) || { echo "STOP: the bundle does not match its MD5SUMS"; exit 1; }
cp "$P29A_HERE"/29a_discover.sh "$P29A_HERE"/29a_inspect.py "$P29A_HERE"/README_29a.txt "$P29A_HERE"/MD5SUMS "$P29A_DIR"/
cp "$P29A_HERE"/tests/test_29a.py "$P29A_DIR"/tests/
source $P29A_SCRATCH/miniforge3/etc/profile.d/conda.sh
conda activate $P29A_SCRATCH/envs/purge
P29A_LOG=$P29A_DIR/results/29a_discover.log

p29a_get() {   # url, output file; up to 6 tries with growing waits
  local i
  for i in 1 2 3 4 5 6; do
    curl -fsSL "$1" -o "$2.part" && mv "$2.part" "$2" && return 0
    echo "  download failed (try $i): $1"; sleep $((30 * i))
  done
  echo "STOP: could not download $1"; return 1
}

{
  echo "== 29a $(date -u +%FT%TZ) on $(hostname)"
  echo "-- synthetic test of 29a_inspect.py"
  p29a_t=$(python "$P29A_DIR"/tests/test_29a.py "$P29A_DIR"/29a_inspect.py "$P29A_SCRATCH"/29_grocap/synth_test_29a | tail -1)
  echo "   $p29a_t"
  [[ "$p29a_t" == *", 0 differ" ]] || { echo "STOP: synthetic test of 29a_inspect.py differs"; exit 1; }
  cd "$P29A_DATA"
  for g in GSM1480321 GSM1480322; do
    u=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1480nnn/$g/suppl/
    p29a_get "$u" ${g}_suppl.html
    grep -o 'href="[^"]*"' ${g}_suppl.html | sed 's/^href="//; s/"$//' | grep -E "^${g}_[A-Za-z0-9._-]+$" > ${g}_files.txt || true
    echo "-- $g: $(wc -l < ${g}_files.txt) supplementary files"
    [ -s ${g}_files.txt ] || { echo "   no file links found; first lines of the listing:"; head -20 ${g}_suppl.html; }
    sed 's/^/   /' ${g}_files.txt
    while read -r f; do
      [ -s "$f" ] || p29a_get "$u$f" "$f"
    done < ${g}_files.txt
  done
  c=hg19ToHg38.over.chain.gz
  [ -s $c ] || p29a_get https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/$c $c
  echo "-- files (md5, bytes, name)"
  for f in $(cat GSM1480321_files.txt GSM1480322_files.txt) $c; do
    echo "   $(md5sum < "$f" | cut -d' ' -f1)  $(stat -c %s "$f")  $f"
  done
  python "$P29A_DIR"/29a_inspect.py $(cat GSM1480321_files.txt GSM1480322_files.txt | sed "s#^#$P29A_DATA/#") "$P29A_DATA/$c"
  echo "== done $(date -u +%FT%TZ)"
} 2>&1 | tee "$P29A_LOG"
echo "log: $P29A_LOG"
