#!/bin/bash
# Step 28d (EXPLORATORY): CpG islands and the strand difference within 5 kb of the gene's start.
# Needs the results of steps 28 and 28b. Submit from the login node:
#   sbatch --job-name=purge28d --mem=16G --time=08:00:00 --output=<dir>/run_28d_%j.log --wrap "bash <dir>/run_28d.sh"
set -e -o pipefail
[ -n "${SLURM_JOB_ID:-}" ] || { echo "STOP: run this on the compute node (sbatch or srun first)"; exit 1; }
P28D_HERE=$(cd "$(dirname "$0")" && pwd)
P28D_PROJ=${PURGE_PROJECT:?}
P28D_SCRATCH=${PURGE_SCRATCH:?}
P28D_DIR=$P28D_PROJ/28d_cpg
P28D_DATA=$P28D_SCRATCH/28d_cpg/data
P28D_R28=$P28D_PROJ/28_k36_position/results
P28D_R28B=$P28D_PROJ/28b_position/results
P28D_S17=$P28D_PROJ/17_observed/17_observed_sites.py
p28d_chk() { [ "$(md5sum < "$1" | cut -d' ' -f1)" = "$2" ] || { echo "STOP: md5 of $1 differs"; exit 1; }; }
export TMPDIR=$P28D_SCRATCH/.tmp XDG_CACHE_HOME=$P28D_SCRATCH/.cache MPLCONFIGDIR=$P28D_SCRATCH/.cache/matplotlib
mkdir -p "$TMPDIR" "$P28D_DATA" "$P28D_DIR/results"
source $P28D_SCRATCH/miniforge3/etc/profile.d/conda.sh
conda activate $P28D_SCRATCH/envs/purge

# 1. bundle and inputs
( cd "$P28D_HERE" && md5sum -c --quiet MD5SUMS ) || { echo "STOP: the bundle does not match its MD5SUMS"; exit 1; }
cp -r "$P28D_HERE"/28d_cpg.py "$P28D_HERE"/28d_download.py "$P28D_HERE"/run_28d.sh "$P28D_HERE"/README_28d.txt "$P28D_HERE"/MD5SUMS "$P28D_HERE"/tests "$P28D_DIR"/
p28d_chk $P28D_S17 659dca46e47fe9a9f35df48cef694889
p28d_chk $P28D_R28/28_K562_genes.tsv 2286575131e7e76cb06812f65fd336ef
p28d_chk $P28D_R28/28_HCT116_genes.tsv 0eb5f2cc15a26096d6355653760fbaae
p28d_chk $P28D_R28B/28b_K562_sites.tsv.gz ac95f13d0aff5be4b1939e12621575c0
p28d_chk $P28D_R28B/28b_HCT116_sites.tsv.gz 5fd83afc5e4c854fb3df008f4d7ec4db
p28d_chk $P28D_R28B/28b_K562_pairs.tsv 8019e1b26beec25a1657264b137d3997
p28d_chk $P28D_R28B/28b_HCT116_pairs.tsv 951d27c8f32c37c62314fce81c04feec
p28d_chk $P28D_R28B/28b_K562_strata.tsv 8ecf4540c3f3bd108ae910a2ea006bb7
p28d_chk $P28D_R28B/28b_HCT116_strata.tsv fbaed9ebac905eec7189d67a27ed2a11
echo "bundle, step-17 code and the step-28 and 28b tables: md5 OK"

# 2. synthetic test: the analysis on synthetic tables (tests/fixtures) and an independent loop-based check
cd "$P28D_DIR"
F=tests/fixtures
{
  python 28d_cpg.py --cgi $F/syn_cgi.txt.gz --r28b $F/r28b --r28 $F/r28 --step17 $P28D_S17 --min-islands 0 \
      --out-prefix "$P28D_SCRATCH"/28d_cpg/synth_test/28d &&
  python tests/naive_check_28d.py $F/syn_cgi.txt.gz $F/r28b $F/r28 "$P28D_SCRATCH"/28d_cpg/synth_test/28d K562 &&
  python tests/naive_check_28d.py $F/syn_cgi.txt.gz $F/r28b $F/r28 "$P28D_SCRATCH"/28d_cpg/synth_test/28d HCT116
} > results/28d_synth_test.log 2>&1 || { echo "STOP: synthetic test failed, see $P28D_DIR/results/28d_synth_test.log"; exit 1; }
[ "$(grep -c '^compared [0-9]* values, 0 differ' results/28d_synth_test.log)" = 2 ] || { echo "STOP: synthetic test: check differs"; exit 1; }
echo "synthetic test passed"

# 3. CpG islands (UCSC hg38 cpgIslandExt) and the run
python 28d_download.py "$P28D_DATA" results 2>&1 | tee results/28d_download.log
python 28d_cpg.py --cgi $P28D_DATA/cpgIslandExt.txt.gz --r28b $P28D_R28B --r28 $P28D_R28 --step17 $P28D_S17 \
    --out-prefix results/28d 2>&1 | tee results/28d.log

# 4. checksums and archive
( cd results && md5sum $(ls | grep -v '^28d_outputs.md5$') > 28d_outputs.md5 )
tar czf $P28D_SCRATCH/28d_results.tar.gz -C $P28D_PROJ 28d_cpg
md5sum $P28D_SCRATCH/28d_results.tar.gz
echo "done: $P28D_DIR/results; archive $P28D_SCRATCH/28d_results.tar.gz"
