#!/bin/bash
# Step 28: H3K36me3, position along the gene and HCT116 for the gene-strand deficit of strong measured start sites.
# Run from the unpacked bundle, inside an srun reservation with at least 32G:   bash run_28.sh
# Writes code and results to project_data/.../PURGE/28_k36_position, downloads and the synthetic test to scratch.
set -e -o pipefail
[ -n "${SLURM_JOB_ID:-}" ] || { echo "STOP: run this on the compute node (srun first)"; exit 1; }
if [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ "$SLURM_MEM_PER_NODE" -lt 31000 ]; then
  echo "STOP: this reservation has $SLURM_MEM_PER_NODE MB; at least 32G needed (exit, then srun with --mem=32G)"; exit 1
fi
P28_HERE=$(cd "$(dirname "$0")" && pwd)
P28_PROJ=${PURGE_PROJECT:?}
P28_SCRATCH=${PURGE_SCRATCH:?}
P28_DIR=$P28_PROJ/28_k36_position
P28_DATA=$P28_SCRATCH/28_k36_position/data
P28_TEST=$P28_SCRATCH/28_k36_position/synth_test
p28_chk() { [ "$(md5sum < "$1" | cut -d' ' -f1)" = "$2" ] || { echo "STOP: md5 of $1 differs"; exit 1; }; }
export TMPDIR=$P28_SCRATCH/.tmp XDG_CACHE_HOME=$P28_SCRATCH/.cache MPLCONFIGDIR=$P28_SCRATCH/.cache/matplotlib
mkdir -p "$TMPDIR" "$P28_DATA" "$P28_TEST" "$P28_DIR/results" "$P28_DIR/tests"
source $P28_SCRATCH/miniforge3/etc/profile.d/conda.sh
conda activate $P28_SCRATCH/envs/purge
python -c "import pyBigWig, pyfaidx, scipy, pandas, numpy" || { echo "STOP: pyBigWig, pyfaidx, scipy, pandas or numpy missing in the purge environment"; exit 1; }

# 1. bundle and inputs from earlier steps
( cd "$P28_HERE" && md5sum -c --quiet MD5SUMS ) || { echo "STOP: the bundle does not match its MD5SUMS"; exit 1; }
cp "$P28_HERE"/28_k36_position.py "$P28_HERE"/28_download.py "$P28_HERE"/run_28.sh "$P28_HERE"/README_28.txt "$P28_HERE"/MD5SUMS "$P28_DIR"/
cp "$P28_HERE"/tests/make_synth.py "$P28_HERE"/tests/naive_check.py "$P28_DIR"/tests/
P28_S10=$P28_PROJ/10_shuffle/10_shuffle_test.py
P28_S17=$P28_PROJ/17_observed/17_observed_sites.py
P28_G17=$P28_PROJ/17_observed/results/17_obs_genes.tsv
P28_GENOME=$P28_SCRATCH/ref/hg38.fa
P28_GTF=$P28_PROJ/02_istart_baseline/ref/gencode.v47.annotation.gtf.gz
P28_KP=$P28_SCRATCH/models/procapnet_K562/ENCFF798GNW.bigWig
P28_KM=$P28_SCRATCH/models/procapnet_K562/ENCFF662SHP.bigWig
for p in $P28_S10 $P28_S17 $P28_G17 $P28_GENOME $P28_GENOME.fai $P28_GTF $P28_KP $P28_KM; do
  [ -s "$p" ] || { echo "STOP: missing $p"; exit 1; }
done
p28_chk $P28_S10 22669ea7a0e7f421e12cd9631081e101
p28_chk $P28_S17 659dca46e47fe9a9f35df48cef694889
p28_chk $P28_G17 f84a1410f272072c825e604ce7e41dd8
p28_chk $P28_KP b36ad5972b742ade1ba41e883a32ff61
p28_chk $P28_KM d1832c22ed40cf9e6f21acc10a81af75
echo "bundle, step-10 and step-17 code, step-17 table and K562 bigWigs: md5 OK"

# 2. synthetic test: step 17 and step 28 on the same synthetic data, and an independent loop-based recount
cd "$P28_TEST"
cp "$P28_HERE"/tests/make_synth.py "$P28_HERE"/tests/naive_check.py .
cp $P28_S10 10_shuffle_test.py
cp $P28_S17 17_observed_sites.py
{
  python make_synth.py
  python 17_observed_sites.py --plus-bw plus_all.bw --minus-bw minus_all.bw --chroms chr1,chr2 --genome synth.fa \
      --gtf synth.gtf.gz --step10 10_shuffle_test.py --out-prefix out17/17_obs --n-boot 200 > out17.log
  python "$P28_DIR"/28_k36_position.py --line SYN --plus-bw plus_rep1.bw,plus_rep2.bw --minus-bw minus_rep1.bw,minus_rep2.bw \
      --peaks peaks.bed.gz --primary A1,A2,B --genes17 out17/17_obs_genes.tsv --chroms chr1,chr2 --genome synth.fa \
      --gtf synth.gtf.gz --step10 10_shuffle_test.py --step17 17_observed_sites.py --out-prefix out28/28_SYN \
      --n-boot 500 --ref-reads 5000
  python naive_check.py
} > "$P28_DIR"/results/28_synth_test.log 2>&1 || { echo "STOP: synthetic test failed, see $P28_DIR/results/28_synth_test.log"; exit 1; }
grep -q "step 17 reproduced exactly" "$P28_DIR"/results/28_synth_test.log || { echo "STOP: synthetic test: step 17 not reproduced"; exit 1; }
grep -q "^compared [0-9]* counts, 0 differ" "$P28_DIR"/results/28_synth_test.log || { echo "STOP: synthetic test: recount differs"; exit 1; }
echo "synthetic test passed: $(grep '^compared' "$P28_DIR"/results/28_synth_test.log)"

# 3. ENCODE downloads (md5 against the metadata of 10 Nov 2026) and choice of the primary HCT116 bigWig set
cd "$P28_DIR"
python 28_download.py "$P28_DATA" results 2>&1 | tee results/28_download.log
P28_HP=$(awk -F'\t' '$1=="primary"{print $3}' results/28_hct116_sets.tsv | sed "s#\([^,]*\)#$P28_DATA/\1.bigWig#g")
P28_HM=$(awk -F'\t' '$1=="primary"{print $4}' results/28_hct116_sets.tsv | sed "s#\([^,]*\)#$P28_DATA/\1.bigWig#g")
P28_AP=$(awk -F'\t' '$1=="sensitivity"{print $3}' results/28_hct116_sets.tsv | sed "s#\([^,]*\)#$P28_DATA/\1.bigWig#g")
P28_AM=$(awk -F'\t' '$1=="sensitivity"{print $4}' results/28_hct116_sets.tsv | sed "s#\([^,]*\)#$P28_DATA/\1.bigWig#g")
[ -n "$P28_HP" ] && [ -n "$P28_HM" ] && [ -n "$P28_AP" ] && [ -n "$P28_AM" ] || { echo "STOP: HCT116 sets not resolved"; exit 1; }

# 4. the three runs
P28_CHR=$(echo chr{1..22} | tr ' ' ,)
P28_ARGS="--chroms $P28_CHR --genome $P28_GENOME --gtf $P28_GTF --step10 $P28_S10 --step17 $P28_S17"
python 28_k36_position.py --line K562 --plus-bw $P28_KP --minus-bw $P28_KM --peaks $P28_DATA/ENCFF520LHY.bed.gz \
    --primary A1,A2 --genes17 $P28_G17 $P28_ARGS --out-prefix results/28_K562 2>&1 | tee results/28_K562.log
python 28_k36_position.py --line HCT116 --plus-bw $P28_HP --minus-bw $P28_HM --peaks $P28_DATA/ENCFF029GQD.bed.gz \
    --primary A1,A2,B $P28_ARGS --out-prefix results/28_HCT116 2>&1 | tee results/28_HCT116.log
python 28_k36_position.py --line HCT116_sensitivity --plus-bw $P28_AP --minus-bw $P28_AM --peaks $P28_DATA/ENCFF029GQD.bed.gz \
    --primary "" $P28_ARGS --out-prefix results/28_HCT116_sensitivity 2>&1 | tee results/28_HCT116_sensitivity.log

# 5. checksums and an archive for transfer (scratch)
( cd results && md5sum $(ls | grep -v '^28_outputs.md5$') > 28_outputs.md5 )
tar czf $P28_SCRATCH/28_results.tar.gz -C $P28_PROJ 28_k36_position
md5sum $P28_SCRATCH/28_results.tar.gz
echo "done: $P28_DIR/results; archive $P28_SCRATCH/28_results.tar.gz"
