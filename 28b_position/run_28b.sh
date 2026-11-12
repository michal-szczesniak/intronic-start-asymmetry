#!/bin/bash
# Step 28b (EXPLORATORY): position along the host gene of the pair bias (Results 1) and of the strong-site deficit
# (Results 4). Needs the results of step 28. Submit from the login node:
#   sbatch --job-name=purge28b --mem=32G --time=08:00:00 --output=<dir>/run_28b_%j.log --wrap "bash <dir>/run_28b.sh"
set -e -o pipefail
[ -n "${SLURM_JOB_ID:-}" ] || { echo "STOP: run this on the compute node (sbatch or srun first)"; exit 1; }
if [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ "$SLURM_MEM_PER_NODE" -lt 31000 ]; then
  echo "STOP: this reservation has $SLURM_MEM_PER_NODE MB; at least 32G needed"; exit 1
fi
P28B_HERE=$(cd "$(dirname "$0")" && pwd)
P28B_PROJ=${PURGE_PROJECT:?}
P28B_SCRATCH=${PURGE_SCRATCH:?}
P28B_DIR=$P28B_PROJ/28b_position
P28B_TEST=$P28B_SCRATCH/28b_position/synth_test
P28B_DATA=$P28B_SCRATCH/28_k36_position/data
p28b_chk() { [ "$(md5sum < "$1" | cut -d' ' -f1)" = "$2" ] || { echo "STOP: md5 of $1 differs"; exit 1; }; }
export TMPDIR=$P28B_SCRATCH/.tmp XDG_CACHE_HOME=$P28B_SCRATCH/.cache MPLCONFIGDIR=$P28B_SCRATCH/.cache/matplotlib
mkdir -p "$TMPDIR" "$P28B_TEST" "$P28B_DIR/results" "$P28B_DIR/tests"
source $P28B_SCRATCH/miniforge3/etc/profile.d/conda.sh
conda activate $P28B_SCRATCH/envs/purge

# 1. bundle and inputs
( cd "$P28B_HERE" && md5sum -c --quiet MD5SUMS ) || { echo "STOP: the bundle does not match its MD5SUMS"; exit 1; }
cp "$P28B_HERE"/28b_position.py "$P28B_HERE"/run_28b.sh "$P28B_HERE"/README_28b.txt "$P28B_HERE"/MD5SUMS "$P28B_DIR"/
cp "$P28B_HERE"/tests/*.py "$P28B_DIR"/tests/
P28B_S10=$P28B_PROJ/10_shuffle/10_shuffle_test.py
P28B_S17=$P28B_PROJ/17_observed/17_observed_sites.py
P28B_S28=$P28B_PROJ/28_k36_position/28_k36_position.py
P28B_SYN=$P28B_PROJ/28_k36_position/tests/make_synth.py
P28B_GK=$P28B_PROJ/28_k36_position/results/28_K562_genes.tsv
P28B_GH=$P28B_PROJ/28_k36_position/results/28_HCT116_genes.tsv
P28B_MP=$P28B_PROJ/23_robustness/results/23b/23b_all_matched_pairs.tsv
P28B_AK=$P28B_PROJ/02_istart_baseline/results/arms_K562.tsv
P28B_AH=$P28B_PROJ/02_istart_baseline/results/arms_HCT116.tsv
P28B_GENOME=$P28B_SCRATCH/ref/hg38.fa
P28B_GTF=$P28B_PROJ/02_istart_baseline/ref/gencode.v47.annotation.gtf.gz
P28B_KP=$P28B_SCRATCH/models/procapnet_K562/ENCFF798GNW.bigWig
P28B_KM=$P28B_SCRATCH/models/procapnet_K562/ENCFF662SHP.bigWig
P28B_HP=$P28B_DATA/ENCFF440UVY.bigWig,$P28B_DATA/ENCFF078CJH.bigWig
P28B_HM=$P28B_DATA/ENCFF103FBH.bigWig,$P28B_DATA/ENCFF907PJM.bigWig
for p in $P28B_S10 $P28B_S17 $P28B_S28 $P28B_SYN $P28B_GK $P28B_GH $P28B_MP $P28B_AK $P28B_AH $P28B_GENOME $P28B_GENOME.fai $P28B_GTF $P28B_KP $P28B_KM; do
  [ -s "$p" ] || { echo "STOP: missing $p"; exit 1; }
done
p28b_chk $P28B_S10 22669ea7a0e7f421e12cd9631081e101
p28b_chk $P28B_S17 659dca46e47fe9a9f35df48cef694889
p28b_chk $P28B_S28 c020e287552668cec7a74a786f2aa5b4
p28b_chk $P28B_SYN 620b1f78a25a986e0d8916f9bd703669
p28b_chk $P28B_GK 2286575131e7e76cb06812f65fd336ef
p28b_chk $P28B_GH 0eb5f2cc15a26096d6355653760fbaae
p28b_chk $P28B_MP 2d7dc179e597adb23f4338821a0f3cf4
p28b_chk $P28B_AK 6b3cf1ed296755f0fa858d4818e62de8
p28b_chk $P28B_AH affec8bc6128e0f37090a90e12523f90
p28b_chk $P28B_KP b36ad5972b742ade1ba41e883a32ff61
p28b_chk $P28B_KM d1832c22ed40cf9e6f21acc10a81af75
p28b_chk $P28B_DATA/ENCFF440UVY.bigWig 9d22de120814a72ce719d8a3274f3694
p28b_chk $P28B_DATA/ENCFF078CJH.bigWig 7aa4cfb8fcaefb8b55f7e860156c83f3
p28b_chk $P28B_DATA/ENCFF103FBH.bigWig 847a7cb26342b4f68bc7d92a20ae22c2
p28b_chk $P28B_DATA/ENCFF907PJM.bigWig b10d937242984712bbb32ca0ba59540f
echo "bundle, code of steps 10, 17 and 28, step-28 gene tables, 23b pairs, arm tables and bigWigs: md5 OK"

# 2. synthetic test: synthetic genome with alternative TSSs, step 28 on it, step 28b, independent recount
cd "$P28B_TEST"
cp $P28B_SYN make_synth.py
cp "$P28B_HERE"/tests/*.py .
cp $P28B_S10 10_shuffle_test.py
cp $P28B_S17 17_observed_sites.py
cp $P28B_S28 28_k36_position.py
{
  python make_synth.py &&
  python make_synth_alt.py &&
  python 28_k36_position.py --line SYN --plus-bw plus_all.bw --minus-bw minus_all.bw --peaks peaks.bed.gz --chroms chr1,chr2 \
      --genome synth.fa --gtf synth2.gtf.gz --step10 10_shuffle_test.py --step17 17_observed_sites.py \
      --out-prefix out28/28_SYN --n-boot 100 > out28.log &&
  P28B_T=$(sed -n 's/.*strong sites: >= \([0-9]*\) reads.*/\1/p' out28.log | head -1) &&
  P28B_EXP=$(python make_pairs.py out28/28_SYN_genes.tsv pairs_id.tsv arms_id.tsv | awk '{print $NF}') &&
  echo "synthetic: strong >= $P28B_T reads; expected overall bias $P28B_EXP" &&
  python "$P28B_DIR"/28b_position.py --line SYN --plus-bw plus_all.bw --minus-bw minus_all.bw --strong $P28B_T \
      --expect-bias $P28B_EXP --genes28 out28/28_SYN_genes.tsv --matched pairs_id.tsv --arms arms_id.tsv --chroms chr1,chr2 --genome synth.fa \
      --gtf synth2.gtf.gz --step10 10_shuffle_test.py --step17 17_observed_sites.py --step28 28_k36_position.py \
      --out-prefix out28b/28b_SYN --n-boot 200 &&
  python naive_check_28b.py synth2.gtf.gz out28/28_SYN_genes.tsv pairs_id.tsv arms_id.tsv out28b/28b_SYN $P28B_T
} > "$P28B_DIR"/results/28b_synth_test.log 2>&1 || { echo "STOP: synthetic test failed, see $P28B_DIR/results/28b_synth_test.log"; exit 1; }
grep -q "^compared [0-9]* values, 0 differ" "$P28B_DIR"/results/28b_synth_test.log || { echo "STOP: synthetic test: recount differs"; exit 1; }
echo "synthetic test passed: $(grep '^compared' "$P28B_DIR"/results/28b_synth_test.log)"

# 3. the two runs (strong-site thresholds and step-23b biases from steps 28 and 23b)
cd "$P28B_DIR"
P28B_CHR=$(echo chr{1..22} | tr ' ' ,)
P28B_ARGS="--matched $P28B_MP --chroms $P28B_CHR --genome $P28B_GENOME --gtf $P28B_GTF --step10 $P28B_S10 --step17 $P28B_S17 --step28 $P28B_S28"
python 28b_position.py --line K562 --plus-bw $P28B_KP --minus-bw $P28B_KM --strong 10 --expect-bias 33.0 --genes28 $P28B_GK --arms $P28B_AK \
    $P28B_ARGS --out-prefix results/28b_K562 2>&1 | tee results/28b_K562.log
python 28b_position.py --line HCT116 --plus-bw $P28B_HP --minus-bw $P28B_HM --strong 24 --expect-bias 20.6 --genes28 $P28B_GH --arms $P28B_AH \
    $P28B_ARGS --out-prefix results/28b_HCT116 2>&1 | tee results/28b_HCT116.log

# 4. checksums and archive
( cd results && md5sum $(ls | grep -v '^28b_outputs.md5$') > 28b_outputs.md5 )
tar czf $P28B_SCRATCH/28b_results.tar.gz -C $P28B_PROJ 28b_position
md5sum $P28B_SCRATCH/28b_results.tar.gz
echo "done: $P28B_DIR/results; archive $P28B_SCRATCH/28b_results.tar.gz"
