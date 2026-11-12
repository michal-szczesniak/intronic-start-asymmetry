#!/bin/bash
# Step 29b: K562 GRO-cap (GSM1480321, with TAP; hg19) lifted over to hg38, then the unchanged step-28 and step-28b
# scripts on it, and the summary of the three pre-specified questions (G1-G3) next to K562 PRO-cap.
# Needs step 29a (downloads in scratch/.../29_grocap/data). Submit from the login node:
#   sbatch --job-name=purge29b --mem=32G --time=08:00:00 --output=<dir>/run_29b_%j.log --wrap "bash <dir>/run_29b.sh"
set -e -o pipefail
[ -n "${SLURM_JOB_ID:-}" ] || { echo "STOP: run this on the compute node (sbatch)"; exit 1; }
if [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ "$SLURM_MEM_PER_NODE" -lt 31000 ]; then
  echo "STOP: this reservation has $SLURM_MEM_PER_NODE MB; at least 32G needed"; exit 1
fi
P29B_HERE=$(cd "$(dirname "$0")" && pwd)
P29B_PROJ=${PURGE_PROJECT:?}
P29B_SCRATCH=${PURGE_SCRATCH:?}
P29B_DIR=$P29B_PROJ/29_grocap
P29B_DATA=$P29B_SCRATCH/29_grocap/data
P29B_HG38=$P29B_SCRATCH/29_grocap/hg38
P29B_TEST=$P29B_SCRATCH/29_grocap/synth_test_29b
p29b_chk() { [ "$(md5sum < "$1" | cut -d' ' -f1)" = "$2" ] || { echo "STOP: md5 of $1 differs"; exit 1; }; }
export TMPDIR=$P29B_SCRATCH/.tmp XDG_CACHE_HOME=$P29B_SCRATCH/.cache MPLCONFIGDIR=$P29B_SCRATCH/.cache/matplotlib
mkdir -p "$TMPDIR" "$P29B_HG38" "$P29B_TEST" "$P29B_DIR/results" "$P29B_DIR/tests"
source $P29B_SCRATCH/miniforge3/etc/profile.d/conda.sh
conda activate $P29B_SCRATCH/envs/purge
echo "== 29b $(date -u +%FT%TZ) on $(hostname), job $SLURM_JOB_ID"

# 1. bundle and inputs
( cd "$P29B_HERE" && md5sum -c --quiet MD5SUMS ) || { echo "STOP: the bundle does not match its MD5SUMS"; exit 1; }
cp "$P29B_HERE"/29b_liftover.py "$P29B_HERE"/29b_summary.py "$P29B_HERE"/run_29b.sh "$P29B_HERE"/README_29b.txt "$P29B_DIR"/
cp "$P29B_HERE"/MD5SUMS "$P29B_DIR"/MD5SUMS_29b
cp "$P29B_HERE"/tests/test_29b_liftover.py "$P29B_DIR"/tests/
P29B_TP=$P29B_DATA/GSM1480321_K562_GROcap_wTAP_plus.bigWig
P29B_TM=$P29B_DATA/GSM1480321_K562_GROcap_wTAP_minus.bigWig
P29B_CHAIN=$P29B_DATA/hg19ToHg38.over.chain.gz
P29B_S10=$P29B_PROJ/10_shuffle/10_shuffle_test.py
P29B_S17=$P29B_PROJ/17_observed/17_observed_sites.py
P29B_S28=$P29B_PROJ/28_k36_position/28_k36_position.py
P29B_S28B=$P29B_PROJ/28b_position/28b_position.py
P29B_MP=$P29B_PROJ/23_robustness/results/23b/23b_all_matched_pairs.tsv
P29B_AK=$P29B_PROJ/02_istart_baseline/results/arms_K562.tsv
P29B_PK=$P29B_SCRATCH/28_k36_position/data/ENCFF520LHY.bed.gz
P29B_R28=$P29B_PROJ/28_k36_position/results/28_K562_strata.tsv
P29B_R28B=$P29B_PROJ/28b_position/results/28b_K562_strata.tsv
P29B_GENOME=$P29B_SCRATCH/ref/hg38.fa
P29B_GTF=$P29B_PROJ/02_istart_baseline/ref/gencode.v47.annotation.gtf.gz
for p in $P29B_TP $P29B_TM $P29B_CHAIN $P29B_S10 $P29B_S17 $P29B_S28 $P29B_S28B $P29B_MP $P29B_AK $P29B_PK $P29B_R28 \
         $P29B_R28B $P29B_GENOME $P29B_GENOME.fai $P29B_GTF; do
  [ -s "$p" ] || { echo "STOP: missing $p"; exit 1; }
done
p29b_chk $P29B_TP 40836b209f344e1de0a9247384b512f0
p29b_chk $P29B_TM ddec9080ce02d2e00a1604f13c1a8181
p29b_chk $P29B_CHAIN 35887f73fe5e2231656504d1f6430900
p29b_chk $P29B_S10 22669ea7a0e7f421e12cd9631081e101
p29b_chk $P29B_S17 659dca46e47fe9a9f35df48cef694889
p29b_chk $P29B_S28 c020e287552668cec7a74a786f2aa5b4
p29b_chk $P29B_S28B 08504653e2a0ea5742c10d0bd714461a
p29b_chk $P29B_MP 2d7dc179e597adb23f4338821a0f3cf4
p29b_chk $P29B_AK 6b3cf1ed296755f0fa858d4818e62de8
p29b_chk $P29B_PK 59c7b1385c96350fd2a7217f597ad495
p29b_chk $P29B_R28 e3b5e9eeb25d487a1d0b346eef5ca664
p29b_chk $P29B_R28B 8ecf4540c3f3bd108ae910a2ea006bb7
echo "bundle, GRO-cap files, chain, code of steps 10, 17, 28 and 28b, pairs, arms, peaks and PRO-cap tables: md5 OK"

# 2. synthetic test of the liftover against an independent position-by-position mapper
python "$P29B_DIR"/tests/test_29b_liftover.py "$P29B_DIR"/29b_liftover.py "$P29B_TEST" > "$P29B_DIR"/results/29b_synth_test.log 2>&1 \
  || { echo "STOP: synthetic test failed, see $P29B_DIR/results/29b_synth_test.log"; exit 1; }
grep -q "^compared [0-9]* values, 0 differ" "$P29B_DIR"/results/29b_synth_test.log || { echo "STOP: synthetic test: values differ"; exit 1; }
echo "synthetic test passed: $(grep '^compared' "$P29B_DIR"/results/29b_synth_test.log)"

# 3. liftover hg19 -> hg38
cd "$P29B_DIR"
python 29b_liftover.py --plus $P29B_TP --minus $P29B_TM --chain $P29B_CHAIN --fai $P29B_GENOME.fai \
    --out-prefix $P29B_HG38/29b_GROcap_hg38 2>&1 | tee results/29b_liftover.log
cp $P29B_HG38/29b_GROcap_hg38_liftover.tsv results/
md5sum $P29B_HG38/29b_GROcap_hg38_plus.bw $P29B_HG38/29b_GROcap_hg38_minus.bw | tee -a results/29b_liftover.log

# 4. step 28 unchanged on GRO-cap (primary: B = G1, A2 = G2)
P29B_CHR=$(echo chr{1..22} | tr ' ' ,)
python $P29B_S28 --line K562_GROcap --plus-bw $P29B_HG38/29b_GROcap_hg38_plus.bw --minus-bw $P29B_HG38/29b_GROcap_hg38_minus.bw \
    --peaks $P29B_PK --primary A2,B --chroms $P29B_CHR --genome $P29B_GENOME --gtf $P29B_GTF --step10 $P29B_S10 \
    --step17 $P29B_S17 --out-prefix results/29_K562_GROcap 2>&1 | tee results/29_K562_GROcap.log
P29B_T=$(sed -n 's/.*strong sites: >= \([0-9]*\) reads.*/\1/p' results/29_K562_GROcap.log | head -1)
[ -n "$P29B_T" ] || { echo "STOP: strong-site threshold not found in the step-28 log"; exit 1; }
echo "strong sites in GRO-cap: >= $P29B_T reads"

# 5. step 28b unchanged on GRO-cap (sites by distance; G3 = 1-5-rest). --line K562 selects the K562 pairs, whose
#    bias must reproduce step 23b (the pairs come from PRO-cap and do not change).
python $P29B_S28B --line K562 --plus-bw $P29B_HG38/29b_GROcap_hg38_plus.bw --minus-bw $P29B_HG38/29b_GROcap_hg38_minus.bw \
    --strong $P29B_T --expect-bias 33.0 --genes28 results/29_K562_GROcap_genes.tsv --matched $P29B_MP --arms $P29B_AK \
    --chroms $P29B_CHR --genome $P29B_GENOME --gtf $P29B_GTF --step10 $P29B_S10 --step17 $P29B_S17 --step28 $P29B_S28 \
    --out-prefix results/29_GROcap_28b 2>&1 | tee results/29_GROcap_28b.log

# 6. summary of G1-G3 next to PRO-cap, with the pre-specified reading
python 29b_summary.py results/29_K562_GROcap_strata.tsv results/29_GROcap_28b_strata.tsv $P29B_R28 $P29B_R28B \
    results/29b_summary 2>&1 | tee results/29b_summary.log

# 7. checksums and archive (code and results; the lifted bigWigs stay in scratch, their md5 is in the log)
( cd results && md5sum $(ls | grep -v '^29b_outputs.md5$') > 29b_outputs.md5 )
tar czf $P29B_SCRATCH/29b_results.tar.gz -C $P29B_PROJ 29_grocap
md5sum $P29B_SCRATCH/29b_results.tar.gz
echo "done $(date -u +%FT%TZ): $P29B_DIR/results; archive $P29B_SCRATCH/29b_results.tar.gz"
