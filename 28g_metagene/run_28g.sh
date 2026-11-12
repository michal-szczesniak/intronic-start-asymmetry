#!/bin/bash
# Step 28g: strong measured sites in 1-kb bins from the 5' end of active genes (Figure 3F), from the site tables of
# steps 28b and 29b. Seconds of CPU; run on a compute node, e.g. inside an existing srun reservation:
#   bash ${PURGE_PROJECT:?}/28g_metagene/run_28g.sh
set -e -o pipefail
[ -n "${SLURM_JOB_ID:-}" ] || { echo "STOP: run this on the compute node (sbatch or srun first)"; exit 1; }
P28G_PROJ=${PURGE_PROJECT:?}
P28G_SCRATCH=${PURGE_SCRATCH:?}
P28G_DIR=$P28G_PROJ/28g_metagene
P28G_OUT=$P28G_SCRATCH/28g_metagene/rerun
p28g_chk() { [ "$(md5sum < "$1" | cut -d' ' -f1)" = "$2" ] || { echo "STOP: md5 of $1 differs"; exit 1; }; }
export TMPDIR=$P28G_SCRATCH/.tmp XDG_CACHE_HOME=$P28G_SCRATCH/.cache
mkdir -p "$TMPDIR" "$P28G_OUT"
source $P28G_SCRATCH/miniforge3/etc/profile.d/conda.sh
conda activate $P28G_SCRATCH/envs/purge
( cd "$P28G_DIR" && md5sum -c --quiet MD5SUMS ) || { echo "STOP: the bundle does not match its MD5SUMS"; exit 1; }
P28G_SK=$P28G_PROJ/28b_position/results/28b_K562_sites.tsv.gz
P28G_SH=$P28G_PROJ/28b_position/results/28b_HCT116_sites.tsv.gz
P28G_SG=$P28G_PROJ/29_grocap/results/29_GROcap_28b_sites.tsv.gz
P28G_GK=$P28G_PROJ/28_k36_position/results/28_K562_genes.tsv
P28G_GH=$P28G_PROJ/28_k36_position/results/28_HCT116_genes.tsv
P28G_GG=$P28G_PROJ/29_grocap/results/29_K562_GROcap_genes.tsv
p28g_chk $P28G_SK ac95f13d0aff5be4b1939e12621575c0
p28g_chk $P28G_SH 5fd83afc5e4c854fb3df008f4d7ec4db
p28g_chk $P28G_SG ddfb3f3c496cef18925ed9f9bdd52056
p28g_chk $P28G_GK 2286575131e7e76cb06812f65fd336ef
p28g_chk $P28G_GH 0eb5f2cc15a26096d6355653760fbaae
p28g_chk $P28G_GG 943264f2710ba4cad030e081466c269f
python "$P28G_DIR"/tests/test_28g.py "$P28G_DIR"/28g_metagene.py "$P28G_OUT"/synth > "$P28G_OUT"/28g_synth_test.log 2>&1
grep -q "^compared [0-9]* values, 0 differ" "$P28G_OUT"/28g_synth_test.log || { echo "STOP: synthetic test failed"; exit 1; }
cd "$P28G_OUT"
python "$P28G_DIR"/28g_metagene.py --label "K562 PRO-cap" --sites $P28G_SK --genes $P28G_GK --expect-genes 10087 \
    --expect-total 1224,1483 --expect-first 232,431 --out 28g_K562_bins.tsv > 28g_K562.log
python "$P28G_DIR"/28g_metagene.py --label "HCT116 PRO-cap" --sites $P28G_SH --genes $P28G_GH --expect-genes 10538 \
    --expect-total 927,1092 --expect-first 111,254 --out 28g_HCT116_bins.tsv > 28g_HCT116.log
python "$P28G_DIR"/28g_metagene.py --label "K562 GRO-cap" --sites $P28G_SG --genes $P28G_GG --expect-genes 9907 \
    --expect-total 2450,2771 --expect-first 412,667 --out 28g_GROcap_bins.tsv > 28g_GROcap.log
for f in 28g_K562_bins.tsv 28g_HCT116_bins.tsv 28g_GROcap_bins.tsv; do
  cmp -s $f "$P28G_DIR"/results/$f && echo "$f: identical to the delivered table" || echo "$f: DIFFERS (versions of pandas or numpy?)"
done
