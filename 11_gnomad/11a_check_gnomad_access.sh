#!/bin/bash
# Step 11a: can this node reach gnomAD v4.1 genome VCFs, how fast, and with which tools?
# Reads only HTTP headers and the first 50 MB of one file (to /dev/null); writes one report.
# Usage: bash 11a_check_gnomad_access.sh results/11a_check_<label>.txt
set -u
OUT=${1:?give the report path}
REL=release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr8.vcf.bgz
URLS="https://storage.googleapis.com/gcp-public-data--gnomad/$REL
https://gnomad-public-us-east-1.s3.amazonaws.com/$REL"
{
echo "== node $(hostname), job ${SLURM_JOB_ID:-none}, $(date '+%F %T')"
echo "== proxy settings"
env | grep -i -E '^(https?|no)_proxy=' || echo "none"
echo "== tools on PATH"
for t in bcftools tabix bgzip curl wget; do printf '%-9s ' "$t"; command -v "$t" || echo "not found"; done
command -v bcftools >/dev/null && bcftools --version | head -2
echo "== HTTP headers: status and size of the chr8 VCF and its index"
for U in $URLS; do
  echo "-- $U"
  curl -sSI --max-time 30 "$U" 2>&1 | grep -i -E '^HTTP|^content-length|^curl' || echo "no response"
  curl -sSI --max-time 30 "$U.tbi" 2>&1 | grep -i -E '^HTTP|^content-length|^curl' || echo "no response (.tbi)"
done
echo "== download speed: first 50 MB of each, discarded"
for U in $URLS; do
  printf '%s\n   ' "$U"
  curl -sS --max-time 120 -r 0-52428799 -o /dev/null -w 'http %{http_code}, %{size_download} bytes in %{time_total} s, %{speed_download} bytes/s\n' "$U" 2>&1
done
} 2>&1 | tee "$OUT"
