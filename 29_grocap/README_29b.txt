Step 29b, PURGE (11 Nov 2026): replication of the strong-site deficit with K562 GRO-cap (pre-specified, step 29).

Input (downloaded and checked in step 29a): GEO GSM1480321, K562 GRO-cap with TAP, hg19 single-base counts
(plus file positive, minus file negative values), and the UCSC chain hg19ToHg38.over.chain.gz.

What run_29b.sh does
  1. checks the bundle, the GRO-cap files, the chain and every input from earlier steps by md5;
  2. synthetic test of the liftover against an independent position-by-position mapper
     (gaps, an inverted chain, overlapping chains, a sequence that is not kept, added reads);
  3. lifts the counts to hg38 (29b_liftover.py; positions in inverted blocks change strand; unmapped and
     ambiguous positions are dropped and counted; the lifted bigWigs stay in scratch);
  4. runs the step-28 script unchanged (G1 = deficit in active genes, G2 = 3' fifth minus 5' fifth);
  5. runs the step-28b script unchanged (G3 = within 5 kb of the 5' end minus beyond; it also re-checks the
     K562 pair bias of step 23b, which does not depend on GRO-cap);
  6. 29b_summary.py: G1-G3 next to K562 PRO-cap and the pre-specified reading.
  7. checksums and archive scratch/.../PURGE/29b_results.tar.gz.

Files: 29b_liftover.py, 29b_summary.py, run_29b.sh, tests/test_29b_liftover.py, MD5SUMS
Submit from the login node (about 30-60 min, 32 GB):
  sbatch --job-name=purge29b --mem=32G --time=08:00:00 \
    --output=<scratch folder>/purge29b/run_29b_%j.log \
    --wrap "bash <scratch folder>/purge29b/run_29b.sh"
