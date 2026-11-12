Step 29a, PURGE (11 Nov 2026): GRO-cap files for the replication of the strong-site deficit (step 29).
The questions of step 29 were recorded in the dated project log before this step (pre-specification, step 29).

What it does
  1. a synthetic test of 29a_inspect.py (known totals must be printed back; "compared 17 values, 0 differ");
  2. lists and downloads the supplementary files of GEO GSM1480321 (K562 GRO-cap with TAP) and GSM1480322
     (K562 GRO-cap without TAP), series GSE60456 (Core et al. 2014, Nat Genet), and the UCSC chain
     hg19ToHg38.over.chain.gz, into scratch/.../PURGE/29_grocap/data;
  3. prints md5 and size of every file and checks the format: sequences, intervals, whole-number counts,
     negative values, totals per chromosome group and the signal near two known starts (GAPDH +, ACTB -).
No strand counts inside genes are computed; the analysis (step 29b) is written after this check.

Files
  29a_discover.sh   run script (bash, in an srun session with internet access)
  29a_inspect.py    format check
  tests/test_29a.py synthetic test
  MD5SUMS

Run
  bash <scratch folder>/purge29a/29a_discover.sh
Log: 29_grocap/results/29a_discover.log
