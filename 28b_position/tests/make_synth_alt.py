"""Copy synth.gtf.gz to synth2.gtf.gz, adding an alternative transcript that starts in the middle of every third
gene (to test distances from the nearest upstream TSS in step 28b)."""
import gzip, re
out = gzip.open("synth2.gtf.gz", "wt")
k = 0
for line in gzip.open("synth.gtf.gz", "rt"):
    out.write(line)
    f = line.split("\t")
    if f[2] == "gene":
        k += 1
        if k % 3 == 0:
            s, e = int(f[3]), int(f[4])
            m = (s + e) // 2
            gid = re.search('gene_id "([^"]+)"', f[8]).group(1)
            a = f[8].rstrip("\n") + f' transcript_id "ALT{gid}";'
            if f[6] == "+":
                out.write(f"{f[0]}\tS\ttranscript\t{m}\t{e}\t.\t+\t.\t{a}\n{f[0]}\tS\texon\t{m}\t{m + 200}\t.\t+\t.\t{a}\n")
            else:
                out.write(f"{f[0]}\tS\ttranscript\t{s}\t{m}\t.\t-\t.\t{a}\n{f[0]}\tS\texon\t{m - 200}\t{m}\t.\t-\t.\t{a}\n")
out.close()
print("written synth2.gtf.gz")
