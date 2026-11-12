#!/usr/bin/env python3
"""Step 28d: download the UCSC CpG-island table for hg38 (cpgIslandExt.txt.gz) with retries after server errors,
and record its size, md5 and date. UCSC publishes no checksum for database tables, so the md5 recorded here
identifies the copy used. Usage: python3 28d_download.py <data_dir> <results_dir>"""
import hashlib, os, sys, time, urllib.error, urllib.request
URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"
WAITS = (30, 60, 120, 240, 480)
data, res = sys.argv[1], sys.argv[2]
os.makedirs(data, exist_ok=True)
os.makedirs(res, exist_ok=True)
dst = os.path.join(data, "cpgIslandExt.txt.gz")
if os.path.isfile(dst) and os.path.getsize(dst) > 0:
    print(f"present: {dst}")
    lastmod = "previously downloaded"
else:
    for k in range(len(WAITS) + 1):
        try:
            with urllib.request.urlopen(URL, timeout=600) as src, open(dst + ".part", "wb") as out:
                lastmod = src.headers.get("Last-Modified", "")
                while True:
                    b = src.read(1 << 22)
                    if not b:
                        break
                    out.write(b)
            os.replace(dst + ".part", dst)
            break
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, ConnectionError) as ex:
            if (isinstance(ex, urllib.error.HTTPError) and ex.code < 500) or k == len(WAITS):
                raise
            print(f"  {ex}; new attempt in {WAITS[k]} s", flush=True)
            time.sleep(WAITS[k])
    print(f"downloaded: {URL} (Last-Modified: {lastmod})")
h = hashlib.md5(open(dst, "rb").read()).hexdigest()
with open(os.path.join(res, "28d_download.tsv"), "w") as fh:
    fh.write("url\tpath\tbytes\tmd5\tlast_modified\tchecked_on\n")
    fh.write(f"{URL}\t{dst}\t{os.path.getsize(dst)}\t{h}\t{lastmod}\t{time.strftime('%Y-%m-%d')}\n")
print(f"md5 {h}, {os.path.getsize(dst):,} bytes")
