#!/usr/bin/env python3
"""Step 28: download the ENCODE files (H3K36me3 peaks of both lines, HCT116 PRO-cap bigWigs), check them against
the ENCODE metadata recorded on 10 Nov 2026, and choose the primary HCT116 bigWig set.

HCT116 PRO-cap (ENCSR165QKR) has two processings of the same two replicates, each with plus- and minus-strand
bigWigs of unique reads:
  set 266IRD  plus ENCFF266IRD (rep 1), ENCFF209TJK (rep 2); minus ENCFF116OZK (rep 1), ENCFF828GNJ (rep 2)
  set 440UVY  plus ENCFF440UVY (rep 1), ENCFF078CJH (rep 2); minus ENCFF103FBH (rep 1), ENCFF907PJM (rep 2)
Rule fixed before any data were seen: the primary set is the one in the experiment's ENCODE default analysis; if
neither or both sets are, the set processed like the K562 data (outputs: bidirectional and unidirectional peaks, no
divergent peaks), i.e. set 266IRD. The other set is a sensitivity analysis.

ENCODE requests are retried up to 6 times (waits 30 s, 1, 2, 4 and 8 min) when the portal answers with a server
error (5xx, e.g. 504 Gateway Time-out) or the connection fails; a 4xx answer stops at once.

Usage: python3 28_download.py <data_dir> <results_dir>
Writes <results_dir>/28_downloads.tsv and <results_dir>/28_hct116_sets.tsv; prints a log.
"""
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request

BASE = "https://www.encodeproject.org"
# accession: (extension, md5 listed by ENCODE on 10 Nov 2026, expected output type)
FILES = {
    "ENCFF520LHY": ("bed.gz", "59c7b1385c96350fd2a7217f597ad495", "replicated peaks"),      # H3K36me3 K562
    "ENCFF029GQD": ("bed.gz", "e8cce2156da0a8c7a8d59bfc7d628ba3", "replicated peaks"),      # H3K36me3 HCT116
    "ENCFF266IRD": ("bigWig", "02d5d1553d06f0fcf45a78501ca64e0c", "plus strand signal of unique reads"),
    "ENCFF209TJK": ("bigWig", "d34834133c0528abe64286e3a9d75182", "plus strand signal of unique reads"),
    "ENCFF116OZK": ("bigWig", "d7d15531778d1702310276065cd285e7", "minus strand signal of unique reads"),
    "ENCFF828GNJ": ("bigWig", "592a3c0a8f7e3ddd6092c582cd7ac2e0", "minus strand signal of unique reads"),
    "ENCFF440UVY": ("bigWig", "9d22de120814a72ce719d8a3274f3694", "plus strand signal of unique reads"),
    "ENCFF078CJH": ("bigWig", "7aa4cfb8fcaefb8b55f7e860156c83f3", "plus strand signal of unique reads"),
    "ENCFF103FBH": ("bigWig", "847a7cb26342b4f68bc7d92a20ae22c2", "minus strand signal of unique reads"),
    "ENCFF907PJM": ("bigWig", "b10d937242984712bbb32ca0ba59540f", "minus strand signal of unique reads"),
}
SETS = {"266IRD": {"plus": ["ENCFF266IRD", "ENCFF209TJK"], "minus": ["ENCFF116OZK", "ENCFF828GNJ"]},
        "440UVY": {"plus": ["ENCFF440UVY", "ENCFF078CJH"], "minus": ["ENCFF103FBH", "ENCFF907PJM"]}}
FALLBACK = "266IRD"


WAITS = (30, 60, 120, 240, 480)


def with_retries(what, fn):
    """Call fn(); on a server error (5xx), a timeout or a failed connection wait and try again (6 attempts)."""
    for k in range(len(WAITS) + 1):
        try:
            return fn()
        except urllib.error.HTTPError as ex:
            if ex.code < 500 or k == len(WAITS):
                raise
            err = f"HTTP {ex.code}"
        except (urllib.error.URLError, TimeoutError, ConnectionError) as ex:
            if k == len(WAITS):
                raise
            err = str(ex)
        print(f"  {what}: {err}; new attempt in {WAITS[k]} s", flush=True)
        time.sleep(WAITS[k])


def get(path):
    def once():
        req = urllib.request.Request(BASE + path, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=180) as r:
            return json.load(r)
    return with_retries(path, once)


def fetch(url, part):
    def once():
        with urllib.request.urlopen(url, timeout=600) as src, open(part, "wb") as out:
            while True:
                b = src.read(1 << 22)
                if not b:
                    break
                out.write(b)
    with_retries(url, once)


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    data, res = sys.argv[1], sys.argv[2]
    os.makedirs(data, exist_ok=True)
    os.makedirs(res, exist_ok=True)
    print("== metadata (ENCODE REST API) ==")
    for acc, (ext, md5, ot) in FILES.items():
        F = get(f"/files/{acc}/?frame=object")
        got = (F.get("md5sum"), F.get("assembly"), F.get("status"), F.get("output_type"))
        print(f"  {acc}: md5 {got[0]}, {got[1]}, {got[2]}, {got[3]}")
        if got != (md5, "GRCh38", "released", ot):
            sys.exit(f"STOP: {acc}: metadata differ from those recorded on 10 Nov 2026")

    E = get("/experiments/ENCSR165QKR/?frame=object")
    da = E.get("default_analysis")
    in_default = set()
    if da:
        in_default = set(get(da + "?frame=object").get("files", []))
    hit = [k for k, v in SETS.items() if all(f"/files/{x}/" in in_default for x in v["plus"] + v["minus"])]
    if len(hit) == 1:
        prim, rule = hit[0], f"ENCODE default analysis ({da})"
    else:
        prim, rule = FALLBACK, (f"default analysis {da or 'not given'} contains {len(hit)} of the two sets; pre-declared "
                                "fallback: the set processed like the K562 data")
    alt = [k for k in SETS if k != prim][0]
    print(f"\n== HCT116 bigWig sets ==\n  primary: set {prim} ({rule})\n  sensitivity: set {alt}")
    with open(os.path.join(res, "28_hct116_sets.tsv"), "w") as fh:
        fh.write("role\tset\tplus\tminus\trule\n")
        for role, k in (("primary", prim), ("sensitivity", alt)):
            fh.write(f"{role}\t{k}\t{','.join(SETS[k]['plus'])}\t{','.join(SETS[k]['minus'])}\t{rule if role == 'primary' else ''}\n")

    print("\n== downloads ==")
    rows = []
    for acc, (ext, md5, _) in FILES.items():
        dst = os.path.join(data, f"{acc}.{ext}")
        if os.path.isfile(dst) and md5_of(dst) == md5:
            print(f"  {acc}: present, md5 OK")
        else:
            t0 = time.time()
            part = dst + ".part"
            fetch(f"{BASE}/files/{acc}/@@download/{acc}.{ext}", part)
            if md5_of(part) != md5:
                os.remove(part)
                sys.exit(f"STOP: {acc}: md5 of the download differs from ENCODE's")
            os.replace(part, dst)
            print(f"  {acc}: downloaded, md5 OK  [{time.time() - t0:.0f} s]")
        rows.append(f"{acc}\t{ext}\t{md5}\t{os.path.getsize(dst)}\t{dst}\t{time.strftime('%Y-%m-%d')}")
    with open(os.path.join(res, "28_downloads.tsv"), "w") as fh:
        fh.write("accession\tformat\tmd5\tbytes\tpath\tchecked_on\n" + "\n".join(rows) + "\n")
    print("written: 28_downloads.tsv, 28_hct116_sets.tsv")


if __name__ == "__main__":
    main()
