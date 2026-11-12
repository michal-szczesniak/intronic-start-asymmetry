#!/usr/bin/env python3
"""Step 16a: fetch the ProCapNet K562 models from ENCODE and check that they load.

ProCapNet (Cochran et al. 2024) was trained with 7-fold cross-validation split by
chromosome (test chromosomes of folds 1..7: chr1,4 | chr2,13,16 | chr5,6,20,21 |
chr7,8,9 | chr10,11,12 | chr3,14,15,17 | chr18,19,22,X,Y). Scoring every iSTART pair
with the model whose test fold contains the pair's chromosome gives sequence
predictions that never saw the measured signal at that locus, for all pairs instead of
the ~150-200 pairs on chr8/9 that Puffin-D leaves us. This step only gets the models:

  1. lists the files of the ENCODE annotation (REST API) and prints their metadata;
  2. downloads the released files of the selected output types (--types; default the
     models and the training/test regions), resumable, md5 checked against ENCODE, and
     unpacks archives;
  3. tries to load every model-like file with torch on the CPU, with the ProCapNet code
     on the import path (the models were saved as whole pickled objects), runs one
     forward pass on a random 2114-nt sequence and prints what comes out.
The printout tells how the files map to folds and what the outputs look like, which the
scoring step (16b) needs.

Note: torch.load(weights_only=False) runs the pickle in the file; use it only on files
from a source you trust (here: the ENCODE portal, Kundaje lab).

Usage:
  python3 16a_get_procapnet.py --accession ENCSR740IPL --outdir $SCRATCH/models/procapnet_K562 \
      --code $SCRATCH/tools/ProCapNet [--base-url https://www.encodeproject.org] [--no-load]
"""
import argparse
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request

MODEL_SUFFIXES = (".model", ".torch", ".pt", ".pth", ".pkl")


def get_json(url):
    if url.startswith("file://"):
        with open(url[len("file://"):].split("?")[0]) as fh:
            return json.load(fh)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode())


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def file_records(ann, base):
    """File objects of the annotation; fetch any that are given only as @id paths."""
    out, seen = [], set()
    for key in ("files", "original_files"):
        for f in ann.get(key, []) or []:
            if isinstance(f, str):
                if f in seen:
                    continue
                try:
                    f = get_json(base + f + ("" if base.startswith("file://") else "?format=json"))
                except Exception as e:                                   # noqa: BLE001
                    print(f"  could not fetch {f}: {e}")
                    continue
            fid = f.get("@id", f.get("accession"))
            if fid in seen:
                continue
            seen.add(fid)
            out.append(f)
    return out


def safe_extract(tar_path, dest):
    with tarfile.open(tar_path) as t:
        for m in t.getmembers():
            p = os.path.realpath(os.path.join(dest, m.name))
            root = os.path.realpath(dest)
            if p != root and not p.startswith(root + os.sep):
                sys.exit(f"STOP: unsafe path in {tar_path}: {m.name}")
            if m.issym() or m.islnk():
                sys.exit(f"STOP: link in {tar_path}: {m.name}")
        try:
            t.extractall(dest, filter="data")
        except TypeError:                            # Python without extraction filters
            t.extractall(dest)


def describe_output(x, depth=0):
    import torch
    pad = "    " * depth
    if isinstance(x, torch.Tensor):
        v = x.detach().float()
        return f"{pad}tensor {tuple(x.shape)} min {v.min():.3g} max {v.max():.3g} mean {v.mean():.3g}"
    if isinstance(x, (list, tuple)):
        return f"{pad}{type(x).__name__} of {len(x)}:\n" + "\n".join(describe_output(y, depth + 1) for y in x)
    return f"{pad}{type(x).__name__}"


def try_load(path, code):
    import numpy as np
    import torch
    for p in (os.path.join(code, "src", "2_train_models"), os.path.join(code, "src", "utils")):
        if os.path.isdir(p) and p not in sys.path:
            sys.path.insert(0, p)
    try:
        obj = torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        obj = torch.load(path, map_location="cpu")
    print(f"    loaded: {type(obj).__module__}.{type(obj).__name__}")
    if isinstance(obj, dict):
        ks = list(obj.keys())
        print(f"    dict with {len(ks)} keys; first: " + ", ".join(f"{k} {tuple(getattr(obj[k], 'shape', ()))}" for k in ks[:8]))
        return
    if isinstance(obj, torch.nn.Module):
        obj.eval()
        n = sum(p.numel() for p in obj.parameters())
        print(f"    parameters: {n:,}; attributes: " + ", ".join(
            f"{a}={getattr(obj, a)}" for a in ("n_filters", "n_layers", "n_outputs", "trimming") if hasattr(obj, a)))
        rng = np.random.default_rng(0)
        x = np.zeros((1, 4, 2114), dtype=np.float32)
        x[0, rng.integers(0, 4, 2114), np.arange(2114)] = 1
        with torch.no_grad():
            y = obj(torch.tensor(x))
        print("    forward on one random 2114-nt sequence:\n" + describe_output(y, 2))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", default="ENCSR740IPL")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--code", required=True, help="clone of github.com/kundajelab/ProCapNet")
    ap.add_argument("--base-url", default="https://www.encodeproject.org")
    ap.add_argument("--types", default="models,training and test regions",
                    help="ENCODE output_type values to download (comma-separated)")
    ap.add_argument("--no-load", action="store_true")
    a = ap.parse_args()
    types = {t.strip() for t in a.types.split(",") if t.strip()}
    base = a.base_url.rstrip("/")
    os.makedirs(a.outdir, exist_ok=True)

    ann_url = f"{base}/annotations/{a.accession}/" + ("index.json" if base.startswith("file://") else "?format=json")
    ann = get_json(ann_url)
    print(f"== {a.accession}: {ann.get('description', '')} | status {ann.get('status')} | "
          f"type {ann.get('annotation_type')} | biosample {ann.get('biosample_ontology', {}).get('term_name', '') if isinstance(ann.get('biosample_ontology'), dict) else ann.get('biosample_ontology')}")
    recs = file_records(ann, base)
    if not recs:
        sys.exit("STOP: no files listed for this annotation (see the JSON above or the portal page)")
    print(f"\n== {len(recs)} files ==")
    for f in recs:
        print(f"  {f.get('accession')}  {f.get('status')}  {f.get('file_format')}  {f.get('output_type')}  "
              f"{(f.get('file_size') or 0) / 1e6:.1f} MB  submitted: {f.get('submitted_file_name', '')}  href: {f.get('href', '')}")

    got = []
    for f in recs:
        if f.get("status") not in ("released", "archived") or not f.get("href") or f.get("output_type") not in types:
            continue
        dest = os.path.join(a.outdir, os.path.basename(f["href"]))
        want = f.get("md5sum")
        if os.path.isfile(dest) and want and md5(dest) == want:
            print(f"  {os.path.basename(dest)} already here, md5 ok")
        else:
            if os.path.isfile(dest):
                os.remove(dest)                      # a finished file that fails the check: fetch again
            part = dest + ".part"
            print(f"  downloading {os.path.basename(dest)} ...", flush=True)
            for attempt in (1, 2):
                r = subprocess.run(["curl", "-sSL", "--retry", "5", "-C", "-", "-o", part, base + f["href"]])
                if r.returncode == 0 and (not want or md5(part) == want):
                    break
                print(f"    attempt {attempt} failed (curl exit {r.returncode}); starting from scratch")
                if os.path.isfile(part):
                    os.remove(part)
            else:
                sys.exit(f"STOP: could not download {dest} with the right md5")
            os.replace(part, dest)
            print(f"    md5 {'ok' if want else 'not given by ENCODE'}")
        got.append(dest)

    for dest in got:
        if tarfile.is_tarfile(dest):
            ex = os.path.join(a.outdir, os.path.basename(dest).split(".")[0])
            done = os.path.join(ex, ".unpacked_ok")
            if not os.path.isfile(done):             # unpack afresh (also after an interrupted attempt)
                if os.path.isdir(ex):
                    shutil.rmtree(ex)
                os.makedirs(ex)
                safe_extract(dest, ex)
                open(done, "w").close()
            print(f"\n== contents of {os.path.basename(dest)} (unpacked in {ex}) ==")
            for root, _, files in sorted(os.walk(ex)):
                for fn in sorted(files):
                    if fn == ".unpacked_ok":
                        continue
                    p = os.path.join(root, fn)
                    print(f"  {os.path.relpath(p, a.outdir)}  {os.path.getsize(p) / 1e6:.2f} MB")

    if a.no_load:
        return
    cands = sorted(p for p in glob.glob(os.path.join(a.outdir, "**", "*"), recursive=True)
                   if os.path.isfile(p) and p.endswith(MODEL_SUFFIXES))
    print(f"\n== loading {len(cands)} model-like files on the CPU ==")
    for p in cands:
        print(f"  {os.path.relpath(p, a.outdir)}")
        try:
            try_load(p, a.code)
        except Exception as e:                                           # noqa: BLE001
            print(f"    FAILED: {type(e).__name__}: {e}")
    for p in sorted(glob.glob(os.path.join(a.outdir, "**", "*.json"), recursive=True))[:7]:
        try:
            with open(p) as fh:
                print(f"\n  {os.path.relpath(p, a.outdir)}: {fh.read(400)}")
        except OSError:
            pass


if __name__ == "__main__":
    main()
