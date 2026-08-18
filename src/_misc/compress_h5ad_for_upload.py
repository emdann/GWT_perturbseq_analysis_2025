#!/usr/bin/env python3
"""
Recompress an uncompressed .assigned_guide.h5ad in place-safe fashion (new file).

The files written by merge_sample_for_upload.py use anndata's default
write_h5ad(), which applies no compression, so X/data and X/indices are stored
raw (~110 GiB/file). Applying gzip + byte-shuffle shrinks them ~6x with no
change to dtypes or structure, and keeps the result readable by any
h5py/anndata without filter plugins.

The matrices are far too large to hold in RAM, so X is streamed slab by slab
with h5py; everything else (obs, var, uns, ...) is copied verbatim.

Usage:
    python compress_h5ad_for_upload.py --sample D4_Rest
    python compress_h5ad_for_upload.py --sample D4_Rest --verify-only
"""

import argparse
import os
import sys
import time
import zlib

import h5py
import numpy as np

# gzip level 4 is roughly throughput-matched to oak's ~75 MiB/s read rate;
# higher levels become CPU-bound without buying much on this data.
GZIP_LEVEL = 4
TARGET_CHUNK_BYTES = 1 << 20   # 1 MiB chunks (source used ~53 KiB -> 1.5M chunks)
SLAB_ELEMENTS = 32_000_000     # ~256 MiB per read for int64


def fmt_size(nbytes):
    return f"{nbytes / 2**30:.1f} GiB"


def fmt_dur(seconds):
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h:d}h{m:02d}m{s:02d}s" if h else f"{m:d}m{s:02d}s"


def stream_dataset(src, dst_group, name, crc_out):
    """Copy a 1-D dataset with compression, accumulating a CRC32 of the source."""
    n = src.shape[0]
    chunk = min(n, max(1, TARGET_CHUNK_BYTES // src.dtype.itemsize))
    dst = dst_group.create_dataset(
        name,
        shape=src.shape,
        dtype=src.dtype,
        chunks=(chunk,),
        compression="gzip",
        compression_opts=GZIP_LEVEL,
        shuffle=True,
    )
    for k, v in src.attrs.items():
        dst.attrs[k] = v

    crc = 0
    t0 = time.time()
    n_slabs = (n + SLAB_ELEMENTS - 1) // SLAB_ELEMENTS
    for i, start in enumerate(range(0, n, SLAB_ELEMENTS)):
        stop = min(start + SLAB_ELEMENTS, n)
        buf = src[start:stop]
        dst[start:stop] = buf
        crc = zlib.crc32(buf.tobytes(), crc)
        if i % 10 == 0 or stop == n:
            elapsed = time.time() - t0
            done = stop * src.dtype.itemsize
            rate = done / 2**20 / elapsed if elapsed else 0
            eta = (n - stop) * src.dtype.itemsize / 2**20 / rate if rate else 0
            print(
                f"    {name}: {100 * stop / n:5.1f}%  {fmt_size(done)}  "
                f"{rate:.0f} MiB/s  ETA {fmt_dur(eta)}",
                flush=True,
            )
    crc_out[name] = crc
    return dst


def attrs_equal(a, b):
    """Compare HDF5 attrs dicts, tolerating array-valued entries (e.g. X 'shape')."""
    if sorted(a.keys()) != sorted(b.keys()):
        return False
    return all(np.array_equal(np.asarray(a[k]), np.asarray(b[k])) for k in a)


def checksum_dataset(dset):
    """CRC32 over a dataset, read slab-wise."""
    crc = 0
    n = dset.shape[0]
    for start in range(0, n, SLAB_ELEMENTS):
        stop = min(start + SLAB_ELEMENTS, n)
        crc = zlib.crc32(dset[start:stop].tobytes(), crc)
    return crc


def compress(src_path, dst_path):
    print(f"source: {src_path}")
    print(f"target: {dst_path}")
    src_bytes = os.path.getsize(src_path)
    print(f"source size: {fmt_size(src_bytes)}")

    crcs = {}
    t0 = time.time()
    with h5py.File(src_path, "r") as fin, h5py.File(dst_path, "w") as fout:
        for k, v in fin.attrs.items():
            fout.attrs[k] = v

        # Everything except X is small -- copy verbatim (preserves encoding attrs).
        for key in fin.keys():
            if key == "X":
                continue
            print(f"  copying /{key} ...", flush=True)
            fin.copy(fin[key], fout, name=key)

        gin = fin["X"]
        if isinstance(gin, h5py.Dataset):
            raise SystemExit("X is dense; this script only handles CSR/CSC X")
        gout = fout.create_group("X")
        for k, v in gin.attrs.items():
            gout.attrs[k] = v

        print(f"  X: {dict(gin.attrs).get('encoding-type')} "
              f"shape={list(gin.attrs.get('shape', []))} nnz={gin['data'].shape[0]:,}")
        for name in ("data", "indices", "indptr"):
            stream_dataset(gin[name], gout, name, crcs)

    dst_bytes = os.path.getsize(dst_path)
    print(f"\nwrote {fmt_size(dst_bytes)} in {fmt_dur(time.time() - t0)} "
          f"({src_bytes / dst_bytes:.2f}x smaller)")
    return crcs


def verify(src_path, dst_path, crcs=None):
    """Full integrity check of the compressed copy against the source."""
    print("\nverifying ...")
    ok = True
    with h5py.File(src_path, "r") as fin, h5py.File(dst_path, "r") as fout:
        # Structure / attrs
        if sorted(fin.keys()) != sorted(fout.keys()):
            print(f"  FAIL root keys: {sorted(fin.keys())} != {sorted(fout.keys())}")
            ok = False
        if not attrs_equal(fin["X"].attrs, fout["X"].attrs):
            print("  FAIL X attrs differ")
            ok = False
        if not attrs_equal(fin.attrs, fout.attrs):
            print("  FAIL root attrs differ")
            ok = False

        # obs/var frames, read via anndata's own reader for a semantic check
        import anndata.io
        for key in ("obs", "var"):
            a = anndata.io.read_elem(fin[key])
            b = anndata.io.read_elem(fout[key])
            same = a.equals(b)
            print(f"  {key}: {'OK' if same else 'FAIL'} ({a.shape[0]:,} rows, "
                  f"{a.shape[1]} cols)")
            ok &= same

        # X payload: CRC32 of the compressed copy vs. the source
        for name in ("data", "indices", "indptr"):
            din, dout = fin["X"][name], fout["X"][name]
            if din.shape != dout.shape or din.dtype != dout.dtype:
                print(f"  FAIL X/{name}: {din.shape}/{din.dtype} vs "
                      f"{dout.shape}/{dout.dtype}")
                ok = False
                continue
            expected = crcs.get(name) if crcs else checksum_dataset(din)
            actual = checksum_dataset(dout)
            same = expected == actual
            print(f"  X/{name}: {'OK' if same else 'FAIL'} crc={actual:08x} "
                  f"({din.dtype}, {din.shape[0]:,} elements)")
            ok &= same

    print("VERIFY PASSED" if ok else "VERIFY FAILED")
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True, help="e.g. D4_Rest")
    p.add_argument("--indir", default="/mnt/oak/users/emma/data/GWT/to_share/")
    p.add_argument("--outdir", default="/mnt/oak/users/emma/data/GWT/to_share_compressed/")
    p.add_argument("--verify-only", action="store_true",
                   help="Only verify an existing compressed file")
    args = p.parse_args()

    src = os.path.join(args.indir, f"{args.sample}.assigned_guide.h5ad")
    dst = os.path.join(args.outdir, f"{args.sample}.assigned_guide.h5ad")
    if not os.path.exists(src):
        sys.exit(f"ERROR: no such file {src}")

    if args.verify_only:
        sys.exit(0 if verify(src, dst) else 1)

    os.makedirs(args.outdir, exist_ok=True)
    if os.path.exists(dst):
        sys.exit(f"ERROR: {dst} exists; remove it or use --verify-only")

    crcs = compress(src, dst)
    if not verify(src, dst, crcs):
        sys.exit(1)
    print(f"\nSUCCESS: {args.sample}")


if __name__ == "__main__":
    main()
