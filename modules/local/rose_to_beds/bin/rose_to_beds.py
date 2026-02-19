#!/usr/bin/env python3
from __future__ import annotations
import argparse, os, sys


def eprint(*a, **k):
    print(*a, file=sys.stderr, **k)


def detect_header_and_rows(path: str):
    """
    Reads a ROSE stitched table that may contain leading '#' comment lines.
    Assumes tab-delimited. First non-# line is header.
    """
    header = None
    rows = []
    with open(path, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            fields = [c.strip() for c in line.split("\t")]
            if header is None:
                header = fields
                continue
            if len(fields) < len(header):
                fields += [""] * (len(header) - len(fields))
            rows.append(fields[: len(header)])
    if header is None:
        raise ValueError(f"Could not find header in: {path}")
    return header, rows


def parse_rose_stitched_table(path: str):
    """
    Returns (SEs_bed6, TEs_bed6) as list-of-lists with BED6 fields:
    chrom start end name score strand
    """
    header, rows = detect_header_and_rows(path)
    idx = {c: i for i, c in enumerate(header)}

    required = ("REGION_ID", "CHROM", "START", "STOP", "isSuper")
    for col in required:
        if col not in idx:
            raise ValueError(
                f"Missing required column '{col}' in {path}. Found columns: {header}"
            )

    rank_col = "stitchedPeakRank" if "stitchedPeakRank" in idx else None

    def to_int(x: str) -> int:
        x = (x or "").strip()
        return int(float(x))

    def to_bool01(x: str) -> int:
        x = (x or "").strip()
        if x in {"1", "True", "true", "YES", "yes"}:
            return 1
        if x in {"0", "False", "false", "NO", "no"}:
            return 0
        return 1 if int(float(x)) != 0 else 0

    SEs, TEs = [], []

    for r in rows:
        try:
            chrom = r[idx["CHROM"]]
            start = to_int(r[idx["START"]])
            end = to_int(r[idx["STOP"]])
            is_super = to_bool01(r[idx["isSuper"]])
            region_id = r[idx["REGION_ID"]] or "."
            score = "0"
            if rank_col is not None:
                score = r[idx[rank_col]] or "0"
        except Exception as ex:
            eprint(f"[WARN] Skipping bad row: {r} ({ex})")
            continue

        if not chrom or end < start:
            continue

        bed6 = [chrom, str(start), str(end), region_id, str(score), "."]
        if is_super == 1:
            SEs.append(bed6)
        else:
            TEs.append(bed6)

    return SEs, TEs


def read_peaks_bed(path: str):
    """
    Reads BED3+ peaks; preserves extra columns.
    Returns list of string lists (entire row).
    """
    peaks = []
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if (
                not line
                or line.startswith("#")
                or line.startswith("track")
                or line.startswith("browser")
            ):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                int(parts[1]); int(parts[2])
            except Exception:
                eprint(f"[WARN] Skipping bad peak row: {line}")
                continue
            peaks.append(parts)
    return peaks


def overlaps_bed3(a_chr, a_s, a_e, b_chr, b_s, b_e) -> bool:
    return a_chr == b_chr and a_s < b_e and b_s < a_e


def build_binned_index(intervals_bed6, bin_size: int):
    """
    intervals_bed6: [chrom, start, end, ...]
    Returns dict[(chrom, bin)] -> list(interval_rows)
    """
    idx = {}
    for iv in intervals_bed6:
        chrom = iv[0]
        s = int(iv[1])
        e = int(iv[2])
        b0 = s // bin_size
        b1 = e // bin_size
        for b in range(b0, b1 + 1):
            idx.setdefault((chrom, b), []).append(iv)
    return idx


def constituent_peaks(peaks_rows, stitched_bed6, bin_size: int):
    idx = build_binned_index(stitched_bed6, bin_size)
    out = []
    for p in peaks_rows:
        chrom = p[0]
        s = int(p[1])
        e = int(p[2])
        b0 = s // bin_size
        b1 = e // bin_size

        hit = False
        for b in range(b0, b1 + 1):
            for iv in idx.get((chrom, b), []):
                if overlaps_bed3(chrom, s, e, iv[0], int(iv[1]), int(iv[2])):
                    hit = True
                    break
            if hit:
                break

        if hit:
            out.append(p)
    return out


def write_rows(path: str, rows):
    with open(path, "w") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rose", required=True, help="ROSE stitched region map table (txt)")
    ap.add_argument("--peaks", required=True, help="Peaks BED (BED3+)")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--prefix", default="", help="Output filename prefix (e.g. sample1)")
    ap.add_argument("--bin-size", type=int, default=100000, help="Index bin size")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    SEs, TEs = parse_rose_stitched_table(args.rose)
    peaks = read_peaks_bed(args.peaks)

    se_const = constituent_peaks(peaks, SEs, args.bin_size)
    te_const = constituent_peaks(peaks, TEs, args.bin_size)

    pfx = f"{args.prefix}_" if args.prefix else ""
    write_rows(os.path.join(args.outdir, f"{pfx}SEs.bed"), SEs)
    write_rows(os.path.join(args.outdir, f"{pfx}TEs.bed"), TEs)
    write_rows(os.path.join(args.outdir, f"{pfx}SEs_constituents.bed"), se_const)
    write_rows(os.path.join(args.outdir, f"{pfx}TEs_constituents.bed"), te_const)


if __name__ == "__main__":
    main()
