#!/usr/bin/env python3
"""
ROSE stitched region map -> SE/TE BEDs + constituent peaks by overlap.

Example:
  python rose_to_beds.py \
    --rose ENCFF493UFO_12KB_STITCHED_REGION_MAP.txt \
    --peaks peaks.bed \
    --outdir .

Outputs:
  SEs.bed
  TEs.bed
  SEs_constituents.bed
  TEs_constituents.bed
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    name: str = "."
    extra: Tuple[str, ...] = ()

    def as_bed6(self) -> List[str]:
        # BED: chrom, start, end, name, score, strand (we use '.' placeholders)
        score = self.extra[0] if len(self.extra) >= 1 else "0"
        strand = self.extra[1] if len(self.extra) >= 2 else "."
        return [self.chrom, str(self.start), str(self.end), self.name, str(score), strand]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def detect_header_and_reader(path: str) -> Tuple[List[str], List[List[str]]]:
    """
    Reads a tab-delimited ROSE stitched file that may contain '#' comment lines.
    Returns (header, rows) where rows are lists of strings aligned to header.
    """
    header: Optional[List[str]] = None
    rows: List[List[str]] = []

    with open(path, "r", newline="") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            # Identify header row by presence of known columns
            if header is None:
                # allow header even if spaces; normalize
                norm = [c.strip() for c in fields]
                must = {"CHROM", "START", "STOP", "isSuper"}
                if must.issubset(set(norm)):
                    header = norm
                    continue
                else:
                    # If file has no explicit header (rare), assume ROSE standard
                    header = norm  # treat first non-# line as header anyway
                    continue

            rows.append([c.strip() for c in fields])

    if header is None:
        raise ValueError(f"Could not find a header in: {path}")

    # Ensure all rows have at least header length (pad if needed)
    fixed_rows = []
    for r in rows:
        if len(r) < len(header):
            r = r + [""] * (len(header) - len(r))
        fixed_rows.append(r[: len(header)])
    return header, fixed_rows


def parse_rose_stitched(path: str) -> Tuple[List[Interval], List[Interval]]:
    header, rows = detect_header_and_reader(path)
    idx: Dict[str, int] = {col: i for i, col in enumerate(header)}

    # Required columns
    for col in ("REGION_ID", "CHROM", "START", "STOP", "isSuper"):
        if col not in idx:
            raise ValueError(
                f"Missing required column '{col}' in {path}. Found columns: {header}"
            )

    # Optional columns to carry into name/score
    # We'll store rank in score if available.
    rank_col = "stitchedPeakRank" if "stitchedPeakRank" in idx else None

    SEs: List[Interval] = []
    TEs: List[Interval] = []

    for r in rows:
        chrom = r[idx["CHROM"]]
        if not chrom:
            continue

        # Coerce ints robustly (strip, allow floats)
        def to_int(x: str) -> int:
            x = (x or "").strip()
            if x == "":
                raise ValueError("empty int")
            # some files may contain floats; cast via float first
            return int(float(x))

        def to_bool01(x: str) -> int:
            x = (x or "").strip()
            if x in {"1", "True", "true", "YES", "yes"}:
                return 1
            if x in {"0", "False", "false", "NO", "no"}:
                return 0
            # fallback: numeric cast
            return 1 if int(float(x)) != 0 else 0

        try:
            start = to_int(r[idx["START"]])
            end = to_int(r[idx["STOP"]])
            is_super = to_bool01(r[idx["isSuper"]])
        except Exception as ex:
            eprint(f"[WARN] Skipping bad row (parse error): {r} ({ex})")
            continue

        if end < start:
            eprint(f"[WARN] Skipping bad interval (end < start): {r}")
            continue

        region_id = r[idx["REGION_ID"]] or "."
        score = "0"
        if rank_col is not None:
            val = r[idx[rank_col]]
            score = val if val != "" else "0"

        iv = Interval(chrom=chrom, start=start, end=end, name=region_id, extra=(score, "."))

        if is_super == 1:
            SEs.append(iv)
        else:
            TEs.append(iv)

    return SEs, TEs


def read_peaks_bed(path: str) -> List[Interval]:
    peaks: List[Interval] = []
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except Exception:
                eprint(f"[WARN] Skipping bad BED row: {line}")
                continue
            name = parts[3] if len(parts) >= 4 and parts[3] else "."
            # Keep all remaining columns to re-emit faithfully
            extra = tuple(parts[4:]) if len(parts) > 4 else ()
            peaks.append(Interval(chrom, start, end, name=name, extra=extra))
    return peaks


def build_interval_index(intervals: List[Interval], bin_size: int = 100_000) -> Dict[Tuple[str, int], List[Interval]]:
    """
    Simple binned index: (chrom, bin_id) -> intervals.
    Good enough for peak x stitched overlaps.
    """
    index: Dict[Tuple[str, int], List[Interval]] = {}
    for iv in intervals:
        b0 = iv.start // bin_size
        b1 = iv.end // bin_size
        for b in range(b0, b1 + 1):
            index.setdefault((iv.chrom, b), []).append(iv)
    return index


def overlaps(a: Interval, b: Interval) -> bool:
    # half-open BED semantics [start, end)
    return a.chrom == b.chrom and a.start < b.end and b.start < a.end


def constituent_peaks(peaks: List[Interval], stitched: List[Interval], bin_size: int = 100_000) -> List[Interval]:
    idx = build_interval_index(stitched, bin_size=bin_size)
    out: List[Interval] = []
    for p in peaks:
        b0 = p.start // bin_size
        b1 = p.end // bin_size
        hit = False
        for b in range(b0, b1 + 1):
            cand = idx.get((p.chrom, b), [])
            for s in cand:
                if overlaps(p, s):
                    hit = True
                    break
            if hit:
                break
        if hit:
            out.append(p)
    return out


def write_bed(path: str, intervals: List[Interval], bed6: bool = False):
    with open(path, "w") as f:
        for iv in intervals:
            if bed6:
                f.write("\t".join(iv.as_bed6()) + "\n")
            else:
                # emit as at least BED3; if peak had extra cols, keep them
                row = [iv.chrom, str(iv.start), str(iv.end)]
                if iv.name != "." or iv.extra:
                    row.append(iv.name)
                if iv.extra:
                    row.extend(list(iv.extra))
                f.write("\t".join(row) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rose", required=True, help="ROSE stitched region map txt (e.g., *_STITCHED_REGION_MAP.txt)")
    ap.add_argument("--peaks", required=True, help="Input peaks BED file (BED3+)")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--bin-size", type=int, default=100_000, help="Bin size for overlap indexing (default: 100000)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    SEs, TEs = parse_rose_stitched(args.rose)
    peaks = read_peaks_bed(args.peaks)

    se_const = constituent_peaks(peaks, SEs, bin_size=args.bin_size)
    te_const = constituent_peaks(peaks, TEs, bin_size=args.bin_size)

    se_bed = os.path.join(args.outdir, "SEs.bed")
    te_bed = os.path.join(args.outdir, "TEs.bed")
    se_const_bed = os.path.join(args.outdir, "SEs_constituents.bed")
    te_const_bed = os.path.join(args.outdir, "TEs_constituents.bed")

    # Stitched outputs as BED6 so you keep a useful name + rank as score
    write_bed(se_bed, SEs, bed6=True)
    write_bed(te_bed, TEs, bed6=True)

    # Constituents: preserve original peak columns as much as possible
    write_bed(se_const_bed, se_const, bed6=False)
    write_bed(te_const_bed, te_const, bed6=False)

    eprint(f"[OK] Wrote {len(SEs)} SE stitched regions -> {se_bed}")
    eprint(f"[OK] Wrote {len(TEs)} TE stitched regions -> {te_bed}")
    eprint(f"[OK] Wrote {len(se_const)} SE constituent peaks -> {se_const_bed}")
    eprint(f"[OK] Wrote {len(te_const)} TE constituent peaks -> {te_const_bed}")


if __name__ == "__main__":
    main()
