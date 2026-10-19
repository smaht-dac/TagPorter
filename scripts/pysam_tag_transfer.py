#!/usr/bin/env python3
"""
In-memory pysam-based tag transfer script.

Used as a comparison baseline for TagPorter benchmarks (B3).
All requested tags from source BAM(s) are loaded into a Python dict keyed by
read name, then transferred to the aligned BAM.

Timing and stats are printed to stderr in the same format as TagPorter
so the same parsing scripts work on both tools.

Usage:
    python pysam_tag_transfer.py \
        --aln aligned.bam \
        --source-bams unaligned.bam \
        --output output.bam \
        --tags MM,ML
"""

import argparse
import array as array_module
import resource
import sys
import time

import pysam


def parse_args():
    p = argparse.ArgumentParser(
        description="In-memory BAM tag transfer (pysam baseline for benchmarking)"
    )
    p.add_argument("--aln", required=True, help="Aligned BAM to add tags to")
    p.add_argument(
        "--source-bams",
        nargs="+",
        required=True,
        metavar="BAM",
        help="Source BAM file(s) containing the tags to transfer",
    )
    p.add_argument("--output", required=True, help="Output BAM path")
    p.add_argument(
        "--tags",
        default="MM,ML",
        help="Comma-separated 2-character tag names (default: MM,ML)",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for pysam BAM I/O (default: 1)",
    )
    return p.parse_args()


def peak_ram_mb():
    """Peak RSS in MB. On Linux ru_maxrss is in KB; on macOS in bytes."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / 1024 / 1024
    return rss / 1024  # Linux: KB -> MB


def main():
    args = parse_args()
    requested = set(args.tags.split(","))

    t0 = time.monotonic()

    # ── Phase 1: collect ──────────────────────────────────────────────────────
    # tag_store: read_name -> list of (tag, value, type_code)
    tag_store: dict[str, list[tuple[str, object, str]]] = {}

    for source_path in args.source_bams:
        n = 0
        with pysam.AlignmentFile(
            source_path, "rb", check_sq=False, threads=args.threads
        ) as src:
            for read in src:
                collected = [
                    (tag, val, vtype)
                    for tag, val, vtype in read.get_tags(with_value_type=True)
                    if tag in requested
                ]
                if collected:
                    tag_store[read.query_name] = collected
                n += 1
                if n % 100_000 == 0:
                    print(
                        f"Collected {n} reads from {source_path}",
                        file=sys.stderr,
                    )
        print(
            f"Finished collecting from {source_path}: {n} reads processed",
            file=sys.stderr,
        )

    t_collect = time.monotonic() - t0

    # ── Phase 2: transfer ─────────────────────────────────────────────────────
    t1 = time.monotonic()
    total = 0
    reads_tagged = 0
    reads_unmatched = 0

    with pysam.AlignmentFile(
        args.aln, "rb", threads=args.threads
    ) as aln, pysam.AlignmentFile(
        args.output, "wb", header=aln.header, threads=args.threads
    ) as out:
        for read in aln:
            total += 1
            stored = tag_store.get(read.query_name)
            if stored is not None:
                for tag, val, vtype in stored:
                    if isinstance(val, array_module.array):
                        read.set_tag(tag, val)
                    else:
                        read.set_tag(tag, val, value_type=vtype)
                reads_tagged += 1
            else:
                reads_unmatched += 1
            out.write(read)

            if total % 100_000 == 0:
                print(f"Processed {total} aligned records", file=sys.stderr)

    pysam.index(args.output, threads=args.threads)

    t_write = time.monotonic() - t1
    t_total = time.monotonic() - t0

    # ── Stats ─────────────────────────────────────────────────────────────────
    ram = peak_ram_mb()
    pct = 100 * reads_tagged / total if total else 0.0
    pct_u = 100 * reads_unmatched / total if total else 0.0

    print(
        f"Phase timing: collect={t_collect:.1f}s  write={t_write:.1f}s  total={t_total:.1f}s",
        file=sys.stderr,
    )
    print(f"Transfer complete: {total} aligned reads processed", file=sys.stderr)
    print(f"  {reads_tagged} reads received tags ({pct:.1f}%)", file=sys.stderr)
    print(
        f"  {reads_unmatched} reads had no matching source record or tags were absent ({pct_u:.1f}%)",
        file=sys.stderr,
    )
    print(f"Peak RAM: {ram:.0f} MB", file=sys.stderr)


if __name__ == "__main__":
    main()
