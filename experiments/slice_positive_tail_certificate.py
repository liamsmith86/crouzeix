#!/usr/bin/env python3
"""Certify the remaining positive-sign compact nome tail.

The two determinant charts use the full deficit-centered certificates from
``slice_positive_face_certificate.py``; bounded coarse re-centering is used
only when binary64 subdivision approaches the exact corner square.  The
three minor charts use the original directed checker.  Exact rational
geometric boxes cover the requested interval without gaps.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import gc
import time

from slice_positive_face_certificate import (
    audit_bernstein_machinery,
    audit_recenter_machinery,
    face_slab_result,
)
from slice_positive_corner_audit import audit_corner_tables
from slice_positive_face_audit import audit_table
from slice_projective_core import load_records
from slice_projective_interval_certificate import (
    certify,
    final_chart_tables,
    interval_tensor_for_chart,
    parse_fraction,
)


DEFAULT_START = Fraction(
    7471344308886696360166308338925321050191281742307131164807,
    13421772800000000000000000000000000000000000000000000000000,
)


def geometric_cover(
    start: Fraction, stop: Fraction, ratio: Fraction
) -> list[tuple[Fraction, Fraction]]:
    """Return a gap-free exact geometric cover of one rational interval."""

    if not 0 < start < stop < 1:
        raise ValueError("require 0 < start < stop < 1")
    if ratio <= 1:
        raise ValueError("geometric ratio must exceed one")
    output: list[tuple[Fraction, Fraction]] = []
    lower = start
    while lower < stop:
        upper = min(stop, lower * ratio)
        output.append((lower, upper))
        lower = upper
    return output


def audit_cover(
    boxes: list[tuple[Fraction, Fraction]], start: Fraction, stop: Fraction
) -> None:
    """Check exact endpoint coverage and adjacency of a generated cover."""

    if not boxes or boxes[0][0] != start or boxes[-1][1] != stop:
        raise AssertionError("geometric cover does not have the requested endpoints")
    if any(left[1] != right[0] for left, right in zip(boxes, boxes[1:])):
        raise AssertionError("geometric cover contains a gap or overlap")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=parse_fraction, default=DEFAULT_START)
    parser.add_argument("--stop", type=parse_fraction, default=Fraction(63, 100))
    parser.add_argument("--ratio", type=parse_fraction, default=Fraction(81, 80))
    parser.add_argument("--degree", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=40)
    parser.add_argument("--max-leaves", type=int, default=100_000)
    parser.add_argument("--recenter-depth", type=int, default=1)
    parser.add_argument("--max-recentered-boxes", type=int, default=64)
    parser.add_argument(
        "--regenerate-records",
        action="store_true",
        help="regenerate the exact transfer-core tables before certification",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.degree < 1:
        raise SystemExit("require degree >= 1")
    if args.max_depth < 0 or args.max_leaves < 1:
        raise SystemExit("require nonnegative depth and positive leaf limit")
    if args.recenter_depth < 0 or args.max_recentered_boxes < 0:
        raise SystemExit("require nonnegative recenter limits")
    records = load_records(force=args.regenerate_records)
    tables = final_chart_tables(records, 1)
    for chart, table in enumerate(tables[:2]):
        audit_table(table, chart)
    audit_corner_tables(records)
    audit_bernstein_machinery()
    audit_recenter_machinery()
    boxes = geometric_cover(args.start, args.stop, args.ratio)
    audit_cover(boxes, args.start, args.stop)
    started = time.monotonic()

    for box_index, box in enumerate(boxes, start=1):
        box_started = time.monotonic()
        print(
            f"BOX {box_index}/{len(boxes)} {box[0]} {box[1]}",
            flush=True,
        )
        for chart, table in enumerate(tables[:2]):
            chart_started = time.monotonic()
            slab = face_slab_result(
                table,
                chart,
                box,
                args.degree,
                Fraction(0),
                Fraction(0),
                Fraction(1),
                args.max_depth,
                args.max_leaves,
                recenter_depth=args.recenter_depth,
                max_recentered_boxes=args.max_recentered_boxes,
            )
            print(
                f"  {table.label}: slab_root={slab.centered_root_lower:.3e}, "
                f"slab_leaves={slab.centered_certificate.leaves}, "
                f"recentered={len(slab.recentered_boxes)}, "
                f"depth={slab.centered_certificate.depth}, "
                f"time={time.monotonic() - chart_started:.2f}s",
                flush=True,
            )
            gc.collect()

        for table in tables[2:]:
            chart_started = time.monotonic()
            tensor = interval_tensor_for_chart(table, box, args.degree)
            result = certify(
                tensor,
                max_depth=args.max_depth,
                max_leaves=args.max_leaves,
            )
            print(
                f"  {table.label}: pass={result.passed}, "
                f"leaves={result.leaves}, depth={result.depth}, "
                f"time={time.monotonic() - chart_started:.2f}s",
                flush=True,
            )
            if not result.passed:
                print(
                    f"    failure_box={result.failure_box}, "
                    f"index={result.failure_index}",
                    flush=True,
                )
                raise SystemExit(1)
            gc.collect()
        print(
            f"BOX PASS {box_index}/{len(boxes)} "
            f"time={time.monotonic() - box_started:.2f}s",
            flush=True,
        )

    print(
        f"PASS positive tail: boxes={len(boxes)}, start={args.start}, "
        f"stop={args.stop}, elapsed={time.monotonic() - started:.2f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
