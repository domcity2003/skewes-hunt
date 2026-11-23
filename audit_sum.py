"""Recompute the Gaussian zero sum from consecutive checked LMFDB files."""

import argparse
from contextlib import contextmanager
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import subprocess
import tempfile
from time import perf_counter

import numpy as np
from mpmath import iv

from platt import manifest_hashes, read_fixed_blocks
from replay_smith import ALPHA, CENTER, CUTOFF
from tapered import error_budget


ARCHIVE = "https://beta.lmfdb.org/data/riemann-zeta-zeros"
SOURCE_COUNT = 103_800_788_359
PUBLISHED_SUM = "1.00288049783641096204"
CONTROL = {"center": "727.951335792", "alpha": "6e12", "cutoff": 10_379_600}
SMITH = {"center": CENTER, "alpha": ALPHA, "cutoff": CUTOFF}
CONTROL_INTERVAL = (1.002889273, 1.002889276)
PHASE_SCALE = 10**50
TAU_TICKS = 628318530717958647692528676655900576839433879875021
TAU = float(Decimal(TAU_TICKS) / PHASE_SCALE)
SIN = tuple((-1)**k / math.factorial(2*k+1) for k in range(13))
COS = tuple((-1)**k / math.factorial(2*k) for k in range(13))
EXP = tuple(1 / math.factorial(k) for k in range(18))


def block_score(start, end, high, low, settings):
    center = settings["center"]
    omega = float(center)
    alpha = float(settings["alpha"])
    if not (14 <= start < end <= CUTOFF and 0 < end - start <= 5000
            and 727 < omega < 728 and std_exponent(end, alpha) <= 10.5
            and 0 < len(high) <= 50_000 and end * (end-start) / alpha <= 0.009):
        raise ValueError("Block exceeds the domain of the arithmetic budget")

    offset = np.ldexp(high.astype(np.float64), -37) + np.ldexp(low.astype(np.float64), -101)
    if not (1e-10 < offset[0] and offset[-1] < end-start-1e-10
            and np.all(np.diff(offset) > 0)):
        raise ValueError("Unordered ordinates or invalid offsets")
    # Reduce the large phase exactly before adding the small block offset.
    numerator, denominator = start.as_integer_ratio()
    omega_ticks = int(Decimal(center) * PHASE_SCALE)
    modulus = denominator * TAU_TICKS
    remainder = (numerator * omega_ticks) % modulus / (denominator * PHASE_SCALE)
    phase = remainder + omega * offset
    phase -= np.rint(phase / TAU) * TAU
    if not np.max(np.abs(phase)) < 3.142:
        raise ArithmeticError("The reduced phase exceeds the polynomial domain")
    square = phase**2
    sine = phase * np.polynomial.polynomial.polyval(square, SIN)
    cosine = np.polynomial.polynomial.polyval(square, COS)

    # The exponential over a block is a scalar factor times a small correction.
    factor = 0.0
    exponent = -std_exponent(start, alpha) / 16
    for coefficient in reversed(EXP):
        factor = factor * exponent + coefficient
    for _ in range(4):
        factor *= factor
    change = -offset * (2*start + offset) / (2*alpha)
    weight = factor * (1 + change * (1 + change * (0.5 + change / 6)))

    gamma = start + offset
    norm = gamma**2 + 0.25
    real = (0.5 + (0.25 - gamma**2)/(omega * norm)) / norm
    imag = -gamma * (1 + 1/(omega * norm)) / norm
    terms = -2 * weight * (real * cosine - imag * sine)
    cutoff = settings["cutoff"]
    if start < cutoff < end:
        if np.any(np.abs(offset - (cutoff-start)) <= 1e-10):
            raise ValueError("An ordinate needs a higher-precision cutoff comparison")
        terms = terms[offset <= cutoff-start]
    return len(terms), float(np.sum(terms))


def std_exponent(height, alpha):
    return height**2 / (2 * alpha)


@contextmanager
def source_file(directory, name, download):
    local = directory / name
    if local.is_file():
        yield local, 0
    elif download:
        with tempfile.NamedTemporaryFile(dir=directory, prefix=name + ".", suffix=".part",
                                         delete=False) as scratch:
            path = Path(scratch.name)
        try:
            start = perf_counter()
            subprocess.run(["curl", "--silent", "--show-error", "--fail", "--location",
                            "--retry", "6", "--retry-all-errors", "--retry-max-time", "1200",
                            "--connect-timeout", "30", "--max-time", "300", "--cookie", "human=1",
                            "--output", str(path), f"{ARCHIVE}/{name}"], check=True)
            yield path, perf_counter() - start
        finally:
            path.unlink(missing_ok=True)
    else:
        raise FileNotFoundError(f"Missing {local}. Add the checked file or use --download.")


def report(settings, rows):
    iv.dps = 55
    b = iv.mpf
    if not abs(b(TAU_TICKS) / PHASE_SCALE - 2 * iv.pi) < b("1e-48"):
        raise ArithmeticError("Phase period failed its enclosure check")
    if not (abs(b(TAU) - 2 * iv.pi) < b("6e-16") and np.finfo(np.float64).nmant == 52):
        raise ArithmeticError("Unsupported binary64 environment")
    name, digest, height, count, _, _, _, _ = rows[-1]
    score = math.fsum(row[4] for row in rows)
    # Revers Lemma 4.8, with 2e-8/gamma computed-error budget per pair.
    error = b("2e-8") * (iv.log(b(height) / (2 * iv.pi))**2 / (4 * iv.pi) + b("0.9321"))
    # Round the displayed interval outwards rather than assuming decimal conversion is exact.
    error += b("1e-14")
    enclosed = b(score) + b([-1, 1]) * error
    result = {
        "settings": settings, "files": len(rows), "ending_file": name, "ending_md5": digest,
        "height": height, "zero_count": count, "score": score,
        "arithmetic_enclosure": str(enclosed),
        "bytes_checked": sum(row[5] for row in rows),
        "compute_seconds": sum(row[6] for row in rows),
        "download_seconds": sum(row[7] for row in rows),
    }
    if settings == CONTROL and height == CONTROL["cutoff"]:
        if count != 22_000_000 or not CONTROL_INTERVAL[0] < score < CONTROL_INTERVAL[1]:
            raise ArithmeticError(f"Control mismatch: {count:,} zeros, score {score}")
        result["control_passed"] = True
    if settings == SMITH:
        result["independent_margin_lower_bound"] = str(
            (enclosed.a - 1 - sum(error_budget(height).values())).a
        )
    if settings == SMITH and height == CUTOFF:
        if count != SOURCE_COUNT:
            raise ArithmeticError(f"Wrong number of zeros at Smith's cutoff: {count:,}")
        result["contains_published_value"] = bool(enclosed.a < b(PUBLISHED_SUM) < enclosed.b)
    return result


def run(args, settings):
    # The archive lists its files lexicographically. Use their numerical starting heights.
    hashes = manifest_hashes(args.manifest)
    catalog = sorted((int(name[6:-4]), name, digest) for name, digest in hashes.items()
                     if name.startswith("zeros_") and name.endswith(".dat"))
    manifest_id = hashlib.sha256(args.manifest.read_bytes()).hexdigest()
    metadata = json.dumps({"version": 1, "manifest_sha256": manifest_id, **settings}, sort_keys=True)
    with sqlite3.connect(args.checkpoint) as db:
        db.execute("CREATE TABLE IF NOT EXISTS config (value TEXT UNIQUE)")
        db.execute("CREATE TABLE IF NOT EXISTS results (name TEXT PRIMARY KEY, md5 TEXT, "
                   "height INTEGER, count INTEGER, score REAL, bytes INTEGER, "
                   "compute REAL, download REAL)")
        existing = db.execute("SELECT value FROM config").fetchone()
        if existing is None:
            db.execute("INSERT INTO config VALUES (?)", (metadata,))
            db.commit()
        elif existing != (metadata,):
            raise ValueError("Checkpoint parameters or manifest changed")
        rows = list(db.execute("SELECT * FROM results ORDER BY height"))
        if [(row[0], row[1]) for row in rows] != [(name, digest)
                                                for _, name, digest in catalog[:len(rows)]]:
            raise ValueError("Checkpoint files do not match the manifest prefix")
        height, count = rows[-1][2:4] if rows else (14, 0)
        if height > args.through:
            raise ValueError("Checkpoint already extends beyond the requested height")
        for slot, (start_height, name, digest) in enumerate(catalog[len(rows):], len(rows)):
            if height >= args.through:
                break
            if start_height != height:
                raise ValueError(f"Missing archive segment after {height:,}")
            with source_file(args.directory, name, args.download) as (path, download_seconds):
                started = perf_counter()
                scores = []
                processed = 0
                end_height = height
                for start, end, last, high, low in read_fixed_blocks(path, digest, height, count):
                    if start < settings["cutoff"]:
                        used, subtotal = block_score(start, end, high, low, settings)
                        scores.append(subtotal)
                        processed += used
                        end_height = min(int(end), settings["cutoff"])
                    if end <= settings["cutoff"] and count + processed != last:
                        raise ValueError(f"Wrong count before height {end}")
                expected_end = catalog[slot+1][0] if slot+1 < len(catalog) else CUTOFF
                if end_height != min(expected_end, settings["cutoff"]):
                    raise ValueError(f"Unexpected file endpoint for {name}")
                if end_height > args.through:
                    raise ValueError("--through must be a file endpoint or the mode's cutoff")
                height, count = end_height, count + processed
                row = (name, digest, height, count, math.fsum(scores), path.stat().st_size,
                       perf_counter() - started, download_seconds)
                db.execute("INSERT INTO results VALUES (?, ?, ?, ?, ?, ?, ?, ?)", row)
                db.commit()
                rows.append(row)
                if download_seconds or len(rows) % 10 == 0 or height == args.through:
                    print(f"{count:,} zeros through {height:,} ({len(rows):,} files)", flush=True)
        if height != args.through or not rows:
            raise ValueError(f"Reached {height:,}, expected {args.through:,}")
        result = report(settings, rows)
        print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--control", action="store_true")
    parser.add_argument("--through", type=int, help="file ending height, default the mode's cutoff")
    parser.add_argument("--download", action="store_true", help="stream missing public archive files")
    args = parser.parse_args()
    settings = CONTROL if args.control else SMITH
    args.through = settings["cutoff"] if args.through is None else args.through
    if not 14 < args.through <= settings["cutoff"]:
        parser.error("--through must lie within the mode's range")
    run(args, settings)


if __name__ == "__main__":
    main()
