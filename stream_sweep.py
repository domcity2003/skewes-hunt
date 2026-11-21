"""Explore the earlier crossover with resumable sums over LMFDB files."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from earlier import ALPHAS, CONTROL, ENDPOINT, TARGET, accumulate
from platt import read_blocks, source_files
from search import error_terms


CONTROL_HEIGHT = 10_379_600
CONTROL_ALPHA = 6e12
CONTROL_SUM = 1.002889274
VERIFIED_HEIGHT = 3e12
CENTERS = np.r_[TARGET + np.linspace(-3e-7, 3e-7, 25), CONTROL]
WIDTHS = np.r_[ENDPOINT - CENTERS[:-1] - 1e-12, 4e-6]
SMOOTHING = np.unique(np.r_[ALPHAS, np.geomspace(4e13, 1e15, 25)])


def write_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def summarize(total, cutoff):
    rows = []
    for i, (center, eta) in enumerate(zip(CENTERS[:-1], WIDTHS[:-1])):
        for j, alpha in enumerate(SMOOTHING):
            errors = error_terms(center, eta, alpha, cutoff, VERIFIED_HEIGHT)
            raw = float(total[i, j])
            base_margin = raw - 1 - errors["R1"]
            rh_margin = base_margin - math.fsum(errors[name] for name in ("R2", "R3", "R4", "R5"))
            rows.append({
                "center": float(center), "eta": float(eta), "alpha": float(alpha), "sum": raw,
                "prime_side_margin": base_margin, "errors": errors,
                "margin_full_rh": rh_margin, "margin_finite_height": rh_margin - errors["R6"],
            })
    return {
        "highest_raw_signal": max(rows, key=lambda row: row["sum"]),
        "best_full_rh": max(rows, key=lambda row: row["margin_full_rh"]),
        "best_finite_height": max(rows, key=lambda row: row["margin_finite_height"]),
    }


def control_result(total, count):
    score = float(total[-1, np.searchsorted(SMOOTHING, CONTROL_ALPHA)])
    errors = error_terms(CONTROL, float(WIDTHS[-1]), CONTROL_ALPHA, CONTROL_HEIGHT)
    margin = score - 1 - math.fsum(errors.values())
    if count != 22_000_000 or abs(score - CONTROL_SUM) > 5e-8 or margin <= 0:
        raise ArithmeticError(f"Failed the Teal Sea control: count={count}, sum={score}, margin={margin}")
    return {"cutoff": CONTROL_HEIGHT, "zeros": count, "sum": score, "margin_full_rh": margin}


def run(directory, manifest, checkpoint, output, through):
    if checkpoint.resolve() == output.resolve():
        raise ValueError("Checkpoint and report need distinct paths")
    files = source_files(directory, manifest)
    files = [(path, digest) for path, digest in files if int(path.stem[6:]) < through]
    configuration = {
        "version": 1, "arithmetic": "float64 screening only",
        "centers": CENTERS.tolist(), "eta": WIDTHS.tolist(), "alpha": SMOOTHING.tolist(),
        "verified_height": VERIFIED_HEIGHT,
    }
    report = {"configuration": configuration, "control": None, "results": []}
    total = np.zeros((len(CENTERS), len(SMOOTHING)))
    height, index = 14, 0
    if checkpoint.exists():
        saved = json.loads(checkpoint.read_text(encoding="utf-8"))
        report = saved["report"]
        if report["configuration"] != configuration:
            raise ValueError("Checkpoint uses a different search configuration")
        results = report["results"]
        sources = [(row["file"], row["md5"]) for row in results]
        if sources != [(path.name, digest) for path, digest in files[:len(sources)]]:
            raise ValueError("Checkpoint does not match the selected input prefix")
        height, index = results[-1]["cutoff"], results[-1]["zeros"]
        total = np.asarray(saved["sums"], dtype=np.float64)
        if total.shape != (len(CENTERS), len(SMOOTHING)) or not np.all(np.isfinite(total)):
            raise ValueError("Invalid sums in checkpoint")

    for path, digest in files[len(report["results"]):]:
        pending = []
        for end, last, ordinates in read_blocks(path, digest, height, index):
            if end >= CONTROL_HEIGHT > height:
                split = np.searchsorted(ordinates, CONTROL_HEIGHT, side="right")
                pending.append(ordinates[:split])
                accumulate(total, np.concatenate(pending), CENTERS, SMOOTHING)
                report["control"] = control_result(total, index + int(split))
                pending = []
                ordinates = ordinates[split:]
            pending.append(ordinates)
            height, index = end, last
        accumulate(total, np.concatenate(pending), CENTERS, SMOOTHING)
        row = {"file": path.name, "md5": digest, "cutoff": height, "zeros": index}
        # Tiny seed files only initialize the sums. Begin comparisons at 446,000.
        if height >= 446_000:
            row.update(summarize(total, height))
        report["results"].append(row)
        write_json(checkpoint, {"report": report, "sums": total.tolist()})
        write_json(output, report)
        print(f"Height {height:,.0f}, {index:,} zeros", flush=True)

    if height != through:
        raise ValueError(f"Input files reached {height:,.0f}, expected {through:,}")
    write_json(output, report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--through", type=int, required=True, help="ending height of the last file")
    args = parser.parse_args()
    run(args.directory, args.manifest, args.checkpoint, args.output, args.through)


if __name__ == "__main__":
    main()
