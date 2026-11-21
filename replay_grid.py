"""Apply Revers's full-RH error budget to the published Teal Sea grid."""

import argparse
from decimal import Decimal
from fractions import Fraction
import json
import math
from pathlib import Path
import re

from search import error_terms


SMITH_ENDPOINT = Decimal("727.9513354023")
ALPHA = "6e12"
SHARD_COUNT = 220
SHARD_SIZE = 100_000
CUTOFF = 10_379_600
ABALL = re.compile(r"\[([^ ]+) \+/- ([^ ]+)\]")
FOURTH_DERIVATIVE_BOUND = 9 * 10**26


def endpoints(ball):
    match = ABALL.fullmatch(ball)
    if match is None:
        raise ValueError(f"Unexpected upstream enclosure: {ball!r}")
    center, radius = map(Fraction, match.groups())
    if radius < 0:
        raise ValueError(f"Negative upstream radius: {ball!r}")
    return center - radius, center + radius


def cubic_weights(at):
    values, derivatives = [], []
    for node in range(4):
        others = [index for index in range(4) if index != node]
        denominator = math.prod(node - index for index in others)
        values.append(Fraction(math.prod(at - index for index in others), denominator))
        derivatives.append(Fraction(sum(math.prod(at - index for index in others if index != skipped)
                                        for skipped in others), denominator))
    return values, derivatives


def continuous_exclusion(grid_config, sums):
    """Bound the sum between uniform grid centers using cubic interpolation."""
    grid = grid_config["descriptive_grid"]
    step = Fraction(grid["step"])
    center = Fraction(grid["center"])
    first, last = grid["k_range"]
    nodes = [center + k * step for k in range(first, last + 1)]
    sums_by_node = {Fraction(key): value for key, value in sums.items()}
    intervals = [sums_by_node[node] for node in nodes]
    if not 700 < nodes[0] < Fraction(SMITH_ENDPOINT) < nodes[-1]:
        raise ValueError("The grid does not cover the requested range")

    # For four equally spaced nodes, the product in the interpolation remainder
    # has absolute value at most step**4 anywhere between the outer nodes.
    remainder = FOURTH_DERIVATIVE_BOUND * step**4 / 24
    ceiling = None
    for segment, (low, high) in enumerate(zip(nodes, nodes[1:])):
        if low >= Fraction(SMITH_ENDPOINT):
            break
        start = min(max(0, segment - 1), len(nodes) - 4)
        values = intervals[start:start + 4]
        at_left, slope_left = cubic_weights(segment - start)
        at_right, slope_right = cubic_weights(segment - start + 1)
        # A cubic's Bernstein basis on the segment is nonnegative and sums to 1.
        controls = (at_left, [a + b / 3 for a, b in zip(at_left, slope_left)],
                    [a - b / 3 for a, b in zip(at_right, slope_right)], at_right)
        polynomial_floor = min(sum(weight * values[index][weight < 0]
                                   for index, weight in enumerate(weights)) for weights in controls)
        covered_high = min(high, Fraction(SMITH_ENDPOINT))
        bound = -1 - polynomial_floor + remainder - 2 / covered_high
        ceiling = bound if ceiling is None else max(ceiling, bound)
    return {"lower_log_x": str(Decimal(nodes[0].numerator) / Decimal(nodes[0].denominator)),
            "upper_log_x": str(SMITH_ENDPOINT),
            "upper_bound_on_best_possible_margin": str(ceiling),
            "upper_bound_approx": float(ceiling),
            "rules_out_positive_lower_bound_at_any_eta": ceiling < 0}


def replay(upstream):
    shard_dir = upstream / "shards_b2"
    published = json.loads((shard_dir / "evaluation.json").read_text(encoding="utf-8"))
    grid_config = json.loads((upstream / "omega_grid.json").read_text(encoding="utf-8"))
    centers = set(grid_config["omegas"])
    shard_paths = sorted(shard_dir.glob("shard_*.json"))
    if len(shard_paths) != SHARD_COUNT or not published["checks"]["complete"]:
        raise ArithmeticError("Incomplete upstream grid")
    sums = {center: (Fraction(0), Fraction(0)) for center in centers}
    last_upper = 0
    for index, path in enumerate(shard_paths):
        shard = json.loads(path.read_text(encoding="utf-8"))
        first_lower, _ = endpoints(shard["first"])
        _, next_upper = endpoints(shard["last"])
        if (shard["n0"] != index * SHARD_SIZE + 1 or shard["num"] != SHARD_SIZE
                or not shard["ordered_disjoint"] or not last_upper < first_lower < next_upper
                or shard["grid_alpha"] != ALPHA or set(shard["grid_sums"]) != centers):
            raise ArithmeticError(f"Inconsistent upstream shard: {path.name}")
        last_upper = next_upper
        for center, ball in shard["grid_sums"].items():
            left, right = endpoints(ball)
            total_left, total_right = sums[center]
            sums[center] = total_left + left, total_right + right
    if last_upper >= CUTOFF or published["checks"]["turing_count_N(T)"] != SHARD_COUNT * SHARD_SIZE:
        raise ArithmeticError("Upstream cutoffs do not agree")
    eta = "0.000004"
    rows = []
    for center, (left, right) in sums.items():
        published_left, published_right = endpoints(published["descriptive_grid"][center])
        if right < published_left or published_right < left:
            raise ArithmeticError(f"Recombined sum disagrees with the upstream grid at {center}")
        errors = error_terms(float(center), float(eta), float(ALPHA), CUTOFF)
        margin = -1 - float(right) - sum(errors.values())
        # For any positive eta, the first summand of R1 is larger than 2/omega.
        ceiling = -1 - left - 2 / Fraction(center)
        endpoint = str(Decimal(center) + Decimal(eta))
        rows.append({"omega": center, "zero_sum_interval_approx": [float(left), float(right)],
                     "eta": eta, "upper_log_x": endpoint, "estimated_margin_full_rh": margin,
                     "rules_out_positive_lower_bound_at_any_eta": ceiling < 0})
    rows.sort(key=lambda row: Decimal(row["omega"]))
    before_record = [row for row in rows if Decimal(row["omega"]) < SMITH_ENDPOINT]
    positive = [row for row in rows if row["estimated_margin_full_rh"] > 0]
    return {"assumption": "full Riemann hypothesis", "smith_upper_log_x": str(SMITH_ENDPOINT),
            "arithmetic": "exact for shard recombination and eta-independent exclusion; float64 for margins",
            "turing_count_source": "published Teal Sea report, not recomputed",
            "grid_centers": len(rows), "recombined_zero_count": SHARD_COUNT * SHARD_SIZE,
            "centers_before_smith": len(before_record),
            "centers_before_smith_ruled_out_at_any_eta": sum(
                row["rules_out_positive_lower_bound_at_any_eta"] for row in before_record),
            "continuous_exclusion": continuous_exclusion(grid_config, sums),
            "earliest_estimated_positive_at_fixed_eta": positive[0] if positive else None,
            "grid": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("upstream", type=Path, help="Path to the upstream hunts/skewes_enclosure directory")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = replay(args.upstream)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
