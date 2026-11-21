"""Render the research figures from the published intervals and taper parameters."""

from decimal import Decimal
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure
from scipy.stats import irwinhall

from tapered import (
    ALPHA, CENTER, DERIVATIVE_ORDER, LEFT_WIDTH, RIGHT_WIDTH, TAPER_WIDTH,
)


FIGURES = Path(__file__).resolve().parent / "figures"
INK = "#162d3d"
MUTED = "#536776"
BLUE = "#126786"
GRID = "#dbe4ea"


def crossover_intervals():
    center = Decimal(CENTER)
    scale = Decimal("1e9")
    smith = [float((Decimal(endpoint) - center) * scale)
             for endpoint in ("727.9513353997", "727.9513354023")]
    tapered = [-float(Decimal(LEFT_WIDTH) * scale), float(Decimal(RIGHT_WIDTH) * scale)]
    figure = Figure(figsize=(8.4, 3.05), layout="constrained")
    ax = figure.subplots()
    ax.set_title("Where a crossover is guaranteed to occur", loc="left", pad=14)
    for y, interval, color in ((1, smith, MUTED), (0, tapered, BLUE)):
        left, right = interval
        ax.hlines(y, left, right, linewidth=3.5, color=color, zorder=3)
        ax.plot([left, right], [y, y], "o", markersize=6, color=color, zorder=4)
        ax.axvline(right, color=color, linewidth=1, linestyle=":", alpha=0.65)
    ax.set(ylim=(-0.58, 1.43), xlim=(-1.83, 1.53), yticks=[1, 0],
           yticklabels=["Smith (2016)\nAssuming RH", "Tapered Gaussian\nFinite verification"],
           xlabel=r"Offset from $\log x = 727.951335401$ (billionths)")
    ax.text(1.19, 0.49, "0.22 earlier", ha="right", va="center", color=BLUE,
            fontsize=9.5, bbox={"facecolor": "white", "edgecolor": "none", "pad": 2})
    ax.grid(axis="x", color=GRID)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    return figure


def taper_profile():
    scale = 1e9
    left, right, width = (float(value) * scale
                          for value in (LEFT_WIDTH, RIGHT_WIDTH, TAPER_WIDTH))
    s = np.linspace(-1.7, 1.19, 1200)
    boxes = DERIVATIVE_ORDER + 1
    # U = width * (Y / boxes - 1/2) for an Irwin–Hall variable Y.
    # The transition regions are disjoint, so at least one CDF is one.
    lower_tail = irwinhall.cdf(np.clip(boxes * (right - s) / width, 0, boxes), boxes)
    upper_tail = irwinhall.cdf(np.clip(boxes * (left + s) / width, 0, boxes), boxes)
    taper = np.minimum(lower_tail, upper_tail)
    gaussian = np.exp(-0.5 * float(ALPHA) * (s / scale)**2)
    floor = 1e-17
    weight = gaussian * taper

    figure = Figure(figsize=(8.4, 3.4), layout="constrained")
    ax = figure.subplots()
    ax.set_title("The taper preserves the center and removes the tails", loc="left", pad=14)
    ax.plot(s, gaussian, color=MUTED, linestyle="--", linewidth=1.6, label="Gaussian")
    ax.plot(s, np.where(weight >= floor, weight, np.nan), color=BLUE, linewidth=2.5,
            label="Tapered Gaussian")
    for endpoint in (-left, right):
        ax.axvline(endpoint, color=MUTED, linewidth=1, linestyle=":", alpha=0.7)
    ax.set(yscale="log", ylim=(floor, 2), xlim=(s[0], s[-1]),
           xlabel=r"Offset from the center $\omega$ (billionths)",
           ylabel="Weight relative to peak")
    ax.set_yticks([1e-16, 1e-12, 1e-8, 1e-4, 1])
    ax.grid(axis="y", color=GRID)
    ax.legend(frameon=False, loc="lower center")
    ax.spines[["top", "right"]].set_visible(False)
    return figure


def main():
    FIGURES.mkdir(exist_ok=True)
    with mpl.rc_context({
        "font.family": "DejaVu Sans", "font.size": 10.5,
        "mathtext.fontset": "stix", "axes.titlesize": 14,
        "axes.titleweight": "bold", "axes.edgecolor": GRID,
        "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
    }):
        crossover_intervals().savefig(FIGURES / "crossover-intervals.png", dpi=180)
        taper_profile().savefig(FIGURES / "taper-profile.png", dpi=180)


if __name__ == "__main__":
    main()
