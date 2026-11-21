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
BLUE = "#006a60"
GRID = "#d9e3df"
PAPER = "#f5f8f6"


def crossover_intervals():
    center = Decimal(CENTER)
    scale = Decimal("1e9")
    smith = [float((Decimal(endpoint) - center) * scale)
             for endpoint in ("727.9513353997", "727.9513354023")]
    tapered = [-float(Decimal(LEFT_WIDTH) * scale), float(Decimal(RIGHT_WIDTH) * scale)]
    background = "#101f24"
    white = "#f1f7f3"
    previous = "#9fb1b8"
    accent = "#80ecc6"
    figure = Figure(figsize=(10, 3.55), facecolor=background)
    figure.text(0.06, 0.875, "An earlier upper endpoint", fontfamily="DejaVu Serif",
                fontsize=24, color=white, va="center")
    ax = figure.add_axes((0.06, 0.24, 0.88, 0.50), facecolor=background)
    for y, interval, color in ((1, smith, previous), (0, tapered, accent)):
        left, right = interval
        ax.hlines(y, left, right, linewidth=3, color=color, zorder=3)
        ax.plot(left, y, "o", markersize=5, color=color, zorder=4)
        ax.plot(right, y, "o", markersize=8, markeredgewidth=1.8,
                markerfacecolor=background, markeredgecolor=color, zorder=4)
    ax.text(-1.70, 1.25, "Smith · 2016", fontsize=12, color=previous)
    ax.text(-1.70, 0.25, "Tapered Gaussian · candidate", fontsize=12, color=accent)
    ax.annotate("", xy=(tapered[1], 0.51), xytext=(smith[1], 0.51),
                arrowprops={"arrowstyle": "->", "color": accent, "lw": 1.4,
                            "shrinkA": 0, "shrinkB": 0})
    for right, y, color in ((tapered[1], 0, accent), (smith[1], 1, previous)):
        ax.vlines(right, min(0.51, y), max(0.51, y), color=color,
                  linewidth=0.8, linestyle=":", alpha=0.8)
    ax.text(tapered[1] - 0.07, 0.51, "0.22 billionths earlier", ha="right",
            va="center", color=accent, fontsize=10.5,
            bbox={"facecolor": background, "edgecolor": "none", "pad": 2})
    ax.set(ylim=(-0.43, 1.51), xlim=(-1.73, 1.49), yticks=[],
           xticks=[-1.5, -1, -0.5, 0, 0.5, 1])
    ax.set_xlabel(r"Offset from $\log x = 727.951335401$ (billionths)",
                  color=previous, labelpad=12)
    ax.grid(axis="x", color="#32474e", alpha=0.7)
    ax.tick_params(axis="x", colors=previous, length=0, pad=6)
    ax.spines[:].set_visible(False)
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

    figure = Figure(figsize=(8.4, 3.4), layout="constrained", facecolor=PAPER)
    figure.get_layout_engine().set(w_pad=0.16, h_pad=0.14)
    ax = figure.subplots()
    ax.set_facecolor(PAPER)
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
