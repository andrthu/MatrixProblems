#!/usr/bin/env python3
"""Plot transmissibility vs. |corresponding linear-system block entry|.

Same input CSV as plot_trans_vs_matrix.py (columns: i, j, trans, and the
NxN block entries m00, m01, ..., m<n-1><n-1>; N is detected from the CSV
header, 3 for the black-oil tool's blocks, 2 for the CO2 one's). Differs
from that script in two ways:
  - the block entry is taken in absolute value, so both axes are plain log
    scale (no symlog needed since there is no sign / zero-crossing left);
  - rows where the block entry is exactly zero are dropped before
    histogramming, per panel, since they carry no magnitude information
    and previously just saturated a horizontal bin near y=0.

Usage: plot_trans_vs_matrix_abs.py <input.csv> <output.png> [x-axis label]
"""
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import LogLocator

# Validated sequential "blue" ramp and chart chrome (light theme), see the
# dataviz skill's references/palette.md. Kept identical to
# plot_trans_vs_matrix.py so the two figures read as one family.
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"


def block_cols(field_names):
    """Block-entry column names ("m00", "m01", ...) present in the CSV, in
    row-major order, plus the block size N. N is detected as one more than
    the largest row/column digit seen across the CSV's "m<r><c>" columns,
    so the same script handles the 3x3 black-oil blocks and the 2x2 CO2
    ones without a hardcoded size."""
    matches = [re.fullmatch(r"m(\d)(\d)", name) for name in field_names]
    digits = [int(d) for m in matches if m for d in m.groups()]
    n = max(digits) + 1 if digits else 0
    return [f"m{r}{c}" for r in range(n) for c in range(n)], n


def make_cmap():
    cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE)
    cmap.set_bad(SURFACE)  # empty bins read as chart surface, not palest blue
    return cmap


def panel_histogram(x_edges, trans, y_abs, n=48):
    y_edges = np.geomspace(y_abs.min(), y_abs.max(), n)
    if y_edges[0] == y_edges[-1]:
        return None
    h, xe, ye = np.histogram2d(trans, y_abs, bins=[x_edges, y_edges])
    return np.ma.masked_equal(h.T, 0), xe, ye


def main():
    if len(sys.argv) < 3:
        print("usage: plot_trans_vs_matrix_abs.py <input.csv> <output.png> [x-axis label]", file=sys.stderr)
        return 1
    csv_path, png_path = sys.argv[1], sys.argv[2]
    x_label = sys.argv[3] if len(sys.argv) > 3 else "transmissibility"

    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    n_total = data.shape[0]
    if n_total == 0:
        print("no data rows in " + csv_path, file=sys.stderr)
        return 1

    trans_all = data["trans"]
    positive = trans_all > 0
    n_dropped_trans = int(n_total - positive.sum())
    if n_dropped_trans:
        print(f"warning: dropped {n_dropped_trans} pairs with {x_label} <= 0 (log x-axis)", file=sys.stderr)
    trans_all = trans_all[positive]

    x_edges = np.geomspace(trans_all.min(), trans_all.max(), 70)

    plt.rcParams.update({
        "font.family": ["system-ui", "DejaVu Sans", "sans-serif"],
        "font.size": 9,
        "text.color": INK_PRIMARY,
        "axes.edgecolor": GRID,
        "axes.labelcolor": INK_SECONDARY,
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    })

    cmap = make_cmap()

    block_columns, n = block_cols(data.dtype.names)

    # First pass: compute all N*N histograms (each on its own nonzero
    # subset) and a shared vmax so the panels are comparable under one
    # colorbar.
    panels = []
    for col in block_columns:
        y = data[col][positive]
        nz = y != 0
        y_abs = np.abs(y[nz])
        trans = trans_all[nz]
        n_zero = int((~nz).sum())
        panels.append((col, n_zero, panel_histogram(x_edges, trans, y_abs)))

    vmax = max((p[2][0].max() for p in panels if p[2] is not None), default=1)

    # Header (title + subtitle) and footer (caption) need roughly the same
    # absolute vertical space (fixed font sizes) no matter how tall the
    # figure is, so reserve them in inches and convert to figure fractions
    # -- fixed fractions like the panel grid uses would leave a shrinking,
    # eventually overlapping, header/footer as N (and so the figure height)
    # shrinks from 3 (black-oil) to 2 (CO2).
    fig_h = 4.0 * n + 1.7
    header_in, footer_in = 1.0, 0.85
    top = 1 - header_in / fig_h
    bottom = footer_in / fig_h

    fig, axes = plt.subplots(n, n, figsize=(4.4 * n, fig_h), sharex=True, facecolor=SURFACE)
    axes = np.atleast_2d(axes)
    fig.suptitle(f"{x_label} vs. |linear-system block entries|",
                  color=INK_PRIMARY, fontsize=15, y=1 - 0.28 / fig_h)
    subtitle = f"{trans_all.size:,} off-diagonal blocks, zero-valued entries dropped per panel"
    if n_dropped_trans:
        subtitle += f" ({n_dropped_trans:,} more dropped, {x_label} <= 0)"
    fig.text(0.5, 1 - 0.62 / fig_h, subtitle, ha="center", color=INK_SECONDARY, fontsize=10)

    x_locator = LogLocator(base=10, numticks=6)
    y_locator = LogLocator(base=10, numticks=6)

    mappable = None
    for idx, (col, n_zero, hist) in enumerate(panels):
        r, c = divmod(idx, n)
        ax = axes[r, c]
        ax.set_facecolor(SURFACE)
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(x_locator)

        if hist is None:
            ax.text(0.5, 0.5, "all values equal\n(constant / zero)",
                     ha="center", va="center", color=INK_MUTED, fontsize=9,
                     transform=ax.transAxes)
        else:
            h, xe, ye = hist
            mappable = ax.pcolormesh(xe, ye, h, cmap=cmap, norm=LogNorm(vmin=1, vmax=vmax),
                                      shading="flat")
            ax.set_yscale("log")
            ax.yaxis.set_major_locator(y_locator)

        ax.grid(True, color=GRID, linewidth=0.6, alpha=0.7)
        title = f"row {r}, col {c}"
        if n_zero:
            title += f"  ({n_zero:,} zeros dropped)"
        ax.set_title(title, color=INK_SECONDARY, fontsize=10, pad=6)
        if r == n - 1:
            ax.set_xlabel(x_label)
        if c == 0:
            ax.set_ylabel("|matrix entry|")

    fig.subplots_adjust(left=0.07, right=0.89, top=top, bottom=bottom,
                         wspace=0.45, hspace=0.35)
    fig.text(0.5, 0.15 / fig_h,
              f"row = equation index, column = unknown index within each {n}×{n} off-diagonal block",
              ha="center", color=INK_MUTED, fontsize=9)

    if mappable is not None:
        cbar_ax = fig.add_axes([0.91, bottom, 0.02, top - bottom])
        cbar = fig.colorbar(mappable, cax=cbar_ax, label="point density (count)")
        cbar.ax.yaxis.label.set_color(INK_SECONDARY)
        cbar.ax.tick_params(colors=INK_MUTED, labelsize=8)

    fig.savefig(png_path, dpi=160, facecolor=SURFACE)
    print(f"pairs plotted (nonzero, per panel varies); trans dropped: {n_dropped_trans}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
