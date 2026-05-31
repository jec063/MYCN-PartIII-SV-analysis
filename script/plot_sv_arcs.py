from pathlib import Path
import argparse
import csv
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import matplotlib.colors as mcolors
import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_pos(value: str) -> int:
    """Return integer genomic position, stripping 'chr' prefix if needed."""
    return int(value.strip())


def arc_patch(ax, x1: float, x2: float, height: float,
              color: str, lw: float, alpha: float = 0.85,
              above: bool = True) -> None:
    """
    Draw a semi-ellipse arc between x1 and x2 on *ax*.

    The arc is drawn above (above=True) or below the x-axis.
    Uses a Bezier-style Path so it works without extra dependencies.
    """
    mid = (x1 + x2) / 2.0
    sign = 1 if above else -1

    # Build arc via parametric ellipse
    theta = np.linspace(0, math.pi, 200)
    xs = mid + (x2 - x1) / 2.0 * np.cos(theta)
    ys = sign * height * np.sin(theta)

    ax.plot(xs, ys, color=color, lw=lw, alpha=alpha, solid_capstyle="round")


def frequency_to_color(freq: int, max_freq: int) -> str:
    """Map frequency count to a colour on a blue→red gradient."""
    cmap = plt.get_cmap("RdYlBu_r")
    norm = mcolors.Normalize(vmin=1, vmax=max_freq)
    return cmap(norm(freq))


def frequency_to_lw(freq: int, max_freq: int,
                    lw_min: float = 1.0, lw_max: float = 4.0) -> float:
    """Map frequency count to line width."""
    if max_freq == 1:
        return lw_max
    return lw_min + (lw_max - lw_min) * (freq - 1) / (max_freq - 1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Draw genomic arc diagrams for top recurrent SVs."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input recurrent SV CSV file (output of count_frequency.py)."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output plot file, e.g. PNG or PDF."
    )
    parser.add_argument(
        "-t", "--top", type=int, default=15,
        help="Number of top SVs to visualise (default: 15)."
    )
    parser.add_argument(
        "--gene-start", type=int, default=15_940_550,
        help="MYCN gene start on chr2, hg19 (default: 15940550)."
    )
    parser.add_argument(
        "--gene-end", type=int, default=15_947_007,
        help="MYCN gene end on chr2, hg19 (default: 15947007)."
    )
    args = parser.parse_args()

    input_csv = Path(args.input)
    output_plot = Path(args.output)
    top_n = args.top
    output_plot.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load SVs
    # ------------------------------------------------------------------
    rows: list[dict] = []
    with open(input_csv, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)

    rows = rows[:top_n]

    # Normalise chromosome names  (some rows use '2', others 'chr2')
    def norm_chr(c: str) -> str:
        c = c.strip()
        if not c.startswith("chr"):
            c = "chr" + c
        return c

    svs = []
    for row in rows:
        svs.append({
            "sv_id":   row["SV_ID"],
            "chr1":    norm_chr(row["chr1"]),
            "pos1":    parse_pos(row["pos1"]),
            "strand1": row["strand1"],
            "chr2":    norm_chr(row["chr2"]),
            "pos2":    parse_pos(row["pos2"]),
            "strand2": row["strand2"],
            "sv_type": row["sv_type"],
            "freq":    int(row["frequency_count"]),
        })

    max_freq = max(sv["freq"] for sv in svs)

    # ------------------------------------------------------------------
    # Derive axis range  (only chr2 SVs in this dataset)
    # ------------------------------------------------------------------
    all_pos = []
    for sv in svs:
        if sv["chr1"] == sv["chr2"]:   # intrachromosomal
            all_pos += [sv["pos1"], sv["pos2"]]

    margin = 500_000
    x_min = min(all_pos) - margin
    x_max = max(all_pos) + margin

    # ------------------------------------------------------------------
    # Figure layout: two axes stacked
    #   top  : full-range view (all 15 SVs)
    #   bottom: zoomed inset on MYCN locus (SV_001–SV_005 region)
    # ------------------------------------------------------------------
    fig, (ax_full, ax_zoom) = plt.subplots(
        2, 1,
        figsize=(14, 8),
        gridspec_kw={"height_ratios": [3, 2]},
    )
    fig.subplots_adjust(hspace=0.45)

    # ── helpers shared by both panels ──────────────────────────────────
    def draw_chromosome_axis(ax, xlim_left, xlim_right, label="chr2"):
        ax.axhline(0, color="black", lw=1.5, zorder=1)
        ax.set_xlim(xlim_left, xlim_right)
        ax.set_xlabel("chr2 position (bp)", fontsize=9)
        ax.set_ylabel("Arc height\n(proportional to span)", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.yaxis.set_visible(False)
        # Chromosome label
        ax.text(
            xlim_left + (xlim_right - xlim_left) * 0.01, -0.05,
            label, fontsize=9, fontweight="bold",
            transform=ax.get_xaxis_transform(), va="top",
        )

    def draw_gene_box(ax, gene_start, gene_end, label="MYCN"):
        """Draw a filled rectangle below the axis to mark the gene."""
        ax.axvspan(gene_start, gene_end, ymin=0, ymax=0.04,
                   color="#e8453c", alpha=0.85, zorder=3)
        mid = (gene_start + gene_end) / 2
        # label just above the rectangle
        ax.annotate(
            label,
            xy=(mid, 0), xycoords=("data", "axes fraction"),
            xytext=(0, -18), textcoords="offset points",
            ha="center", fontsize=8, color="#e8453c", fontweight="bold",
        )

    def draw_svs_on_ax(ax, sv_list, x_left, x_right,
                       gene_start, gene_end, arc_scale=1.0):
        draw_chromosome_axis(ax, x_left, x_right)
        draw_gene_box(ax, gene_start, gene_end)

        max_span = max(
            abs(sv["pos2"] - sv["pos1"])
            for sv in sv_list
            if sv["chr1"] == sv["chr2"]
        )

        for sv in sv_list:
            if sv["chr1"] != sv["chr2"]:
                continue  # skip interchromosomal for linear view

            p1, p2 = sv["pos1"], sv["pos2"]
            if p1 > p2:
                p1, p2 = p2, p1

            # Skip if both endpoints outside view
            if p2 < x_left or p1 > x_right:
                continue

            span = p2 - p1
            height = arc_scale * (span / max_span)

            color = frequency_to_color(sv["freq"], max_freq)
            lw = frequency_to_lw(sv["freq"], max_freq)

            arc_patch(ax, p1, p2, height, color=color, lw=lw)

            # Tick marks at breakpoints
            for bp in (p1, p2):
                ax.plot([bp, bp], [0, 0.015], color=color, lw=1.2,
                        transform=ax.get_xaxis_transform(), clip_on=False)

            # Label with SV_ID near arc apex
            mid = (p1 + p2) / 2
            ax.text(
                mid, height + arc_scale * 0.03,
                sv["sv_id"],
                ha="center", va="bottom", fontsize=6.5,
                color=frequency_to_color(sv["freq"], max_freq),
            )

        ax.set_ylim(-0.15, arc_scale * 1.25)

    # ── Panel 1: full range ─────────────────────────────────────────────
    draw_svs_on_ax(
        ax_full, svs, x_min, x_max,
        args.gene_start, args.gene_end,
        arc_scale=1.0,
    )
    ax_full.set_title(
        f"Top {top_n} recurrent SVs in MYCN-amplified samples — chr2 overview",
        fontsize=11, fontweight="bold", pad=8,
    )

    # ── Panel 2: zoom on MYCN locus ─────────────────────────────────────
    zoom_margin = 200_000
    zoom_min = args.gene_start - zoom_margin
    zoom_max = 16_500_000          # covers SV_001–SV_005 cluster

    draw_svs_on_ax(
        ax_zoom, svs, zoom_min, zoom_max,
        args.gene_start, args.gene_end,
        arc_scale=0.8,
    )
    ax_zoom.set_title(
        "Zoomed view — MYCN locus (chr2 ~15.7–16.5 Mb)",
        fontsize=10, fontweight="bold", pad=8,
    )

    # ── Shared colour-bar legend ─────────────────────────────────────────
    cmap = plt.get_cmap("RdYlBu_r")
    norm = mcolors.Normalize(vmin=1, vmax=max_freq)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(
        sm, ax=[ax_full, ax_zoom],
        orientation="vertical",
        fraction=0.015, pad=0.02, shrink=0.8,
    )
    cbar.set_label("Frequency\n(# samples)", fontsize=9)
    cbar.set_ticks(range(1, max_freq + 1))

    # ── Frequency legend (line width) ───────────────────────────────────
    legend_handles = []
    for f in sorted(set(sv["freq"] for sv in svs), reverse=True):
        lw = frequency_to_lw(f, max_freq)
        color = frequency_to_color(f, max_freq)
        legend_handles.append(
            mpatches.Patch(
                facecolor=color, edgecolor="none",
                label=f"freq = {f}  (n={sum(1 for s in svs if s['freq']==f)} SVs)",
            )
        )
    ax_full.legend(
        handles=legend_handles,
        title="Recurrence", title_fontsize=8,
        fontsize=7.5, loc="upper right",
        framealpha=0.8,
    )

    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Wrote: {output_plot}")


main()