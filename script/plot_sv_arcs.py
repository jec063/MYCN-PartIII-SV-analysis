from pathlib import Path
import argparse
import csv
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_pos(value: str) -> int:
    return int(value.strip())


def arc_patch(ax, x1: float, x2: float, height: float,
              color, lw: float, alpha: float = 0.9) -> None:
    """Draw a filled semi-ellipse arc between x1 and x2."""
    mid = (x1 + x2) / 2.0
    theta = np.linspace(0, math.pi, 300)
    xs = mid + (x2 - x1) / 2.0 * np.cos(theta)
    ys = height * np.sin(theta)
    ax.plot(xs, ys, color=color, lw=lw, alpha=alpha, solid_capstyle="round")


def frequency_to_color(freq: int, max_freq: int):
    # Distinct, saturated colours per frequency tier:
    # freq=2 → deep blue, freq=3 → vivid orange, freq=4 → crimson red
    palette = {
        1: "#4a90d9",
        2: "#2166ac",
        3: "#f4a82e",
        4: "#d62728",
    }
    if freq in palette:
        return palette[freq]
    # fallback gradient for any frequency outside the palette
    cmap = plt.get_cmap("RdYlBu_r")
    norm = mcolors.Normalize(vmin=1, vmax=max_freq)
    return cmap(norm(freq))


def frequency_to_lw(freq: int, max_freq: int,
                    lw_min: float = 1.5, lw_max: float = 5.0) -> float:
    if max_freq == 1:
        return lw_max
    return lw_min + (lw_max - lw_min) * (freq - 1) / (max_freq - 1)


def draw_panel(ax, sv_list, x_min, x_max,
               gene_start, gene_end, gene_label,
               title, max_freq, label_svs=True):
    """Draw one panel: chromosome axis + arcs for SVs whose breakpoints fall in range."""

    # Filter SVs: at least one breakpoint inside the view
    visible = []
    for sv in sv_list:
        if sv["chr1"] != sv["chr2"]:
            continue
        p1, p2 = sorted([sv["pos1"], sv["pos2"]])
        # include if either endpoint is in range, or arc spans the range
        if p2 >= x_min and p1 <= x_max:
            visible.append((p1, p2, sv))

    if not visible:
        ax.set_visible(False)
        return

    max_span = max(p2 - p1 for p1, p2, _ in visible) if visible else 1

    for p1, p2, sv in visible:
        span = p2 - p1
        height = span / max_span  # normalised 0–1
        color = frequency_to_color(sv["freq"], max_freq)
        lw = frequency_to_lw(sv["freq"], max_freq)
        arc_patch(ax, p1, p2, height, color=color, lw=lw)

        # Tick marks at breakpoints
        for bp in (p1, p2):
            if x_min <= bp <= x_max:
                ax.plot([bp, bp], [0, 0.02], color=color, lw=1.5,
                        transform=ax.get_xaxis_transform(), clip_on=False)

        if label_svs:
            mid = (p1 + p2) / 2.0
            ax.text(mid, height + 0.04, sv["sv_id"],
                    ha="center", va="bottom", fontsize=7,
                    color=color, fontweight="bold")

    # Gene annotation box
    if gene_end >= x_min and gene_start <= x_max:
        ax.axvspan(gene_start, gene_end, ymin=0, ymax=0.05,
                   color="#e8453c", alpha=0.9, zorder=4)
        mid_gene = (gene_start + gene_end) / 2
        ax.annotate(
            gene_label,
            xy=(mid_gene, 0),
            xycoords=("data", "axes fraction"),
            xytext=(0, -20), textcoords="offset points",
            ha="center", fontsize=8.5, color="#e8453c", fontweight="bold",
        )

    # Axis styling
    ax.axhline(0, color="black", lw=2, zorder=2)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-0.18, 1.35)
    ax.set_title(title, fontsize=10, fontweight="bold", pad=6)
    ax.set_xlabel("chr2 position (bp)", fontsize=8.5)
    ax.yaxis.set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Chromosome label bottom-left
    ax.text(0.01, -0.08, "chr2", transform=ax.transAxes,
            fontsize=8, fontweight="bold", va="top")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Draw genomic arc diagrams for top recurrent SVs."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input recurrent SV CSV file.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output plot file, e.g. PNG or PDF.")
    parser.add_argument("-t", "--top", type=int, default=15,
                        help="Number of top SVs to visualise (default: 15).")
    parser.add_argument("--gene-start", type=int, default=16_225_000,
                        help="MYCN gene start on chr2 in your data.")
    parser.add_argument("--gene-end", type=int, default=16_227_000,
                        help="MYCN gene end on chr2 in your data.")
    args = parser.parse_args()

    input_csv  = Path(args.input)
    output_plot = Path(args.output)
    top_n = args.top
    output_plot.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load SVs
    # ------------------------------------------------------------------
    def norm_chr(c: str) -> str:
        c = c.strip()
        return c if c.startswith("chr") else "chr" + c

    svs = []
    with open(input_csv, "r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
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

    svs = svs[:top_n]
    max_freq = max(sv["freq"] for sv in svs)

    # ------------------------------------------------------------------
    # Define three panels based on observed position clusters
    #   Panel 1 (overview) : full range of all breakpoints
    #   Panel 2 (mid-range): ~8.5 Mb – 15 Mb cluster
    #   Panel 3 (MYCN zoom): tight window around MYCN locus
    # ------------------------------------------------------------------
    all_pos = sorted(
        p for sv in svs if sv["chr1"] == sv["chr2"]
        for p in (sv["pos1"], sv["pos2"])
    )
    pad = 500_000
    full_min = all_pos[0]  - pad
    full_max = all_pos[-1] + pad

    fig, axes = plt.subplots(
        3, 1, figsize=(16, 13),
        gridspec_kw={"height_ratios": [1, 1, 1]},
    )
    fig.subplots_adjust(
        left=0.05, right=0.86,
        top=0.94, bottom=0.07,
        hspace=0.65,
    )
    fig.suptitle(
        f"Top {top_n} recurrent SVs in MYCN-amplified samples",
        fontsize=13, fontweight="bold", y=0.98,
    )

    # Panel 1: full overview
    draw_panel(
        axes[0], svs,
        x_min=full_min, x_max=full_max,
        gene_start=args.gene_start, gene_end=args.gene_end,
        gene_label="MYCN",
        title="Overview — full chr2 range",
        max_freq=max_freq,
        label_svs=True,
    )

    # Panel 2: mid-range cluster (~8.5–15 Mb)
    draw_panel(
        axes[1], svs,
        x_min=8_000_000, x_max=15_200_000,
        gene_start=args.gene_start, gene_end=args.gene_end,
        gene_label="MYCN",
        title="Mid-range cluster (chr2 ~8.5–15.2 Mb)",
        max_freq=max_freq,
        label_svs=True,
    )

    # Panel 3: tight zoom on MYCN breakpoint cluster (~16.22–16.23 Mb)
    draw_panel(
        axes[2], svs,
        x_min=16_220_000, x_max=16_235_000,
        gene_start=args.gene_start, gene_end=args.gene_end,
        gene_label="MYCN",
        title="Zoomed — MYCN breakpoint cluster (chr2 ~16.220–16.235 Mb)",
        max_freq=max_freq,
        label_svs=True,
    )

    # ------------------------------------------------------------------
    # Shared colourbar — discrete palette
    # ------------------------------------------------------------------
    palette = {1: "#4a90d9", 2: "#2166ac", 3: "#f4a82e", 4: "#d62728"}
    freq_levels = sorted(set(sv["freq"] for sv in svs))
    colors_list = [palette.get(f, "#888888") for f in freq_levels]
    cmap_disc = mcolors.ListedColormap(colors_list)
    bounds = [f - 0.5 for f in freq_levels] + [freq_levels[-1] + 0.5]
    norm_disc = mcolors.BoundaryNorm(bounds, cmap_disc.N)
    sm = plt.cm.ScalarMappable(cmap=cmap_disc, norm=norm_disc)
    sm.set_array([])
    cbar = fig.colorbar(
        sm, ax=axes,
        orientation="vertical",
        fraction=0.018, pad=0.015, shrink=0.85, aspect=30,
        ticks=freq_levels,
    )
    cbar.set_label("Frequency\n(# samples)", fontsize=9)
    cbar.set_ticklabels([str(f) for f in freq_levels])

    # ------------------------------------------------------------------
    # Legend for frequency tiers
    # ------------------------------------------------------------------
    legend_handles = []
    for f in sorted(set(sv["freq"] for sv in svs), reverse=True):
        count = sum(1 for s in svs if s["freq"] == f)
        legend_handles.append(
            mpatches.Patch(
                facecolor=frequency_to_color(f, max_freq),
                edgecolor="none",
                label=f"freq = {f}  ({count} SV{'s' if count > 1 else ''})",
            )
        )
    axes[0].legend(
        handles=legend_handles,
        title="Recurrence", title_fontsize=8.5,
        fontsize=8, loc="upper left",
        framealpha=0.85,
    )

    plt.savefig(output_plot, dpi=150)
    plt.close()
    print(f"Wrote: {output_plot}")


main()