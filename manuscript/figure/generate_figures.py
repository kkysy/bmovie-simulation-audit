"""Generate the eight frozen-specification Bmovie manuscript figures.

Six main-text figures (Figure_1..Figure_6) and two supplementary figures
(Figure_S1, Figure_S2).  This script reads only completed summary/output
tables.  It performs no simulation, statistical re-estimation, or manuscript
editing.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "manuscript" / "figure"
DATA = ROOT / "processed" / "subject" / "group" / "ieeg_methods_paper"
CENSUS = (ROOT / "processed" / "subject" / "group" / "ieeg_acc00_sim" / "BangYoureDead" /
          "checkpoints" / "edgeguard_fix_20260829" / "validator" /
          "task-bangyouredead_desc-acc00sim-feasible-domain-census.tsv")

BLUE = "#2C6DB2"       # Mpower
ORANGE = "#D97824"     # Mphase
INK = "#24313D"
GRAY = "#6E7781"
LIGHT_GRAY = "#E8EDF2"
GREEN = "#3E8E6F"
RED = "#B94A48"
PURPLE = "#7B5EA7"     # Circular-modulo row (distinct from RED, BLUE, ORANGE)


def setup():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 10,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    OUT.mkdir(parents=True, exist_ok=True)


def add_panel_labels(fig, panel_rows, offsets=None):
    """Place panel letters at the upper-left of each axes, above its title."""
    offsets = offsets or {}
    fig.canvas.draw()
    labels = []
    for row in panel_rows:
        for letter, ax in row:
            pos = ax.get_position()
            dx, dy = offsets.get(letter, (-0.060, 0.075))
            labels.append(fig.text(pos.x0 + dx, pos.y1 + dy, letter,
                                   transform=fig.transFigure, fontsize=11,
                                   fontweight="bold", va="bottom", ha="left", color=INK,
                                   gid=f"panel-letter-{letter}"))
    return labels


def check_layout(fig, name, panel_rows, labels):
    """Render-time checks for panel-letter alignment and unintended text overlap.

    Tick labels, axis labels, and data-attached annotations are an explicit
    whitelist: their proximity is controlled by Matplotlib's axis layout or by
    leader lines.  Panel letters, titles, legends, free annotations, and figure
    text are checked pairwise.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    report = []
    for row in panel_rows:
        row_labels = [label for label in labels if label.get_text() in {letter for letter, _ in row}]
        centers = [label.get_window_extent(renderer).y0 + label.get_window_extent(renderer).height / 2
                   for label in row_labels]
        if len(centers) > 1:
            delta = max(centers) - min(centers)
            assert delta < 2.0, f"{name}: panel-letter row misalignment {delta:.2f}px"
            report.append(f"letter-row Δ={delta:.2f}px")

    checked = list(labels)
    for ax in fig.axes:
        checked.extend([ax.title, ax._left_title, ax._right_title])
        legend = ax.get_legend()
        if legend is not None:
            checked.extend(legend.get_texts())
        checked.extend(text for text in ax.texts if text.get_gid() != "data-label")
    for legend in fig.legends:
        checked.extend(legend.get_texts())
    checked.extend(text for text in fig.texts if text.get_gid() is None)
    checked = [text for text in checked if isinstance(text, Text) and text.get_visible() and text.get_text()]
    unexpected = []
    for i, left in enumerate(checked):
        left_box = left.get_window_extent(renderer)
        for right in checked[i + 1:]:
            if left.axes is right.axes and left.get_text() == right.get_text():
                continue
            if left_box.overlaps(right.get_window_extent(renderer)):
                unexpected.append((left.get_text(), right.get_text()))
    assert not unexpected, f"{name}: unexpected text overlap {unexpected[:3]}"
    # Free plot annotations must not cover plotted lines or point collections.
    # Axis/tick labels and explicitly whitelisted contextual labels are excluded.
    text_data_overlaps = []
    for ax in fig.axes:
        free_text = [text for text in ax.texts
                     if text.get_visible() and text.get_text() and
                     text.get_gid() not in {"data-label", "text-data-whitelist"}]
        for text in free_text:
            text_box = text.get_window_extent(renderer)
            for artist in [*ax.lines, *ax.collections]:
                if not artist.get_visible():
                    continue
                artist_box = artist.get_window_extent(renderer)
                if (np.isfinite(artist_box.bounds).all() and text_box.overlaps(artist_box)):
                    text_data_overlaps.append((text.get_text(), type(artist).__name__))
    assert not text_data_overlaps, (
        f"{name}: unexpected text-data overlap {text_data_overlaps[:3]}")
    box_pairs = getattr(fig, "_box_text_pairs", [])
    for patch, text in box_pairs:
        patch_box = patch.get_window_extent(renderer)
        text_box = text.get_window_extent(renderer)
        assert (patch_box.contains(text_box.x0, text_box.y0) and
                patch_box.contains(text_box.x1, text_box.y1)), (
                    f"{name}: box {patch_box.bounds} does not contain text "
                    f"{text.get_text()!r} at {text_box.bounds}")
    for arrow, text in getattr(fig, "_arrow_text_pairs", []):
        arrow_box = arrow.arrow_patch.get_window_extent(renderer)
        assert not arrow_box.overlaps(text.get_window_extent(renderer)), (
            f"{name}: arrow intersects text {text.get_text()!r}")
    print(f"LAYOUT CHECK {name}: PASS; {', '.join(report) or 'one panel'}; text-text=0; text-data=0 (data-label whitelist); box-text={len(box_pairs)} contained")


def finish(fig, name, panel_rows=(), panel_label_offsets=None):
    labels = add_panel_labels(fig, panel_rows, panel_label_offsets) if panel_rows else []
    check_layout(fig, name, panel_rows, labels)
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def errorbar_from_bounds(ax, x, y, low, high, **kwargs):
    return ax.errorbar(x, y, yerr=np.vstack((y - low, high - y)), **kwargs)


def _legend_scale(ax):
    """Scale legend typography/handles to the physical axes size."""
    ax.figure.canvas.draw()
    bbox = ax.get_window_extent()
    return float(np.clip(min(bbox.width, bbox.height) / 300.0, 0.72, 1.10))


def _legend_kwargs(ax):
    scale = _legend_scale(ax)
    return dict(loc="upper left", ncol=1, fontsize=8.0 * scale,
                frameon=True, fancybox=False, facecolor="white",
                edgecolor="black", framealpha=0.96,
                borderpad=0.70 * scale, labelspacing=0.45 * scale,
                handletextpad=0.65 * scale, handlelength=1.35 * scale,
                alignment="center")


def _style_legend(legend):
    legend.get_frame().set_linewidth(0.5)
    return legend


def route_legend(ax, color, **kwargs):
    """R1/R2 legend: run is encoded by fill state, not by colour."""
    extra = list(kwargs.pop("extra", []))
    scale = _legend_scale(ax)
    handles = [
        Line2D([], [], marker="o", linestyle="None", markersize=6.0 * scale,
               markerfacecolor=color, markeredgecolor=color, color=color, label="R1"),
        Line2D([], [], marker="s", linestyle="None", markersize=6.0 * scale,
               markerfacecolor="white", markeredgecolor=color, color=color, label="R2"),
    ]
    legend_kwargs = _legend_kwargs(ax)
    legend_kwargs.update(kwargs)
    return _style_legend(ax.legend(handles=handles + extra, **legend_kwargs))


def bar_legend(ax, color=INK, **kwargs):
    """R1/R2 legend using rectangular swatches for bar charts."""
    handles = [
        patches.Patch(facecolor=color, edgecolor=color, label="R1"),
        patches.Patch(facecolor="white", edgecolor=color, label="R2"),
    ]
    legend_kwargs = _legend_kwargs(ax)
    legend_kwargs.update(kwargs)
    return _style_legend(ax.legend(handles=handles, **legend_kwargs))


def fmt_signed(value, digits=2):
    return f"{value:+.{digits}f}"


def figure1():
    fig, ax = plt.subplots(figsize=(7.48, 6.4))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    def box(x, y, w, h, text, color, fs=6.0, ls="-"):
        rect = patches.FancyBboxPatch((x, y), w, h,
                                      boxstyle="round,pad=0.012,rounding_size=0.018",
                                      facecolor=color, edgecolor=INK, linewidth=0.9,
                                      linestyle=ls)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                color=INK)

    def arrow(x1, y1, x2, y2):
        return ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(
            arrowstyle="->", color=INK, lw=1.1, shrinkA=3, shrinkB=3))

    ax.text(0.50, 0.985, "Design-faithful audit framework for event-related inference",
            fontsize=10.5, fontweight="bold", va="top", ha="center", color=INK)
    ax.text(0.98, 0.952, "Schematic, not a data result", fontsize=7.2,
            va="top", ha="right", color=GRAY, style="italic")
    # Two bands: the upper band is the reusable audit sequence (general
    # framework), the lower band is the Bmovie implementation example that
    # realizes it. Frozen stage tags (P0-1..P0-7) are carried in parentheses
    # on the framework boxes; numeric constants live in Methods. Boxes are
    # kept inside the axes (the round pad extends past the given rect).
    ax.text(0.02, 0.906, "General framework: reusable audit sequence",
            fontsize=7.6, color=GRAY, va="bottom")
    ax.plot([0.700, 0.740], [0.9125, 0.9125], color=INK, lw=0.9, ls=(0, (4, 2)))
    ax.text(0.98, 0.906, "dashed = conditional steps", fontsize=7.0,
            color=GRAY, va="bottom", ha="right")

    # Upper band, row 1: design -> worlds -> null calibration -> power.
    box(0.02, 0.745, 0.205, 0.125,
        "1. Observed design\nevent table · support\nhierarchy · frozen path", "#F4F7F9")
    box(0.265, 0.745, 0.205, 0.125, "2. Synthetic worlds\ndesign-faithful,\nnull and injected", "#F4F7F9")
    box(0.510, 0.745, 0.205, 0.125, "3. Null calibration\nempirical FPR + CI\n(P0-1)", "#F4F7F9")
    box(0.755, 0.745, 0.205, 0.125, "4. Estimand-matched power\ninjection · power · bias\n· coverage (P0-2)", "#F4F7F9")
    arrow(0.225, 0.8075, 0.265, 0.8075)
    arrow(0.470, 0.8075, 0.510, 0.8075)
    arrow(0.715, 0.8075, 0.755, 0.8075)

    # Upper band, row 2 (right to left): mechanism -> repair/revalidate ->
    # blinded -> surrogate -> perturbation. Steps 5, 6, and 7 are conditional
    # (Table 2), so their boxes carry dashed borders and trigger conditions.
    arrow(0.8745, 0.745, 0.8745, 0.635)
    box(0.792, 0.51, 0.165, 0.125, "5. Mechanism diagnosis\n(P0-3) if calibration\nor sensitivity fails", "#F4F7F9",
        fs=5.6, ls=(0, (4, 2)))
    box(0.599, 0.51, 0.165, 0.125, "6. Repair & revalidate\n(P0-7) on fresh null and\ninjection-matched worlds", "#F4F7F9",
        fs=5.6, ls=(0, (4, 2)))
    box(0.406, 0.51, 0.165, 0.125, "7. Blinded implementation\nstress test (P0-4) when\ndecision-relevant", "#F4F7F9",
        fs=5.6, ls=(0, (4, 2)))
    box(0.213, 0.51, 0.165, 0.125, "8. Surrogate audit (P0-5)\ninvariant · structure ·\nfeasible domain", "#F4F7F9",
        fs=5.6)
    box(0.02, 0.51, 0.165, 0.125, "9. Design perturbation\nrobustness (P0-6)", "#F4F7F9",
        fs=5.6)
    arrow(0.792, 0.5725, 0.764, 0.5725)
    arrow(0.599, 0.5725, 0.571, 0.5725)
    arrow(0.406, 0.5725, 0.378, 0.5725)
    arrow(0.213, 0.5725, 0.185, 0.5725)

    # Upper band, row 3: interpretation gate, entered from the perturbation box.
    arrow(0.1025, 0.51, 0.1025, 0.405)
    box(0.02, 0.28, 0.45, 0.125,
        "10. Interpretation\nafter calibration, sensitivity, support, perturbation\nrobustness, and any required repair–revalidation", "#EEF6F1")

    ax.text(0.02, 0.215, "Bmovie analysis paths audited by the framework", fontsize=7.6,
            color=GRAY, va="bottom")

    # Lower band: the observed design (schedule, support, hierarchy) feeds
    # design-faithful synthetic worlds, which are scored by two stacked routes
    # converging into the shared inference at the far right.
    box(0.02, 0.075, 0.21, 0.12,
        "Observed Bmovie design\n1,901 events · bad support\nparticipant → session → pair", "#EAF2FB", 5.8)
    box(0.27, 0.075, 0.21, 0.12,
        "Design-faithful synthetic worlds\npreserve schedule, clean\nsupport, and hierarchy", "#EAF2FB", 5.8)
    box(0.52, 0.125, 0.20, 0.055, "Power change (Mpower)", "#EAF2FB", 6.2)
    box(0.52, 0.05, 0.20, 0.055, "Phase consistency (Mphase)", "#FDEFE2", 6.2)
    box(0.76, 0.075, 0.22, 0.12, "Shared inference\nexact sign-flip max-T,\nFWER 0.05",
        "#F2F5F7", 6.0)
    arrow(0.23, 0.135, 0.27, 0.135)
    arrow(0.48, 0.15, 0.52, 0.1525)
    arrow(0.48, 0.10, 0.52, 0.0775)
    arrow(0.72, 0.1525, 0.76, 0.145)
    arrow(0.72, 0.0775, 0.76, 0.095)
    finish(fig, "Figure_1")


def figure2():
    cal = pd.read_csv(DATA / "p0_1_calibration_summary" / "p0_1_null_calibration.tsv", sep="\t")
    cal = cal[cal["estimator_family"].isin(["Mpower", "Mphase"])].copy()
    fig = plt.figure(figsize=(7.48, 3.15))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.1], wspace=0.55)
    axa, axb = [fig.add_subplot(gs[0, i]) for i in range(2)]

    order = [("Mpower", "R1"), ("Mpower", "R2"), ("Mphase", "R1"), ("Mphase", "R2")]
    labels = [f"{fam} {run}" for fam, run in order]
    for y, (fam, run) in enumerate(order):
        row = cal[(cal.estimator_family == fam) & (cal.run == run)].iloc[0]
        color = BLUE if fam == "Mpower" else ORANGE
        marker = "o" if run == "R1" else "s"
        axa.scatter(row.null_mean_over_sd, y, s=44, color=color, marker=marker, zorder=3)
        # Unit is carried by the x-axis label, so no suffix on either route.
        axa.text(row.null_mean_over_sd + 0.07, y,
                 fmt_signed(row.null_mean_over_sd, 2),
                 va="center", fontsize=7.2)
    axa.axvline(0, color=GRAY, lw=0.9, zorder=0)
    axa.set_yticks(range(4), labels)
    axa.invert_yaxis()
    axa.set_xlim(-1.95, 0.35)
    axa.set_xlabel("Null-center displacement (SD)")
    axa.set_title("Null center", loc="center", x=0.50, pad=10)

    xbase = {"Mpower": 0.0, "Mphase": 1.0}
    bar_offsets = {"R1": -0.18, "R2": 0.18}
    bar_width = 0.30
    for fam in ["Mpower", "Mphase"]:
        color = BLUE if fam == "Mpower" else ORANGE
        for run in ["R1", "R2"]:
            row = cal[(cal.estimator_family == fam) & (cal.run == run)].iloc[0]
            x = xbase[fam] + bar_offsets[run]
            axb.bar(x, row.fpr, width=bar_width, color=color if run == "R1" else "white",
                    edgecolor=color, linewidth=1.1, zorder=2)
            axb.errorbar(x, row.fpr, yerr=[[row.fpr - row.ci95_low],
                                           [row.ci95_high - row.fpr]],
                         fmt="none", ecolor=INK, capsize=2.5, lw=1.1, zorder=4)
    axb.axhline(0.05, color=GRAY, lw=0.9, ls="--")
    axb.text(0.5, 0.058, "nominal 0.05", fontsize=7, color=GRAY, ha="center")
    axb.set_xticks([0, 1], ["Mpower", "Mphase"])
    axb.set_xlim(-0.45, 1.45)
    axb.set_ylim(0, 0.75)
    axb.set_ylabel("False-positive rate (95% CP interval)")
    axb.set_title("Null FPR", loc="center", x=0.50, pad=10)
    bar_legend(axb)
    finish(fig, "Figure_2", [[("A", axa), ("B", axb)]])


def figure3():
    d = pd.read_csv(DATA / "p0_2_sim_powercurve" / "tables" /
                    "p0_2_powercurve_summary.tsv", sep="\t")
    grid = [f"G{i:02d}" for i in range(1, 8)]
    fig = plt.figure(figsize=(7.48, 6.55))
    gs = GridSpec(2, 2, figure=fig, wspace=0.42, hspace=0.92,
                  height_ratios=[1.0, 0.92])
    axa, axb = [fig.add_subplot(gs[0, i]) for i in range(2)]
    axc = fig.add_subplot(gs[1, :])
    offsets = {"R1": -0.10, "R2": 0.10}
    markers = {"R1": "o", "R2": "s"}
    for run, color in [("R1", BLUE), ("R2", BLUE)]:
        s = d[d.run == run].set_index("grid_id").loc[grid].reset_index()
        x = np.arange(7) + offsets[run]
        errorbar_from_bounds(axa, x, s.power_hat.to_numpy(), s.power_cp95_low.to_numpy(),
                             s.power_cp95_high.to_numpy(), fmt=markers[run] + "-", color=color,
                             mfc="white" if run == "R2" else color, capsize=2.5, lw=1.15,
                             markersize=4.2, label=run)
        errorbar_from_bounds(axb, x, s.coverage_rate.to_numpy(), s.coverage_cp95_low.to_numpy(),
                             s.coverage_cp95_high.to_numpy(), fmt=markers[run] + "-", color=color,
                             mfc="white" if run == "R2" else color, capsize=2.5, lw=1.15,
                             markersize=4.2, label=run)
        axc.errorbar(x, s.observed_slope_mean.to_numpy(),
                     yerr=s.observed_slope_sample_sd.to_numpy(),
                     fmt=markers[run] + "-", color=color,
                     mfc="white" if run == "R2" else color, capsize=2.5,
                     lw=1.15, markersize=4.2, label=run)
    for ax in [axa, axb]:
        ax.set_xticks(range(7), grid)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="center")
        for x, label in zip([4, 5, 6], ["0.5×", "1×", "2×"]):
            ax.text(x, -0.18, label, transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=7.0)
        # Sit the axis label just below the 0.5x/1x/2x second-row annotations;
        # a larger labelpad drops it into panel C's letter/title band.
        ax.set_xlabel("Additive grid cell", labelpad=12)
    for ax in [axa, axb, axc]:
        # Separate the near-null cells (G01-G04) from the reference-scale
        # injections (G05-G07) without adding text.
        ax.axvline(3.5, color=GRAY, lw=0.7, ls=":", zorder=0)
    axa.axhline(0.05, color=GRAY, ls="--", lw=0.9)
    axa.text(6.12, 0.06, "0.05", fontsize=7, color=GRAY, gid="data-label")
    # G07's R1 interval reaches ~0.26; the ceiling adds headroom above it for
    # the in-panel region labels flanking the dotted separator.
    axa.set_ylim(0, 0.34)
    axa.set_ylabel("Empirical rejection rate (95% CP interval)")
    # pad clears the tall y-axis labels, whose tops reach above the axes into
    # the title band (panel B's title otherwise collides with its ylabel).
    axa.set_title("Empirical rejection rate remained below 0.12", loc="center", x=0.50, pad=24, fontsize=9)
    route_legend(axa, BLUE)

    # Leave an upper band for the common upper-left R1/R2 legend and the
    # in-panel region labels flanking the dotted separator.
    axb.set_ylim(0, 1.52)
    axb.set_ylabel("Coverage proportion (95% CP interval)")
    axb.set_title("Coverage collapsed to 0.00 at the largest injections", loc="center", x=0.50, pad=24, fontsize=9)
    route_legend(axb, BLUE)

    truth = d[d.run == "R1"].set_index("grid_id").loc[grid, "true_bridge_slope"].to_numpy()
    axc.plot(np.arange(7), truth, "D--", color=INK, lw=1.15, markersize=4.0,
             markerfacecolor="white", label="Injected truth")
    axc.axhline(0, color=GRAY, lw=0.9)
    axc.set_xticks(range(7), grid)
    for x, label in zip([4, 5, 6], ["0.5×", "1×", "2×"]):
        axc.text(x, -0.16, label, transform=axc.get_xaxis_transform(),
                 ha="center", va="top", fontsize=7.0)
    axc.set_xlabel("Additive grid cell", labelpad=31)
    axc.set_ylabel("Slope (log10 power per z)")
    axc.ticklabel_format(axis="y", style="sci", scilimits=(-3, -3))

    # Region labels live inside each panel at the top, flanking the dotted
    # separator (two-line form so they clear the line and the data traces);
    # panel y ceilings are raised to make room for them.
    def region_labels(ax, y):
        ax.text(3.30, y, "near-null\ninjections", ha="right", va="top",
                fontsize=7.0, color=GRAY, style="italic", linespacing=1.25)
        ax.text(3.70, y, "reference-scale\ninjections", ha="left", va="top",
                fontsize=7.0, color=GRAY, style="italic", linespacing=1.25)

    region_labels(axa, 0.335)
    region_labels(axb, 1.505)
    c_lo, c_hi = axc.get_ylim()
    axc.set_ylim(c_lo, c_hi * 1.50)
    region_labels(axc, c_hi * 1.44)
    axc.set_title("Observed slopes attenuated as injected truth increased",
                  loc="center", x=0.50, pad=10, fontsize=9)
    handles = [
        Line2D([], [], marker="o", linestyle="-", markersize=4.5,
               markerfacecolor=BLUE, markeredgecolor=BLUE, color=BLUE, label="R1 observed"),
        Line2D([], [], marker="s", linestyle="-", markersize=4.5,
               markerfacecolor="white", markeredgecolor=BLUE, color=BLUE, label="R2 observed"),
        Line2D([], [], marker="D", linestyle="--", markersize=4.2,
               markerfacecolor="white", markeredgecolor=INK, color=INK, label="Injected truth"),
    ]
    _style_legend(axc.legend(handles=handles, **_legend_kwargs(axc)))
    finish(fig, "Figure_4", [[("A", axa), ("B", axb)], [("C", axc)]])


def figure4():
    contrasts = pd.read_csv(DATA / "p0_3_ablation" / "tables" /
                            "p0_3_mpower_factor_contrasts.tsv", sep="\t")
    diagnostic = pd.read_csv(DATA / "p0_3_ablation" / "tables" /
                             "p0_3_mpower_pair_diagnostics.tsv.gz", sep="\t")
    replay = pd.read_csv(DATA / "p0_3_ablation" / "tables" /
                         "p0_3_mphase_replay_shape_summary.tsv", sep="\t")
    hist = pd.read_csv(DATA / "p0_3_ablation" / "tables" /
                       "p0_3_mphase_historical_guard_summary.tsv", sep="\t")
    fig = plt.figure(figsize=(7.48, 6.6))
    gs = GridSpec(2, 2, figure=fig, wspace=0.38, hspace=0.62,
                  height_ratios=[1.05, 1.0])
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    sub = gs[1, 0].subgridspec(2, 1, hspace=0.82, height_ratios=[1.2, 1.0])
    axc1 = fig.add_subplot(sub[0, 0])
    axc2 = fig.add_subplot(sub[1, 0])
    axd = fig.add_subplot(gs[1, 1])

    corder = list(dict.fromkeys(contrasts.contrast.tolist()))
    # Machine-readable contrast names from the summary table are not readable
    # axis labels; print the direction and factor instead and keep the verbatim
    # names in the caption table.
    READABLE = {
        "pre_minus_post_canonical": "Pre − post (canonical)",
        "pre_minus_post_half": "Pre − post (half)",
        "pre_minus_post_quarter": "Pre − post (quarter)",
        "symmetric_minus_post_canonical": "Symmetric − post (canonical)",
        "symmetric_minus_post_half": "Symmetric − post (half)",
        "symmetric_minus_post_quarter": "Symmetric − post (quarter)",
        "post_half_minus_canonical": "Half − canonical (post)",
        "post_quarter_minus_canonical": "Quarter − canonical (post)",
        "P10_minus_P01": "Theta carrier − broadband (P10 − P01)",
    }
    for run, marker, offset, face in [("R1", "o", -0.22, BLUE), ("R2", "s", 0.22, "white")]:
        s = contrasts[contrasts.run == run].set_index("contrast").loc[corder].reset_index()
        y = np.arange(len(corder)) + offset
        axa.errorbar(s.paired_difference_mean, y, xerr=s.paired_difference_sample_sd,
                     fmt=marker, color=BLUE, mfc=face, mec=BLUE, capsize=2,
                     markersize=4.4, lw=1.0, label=run)
    axa.set_yticks(np.arange(len(corder)), [READABLE[c] for c in corder], fontsize=6.2)
    axa.invert_yaxis()
    # Reserve a blank band above the first contrast so the common legend
    # does not cover the near-zero estimates at the top of the panel.
    axa.set_ylim(len(corder) - 0.5, -2.8)
    # P10_minus_P01 is about two orders of magnitude larger than every
    # placement/retention contrast, so the panel reads as uniformly null unless
    # that one contrast is labelled numerically.
    carrier = "P10_minus_P01"
    if carrier in corder:
        sub = contrasts[contrasts.contrast == carrier]
        tip = sub.paired_difference_mean.max() + sub.paired_difference_sample_sd.max()
        axa.text(tip + 0.006, corder.index(carrier), f"+{sub.paired_difference_mean.max():.3f}",
                 va="center", ha="left", fontsize=6.6, color=BLUE)
        axa.set_xlim(-0.0035, tip + 0.028)
    axa.set_xlabel("Paired mean difference (sample-SD error bars)")
    axa.set_title("Mpower factor contrasts", loc="center", x=0.50, pad=10)
    route_legend(axa, BLUE)

    # P04-P06 are pre-only cells. The post-component overlap metric is
    # undefined for them and is stored as NaN; exclude those cells rather than
    # presenting empty categories as if they were missing plotted values.
    cells = [cell for cell in sorted(diagnostic.cell_id.unique())
             if cell not in {"P04", "P05", "P06"}]
    vals = [diagnostic.loc[diagnostic.cell_id == c,
                           "predecessor_in_pre_overlap_fraction_post_component"].to_numpy()
            for c in cells]
    violin = axb.violinplot(vals, positions=np.arange(len(cells)) + 1,
                            showmeans=False, showmedians=False, showextrema=False, widths=0.82)
    for body in violin["bodies"]:
        body.set_facecolor(BLUE)
        body.set_edgecolor(BLUE)
        body.set_alpha(0.55)
    jitter_rng = np.random.default_rng(0)
    for ix, values in enumerate(vals, start=1):
        axb.scatter(ix + jitter_rng.uniform(-0.23, 0.23, size=len(values)), np.clip(values, 0, 1),
                    s=0.18, color=BLUE, alpha=0.035, edgecolors="none", rasterized=True, zorder=2)
        # Full-retention cells are constant at 0.9995, not 1.0, so the
        # saturation strip must trigger above 0.99 rather than at isclose(1).
        saturated = values[values > 0.99]
        if len(saturated):
            # Deterministic in-cell spread keeps the observed saturation visible.
            axb.scatter(ix + np.linspace(-0.24, 0.24, len(saturated)), saturated,
                        s=1.4, color=BLUE, alpha=0.22, edgecolors="none",
                        rasterized=True, zorder=3, clip_on=False)
    # Human-readable two-row tick labels: cell ID over retention. The
    # post-versus-symmetric placement split rides on the axis label because
    # full "placement·retention" words do not fit the narrow category width.
    CELL2ROWS = {"P01": "P01\nfull", "P02": "P02\nhalf", "P03": "P03\nquarter",
                 "P07": "P07\nfull", "P08": "P08\nhalf", "P09": "P09\nquarter",
                 "P10": "P10\ntheta"}
    axb.set_xticks(np.arange(len(cells)) + 1, [CELL2ROWS[c] for c in cells], fontsize=6.4)
    axb.set_ylim(0, 1)
    # Full-retention cells (P01, P07, P10) sit at 0.9995, indistinguishable
    # from a clipped axis bound; print the medians so the saturation reads as
    # a value rather than a drawing artifact.
    for c in cells:
        med = float(np.nanmedian(diagnostic.loc[diagnostic.cell_id == c,
                                                "predecessor_in_pre_overlap_fraction_post_component"]))
        if med > 0.99:
            axb.text(cells.index(c) + 1, 1.01, f"{med:.4f}", ha="center", va="bottom",
                     fontsize=6.0, color=BLUE, clip_on=False)
    axb.set_xlabel("P0-3 Mpower cell\n(P01–03 post, P07–09 symmetric placement)", fontsize=7.2)
    axb.set_ylabel("Predecessor post-component\noverlap fraction in pre window")
    axb.set_title("Predecessor-overlap diagnostic", loc="center", x=0.50, pad=10)

    cat = [("M01", "R1"), ("M01", "R2"), ("M02", "R1"), ("M02", "R2")]
    for xi, (cell, run) in enumerate(cat):
        s = replay[(replay.cell_id == cell) & (replay.run == run)]
        axc1.scatter(np.full(len(s), xi), s.observed_standardized_score,
                     s=10, color=ORANGE, alpha=0.45, edgecolors="none")
    axc1.axhline(0, color=GRAY, lw=0.8)
    # Mark the nominal two-sided sign-flip rejection band so the displaced
    # replay mass reads against the thresholds rather than by inspection.
    for thr in (-1.96, 1.96):
        axc1.axhline(thr, color=GRAY, ls="--", lw=0.9)
    axc1.set_ylim(-3.4, 2.5)
    # Lead with the scientific variable (guard durations); the M01/M02 cell
    # IDs are carried by the caption instead of the axis.
    GUARD2ROWS = {"M01": "0.75/0.756 s", "M02": "1.5/1.5 s"}
    axc1.set_xticks(range(4), [f"{GUARD2ROWS[c]}\n{r}" for c, r in cat], fontsize=6.4)
    axc1.set_ylabel("Observed standardized score")
    axc1.set_title("Mphase guard-replay standardized scores", loc="center", x=0.50, fontsize=7.2, pad=12)

    guards = [("Mphase_guard_0p75", "0.75 / 0.756 s"), ("Mphase_guard_1p5", "1.5 / 1.5 s")]
    guard_x = [0.0, 1.0]
    for xi, (cell, name) in enumerate(guards):
        for run, face in [("R1", ORANGE), ("R2", "white")]:
            row = hist[(hist.reference_cell == cell) & (hist.run == run)].iloc[0]
            x = guard_x[xi] + (-0.18 if run == "R1" else 0.18)
            axc2.bar(x, row.fpr, width=0.30, color=face, edgecolor=ORANGE,
                     linewidth=1.1, zorder=2)
    axc2.set_xticks(guard_x, [x[1] for x in guards])
    axc2.set_xlim(-0.40, 1.40)
    axc2.set_ylim(0, 1.00)
    axc2.set_ylabel("FPR")
    axc2.set_xlabel("Historical guard configuration")
    axc2.set_title("Separate 400-world historical references", loc="center", x=0.50, fontsize=7.2, pad=10)
    bar_legend(axc2, color=ORANGE, loc="upper right")

    predicted = np.array([0.37, 0.30])
    observed = np.array([0.3875, 0.3275])
    for run, x, y, marker in zip(["R1", "R2"], predicted, observed, ["o", "s"]):
        axd.scatter(x, y, s=47, marker=marker, facecolor="white" if run == "R2" else ORANGE,
                    edgecolor=ORANGE, zorder=3)
    axd.plot([0.26, 0.41], [0.26, 0.41], color=GRAY, ls="--", lw=0.9)
    axd.set_xlim(0.26, 0.41)
    axd.set_ylim(0.26, 0.41)
    axd.set_xlabel("Normal-shift predicted FPR")
    axd.set_ylabel("Observed FPR")
    axd.set_title("Mphase center-displacement check", loc="center", x=0.50, pad=10)
    route_legend(axd, ORANGE)
    finish(fig, "Figure_5", [[("A", axa), ("B", axb)], [("C", axc1), ("D", axd)]])


def figure5():
    summary = pd.read_csv(DATA / "p0_4_diagnostic_panel" / "tables" /
                          "p0_4_member_fpr_summary.tsv", sep="\t")
    paired = pd.read_csv(DATA / "p0_4_diagnostic_panel" / "tables" /
                         "p0_4_paired_contrast_MD2.tsv", sep="\t")
    worlds = pd.read_csv(DATA / "p0_4_diagnostic_panel" / "tables" /
                         "p0_4_member_worlds.tsv", sep="\t")
    decoded = {
        "PanelMember-01": "P-REF",
        "PanelMember-02": "M-REF",
        "PanelMember-03": "M-D2",
        "PanelMember-04": "M-D1\n(exploratory)",
    }
    ids = list(decoded)
    fig = plt.figure(figsize=(7.48, 6.5))
    fig.subplots_adjust(bottom=0.16)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 1.2], hspace=0.58, wspace=0.52)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    sub = gs[1, :].subgridspec(1, 2, wspace=0.09)
    axc1, axc2 = [fig.add_subplot(sub[0, i]) for i in range(2)]
    # Shrink the raster axes vertically (bottom-anchored) to open a legend
    # band between the R1/R2 titles and the axes tops; titles and the panel
    # letter keep their original heights.
    for ax in (axc1, axc2):
        b = ax.get_position()
        ax.set_position([b.x0, b.y0, b.width, b.height * 0.87])

    member_x = [0.0, 1.25, 2.50, 3.75]
    for xi, opaque in enumerate(ids):
        color = BLUE if opaque == "PanelMember-01" else ORANGE
        for run, offset, marker, face in [("R1", -0.18, "o", color), ("R2", 0.18, "s", "white")]:
            row = summary[(summary.opaque_id == opaque) & (summary.run == run)].iloc[0]
            x = member_x[xi] + offset
            errorbar_from_bounds(axa, x, row.fpr, row.cp95_low, row.cp95_high,
                                 fmt=marker, color=color, mfc=face, mec=color,
                                 capsize=2.5, lw=1.1, markersize=5.0, zorder=3)
    axa.axhline(0.05, color=GRAY, lw=0.9, ls="--")
    axa.text(4.15, 0.061, "nominal 0.05", ha="right", fontsize=6.8, color=GRAY,
             gid="data-label")
    axa.set_xticks(member_x, [decoded[x] for x in ids])
    axa.set_xlim(-0.58, 4.33)
    axa.set_ylim(0, 1.10)
    axa.set_ylabel("False-positive rate (95% CP interval)")
    # The long y-label's top end grazes the panel title's left end; shift the
    # label outward to clear it.
    axa.yaxis.set_label_coords(-0.185, 0.5)
    axa.set_title("Member FPRs showed distinct null patterns", loc="center", x=0.50, pad=10)
    route_legend(axa, INK)

    for xi, run in enumerate(["R1", "R2"]):
        row = paired[paired.run == run].iloc[0]
        axb.scatter(xi, row.risk_difference, s=52, color=ORANGE, zorder=3)
        axb.vlines(xi, 0, row.risk_difference, color=ORANGE, lw=1.25)
        axb.text(xi, row.risk_difference + 0.015,
                 f"K10 = {int(row.K10)}, K01 = {int(row.K01)}\n"
                 f"+{row.risk_difference:.4f}; p = {row.exact_mcnemar_p_two_sided:.5f}",
                 ha="center", va="bottom", fontsize=7)
    axb.axhline(0, color=GRAY, lw=0.9)
    axb.set_xlim(-0.55, 1.55)
    axb.set_ylim(-0.025, 0.30)
    axb.set_xticks([0, 1], ["R1", "R2"])
    axb.set_ylabel("Paired risk difference (M-D2 − M-REF)")
    axb.set_title("M-D2 exceeded M-REF in both runs", loc="center", x=0.50, fontsize=8.5, pad=10)

    for ax, run in [(axc1, "R1"), (axc2, "R2")]:
        subset = worlds[worlds.run == run].copy()
        subset["y"] = subset.opaque_id.map({opaque: i for i, opaque in enumerate(ids)})
        md1_row = subset.y == 3
        for rejected, color, label in [(0, LIGHT_GRAY, "not rejected"), (1, RED, "rejected")]:
            # M-D1 is exploratory; render its full row at reduced saturation so
            # the raster's visual weight stays with the decision members.
            for mask, a in [(~md1_row, 1.0), (md1_row, 0.42)]:
                s = subset[(subset.reject_alpha_0_05 == rejected) & mask]
                ax.scatter(s.seed_index, s.y, s=12, marker="s", color=color,
                           alpha=a, edgecolor="none",
                           label=label if (a == 1.0) else None, rasterized=True)
        ax.set_xlim(0, 121)
        ax.set_ylim(3.5, -0.5)
        ax.set_xticks([1, 30, 60, 90, 120])
        ax.set_xlabel("Common-random-number seed index")
        # Raise the titles clear of the shared legend band below them; the
        # panel-C letter bounds how far they may go.
        ax.set_title(run, pad=24)
        if ax is axc1:
            ax.set_yticks(range(4), [decoded[x] for x in ids])
            ax.set_ylabel("Panel member", labelpad=1)
            # Long member tick labels occupy the space immediately left of the
            # spine; keep the title close while leaving a clean separation.
            ax.yaxis.set_label_coords(-0.16, 0.5)
        else:
            ax.set_yticks(range(4), [])
    handles, labels = axc1.get_legend_handles_labels()
    fig.text(0.5, 0.080, "World-level rejection raster", ha="center", fontsize=9.2)
    # Anchor the legend between the R1/R2 subplot titles and the shrunk axes
    # tops, clearing the R2 y-axis.
    fig.legend(handles, labels, frameon=False, loc="center", ncol=2,
               bbox_to_anchor=(0.5, 0.442), fontsize=8.5,
               handlelength=1.0,
               handletextpad=0.45, columnspacing=0.75)
    finish(fig, "Figure_S1", [[("A", axa), ("B", axb)], [("C", axc1)]],
           panel_label_offsets={"C": (-0.085, 0.084)})


def _merge_intervals(intervals):
    """Merge overlapping or touching intervals on [0, 1]; input need not be sorted."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1e-12:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _complement_intervals(intervals):
    """Complement of already-merged intervals within [0, 1]."""
    gaps = []
    previous = 0.0
    for start, end in intervals:
        if start > previous + 1e-12:
            gaps.append((previous, start))
        previous = max(previous, end)
    if previous < 1.0 - 1e-12:
        gaps.append((previous, 1.0))
    return gaps


def figure6():
    census = pd.read_csv(CENSUS, sep="\t")
    frac = census.domain_fraction.to_numpy() * 100.0
    med = float(np.median(frac))

    fig = plt.figure(figsize=(7.48, 3.95))
    fig.subplots_adjust(top=0.80, bottom=0.17, left=0.10, right=0.98)
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.05, 1.0], wspace=0.40)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])

    axa.set_xlim(0, 1)
    axa.set_ylim(-0.1, 3.3)
    # Tick marks sit at band centers: each band is 0.38 tall, so the center is
    # its bottom edge plus 0.19.
    axa.set_yticks([0.44, 1.44, 2.44],
                   ["Prohibited", "Circular modulo", "Feasible intervals"])
    axa.set_xticks([0, 1.0], ["0", "movie schedule"])
    axa.set_xlabel("Circular shift δ")
    axa.set_title("Feasible-domain construction", loc="center", x=0.50,
                  fontsize=10.5, fontweight="bold", color=INK, pad=10)
    # Schematic interval construction only; no selected empirical pair. The
    # three bands are mutually consistent and non-overlapping: the lower band
    # shows the dilated prohibitions after circular mapping and merging, the
    # middle band highlights the two slivers that the wrap-around produced,
    # and the upper band is the exact complement of the lower band. Widths remain illustrative: the drawn feasible slivers are
    # far wider than the audited median in Panel B (0.96% of the schedule).
    raw_lower = [(0.05, 0.30), (0.52, 0.77)]
    raw_upper = [(0.25, 0.48), (0.70, 0.93)]
    raw_wrap = (0.95, 1.02)   # runs past the right edge; continues at the left
    mapped = []
    for start, end in [*raw_lower, *raw_upper, raw_wrap]:
        if end <= 1.0:
            mapped.append((start, end))
        else:
            mapped.append((start, 1.0))
            mapped.append((0.0, end - 1.0))
    prohibited = _merge_intervals(mapped)
    feasible = _complement_intervals(prohibited)
    for start, end in prohibited:
        axa.add_patch(patches.Rectangle((start, 0.25), end - start, 0.38,
                                        facecolor=RED, edgecolor="none", alpha=0.78))
    wrap_left = (0.0, raw_wrap[1] - 1.0)        # the piece that wrapped to t=0
    wrap_right = (raw_wrap[0], 1.0)             # the tail of the original piece
    for start, end in (wrap_left, wrap_right):
        axa.add_patch(patches.Rectangle((start, 1.25), end - start, 0.38,
                                        facecolor=PURPLE, edgecolor="none", alpha=0.78))
    for start, end in feasible:
        axa.add_patch(patches.Rectangle((start, 2.25), end - start, 0.38,
                                        facecolor=GREEN, edgecolor="none", alpha=0.85))
    axa.text(0.98, 3.08, "Schematic; widths not to scale, no data-derived pair",
             ha="right", va="top", fontsize=7, color=GRAY)
    axa.text(0.5, 0.04, "edge/bad dilation, circularly merged", ha="center",
             fontsize=8.4, color=RED)
    axa.text(0.5, 1.04, "circular wrap-around", ha="center",
             fontsize=8.4, color=PURPLE)
    axa.text(0.5, 2.04, "enumerated feasible intervals", ha="center",
             fontsize=8.4, color=GREEN)

    # The census result itself: every eligible pair's feasible domain is a
    # fraction of a percent of the schedule.
    axb.hist(frac, bins=np.linspace(0.0, 1.0, 21), color=GREEN, edgecolor="white",
             linewidth=0.4)
    axb.axvline(med, color=INK, ls="--", lw=1.0)
    axb.text(med - 0.03, axb.get_ylim()[1] * 0.94, f"median {med:.2f}%",
             ha="right", va="top", fontsize=7.2, color=INK, gid="data-label")
    axb.text(0.04, 0.96, "351/351 pairs nonempty", transform=axb.transAxes,
             ha="left", va="top", fontsize=7.6, color=INK, fontweight="bold")
    axb.set_xlim(0, 1)
    axb.set_xlabel("Feasible-domain fraction of the movie schedule (%)")
    axb.set_ylabel("Number of pairs")
    axb.set_title("The audited feasible domains", loc="center", x=0.50,
                  fontsize=10.5, fontweight="bold", color=INK, pad=10)

    finish(fig, "Figure_6", [[("A", axa), ("B", axb)]])


def figure7():
    """P0-6 event-density perturbation: ordered Mphase FPR decline, flat Mpower.

    Panel A plots per-density FPR point estimates with Clopper-Pearson 95%
    intervals from paired_fpr_contrast.tsv (Mphase at 0.5x/1x/1.5x, Mpower at
    1x/1.5x, where paired rows exist). Panel B is the paired-difference forest
    plot for the prespecified 1.5x-versus-1x contrast with Newcombe method-10
    intervals and the preregistered attenuation threshold.
    """
    P06 = ROOT / "processed" / "subject" / "group" / "ieeg_p0_6_perturbation"
    pc = pd.read_csv(P06 / "paired_fpr_contrast.tsv", sep="\t")

    fig = plt.figure(figsize=(7.48, 3.55))
    gs = GridSpec(1, 2, figure=fig, wspace=0.5)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])

    def density_value(tag):
        return float(tag.replace("p", ".").rstrip("x"))

    # Panel A: per-density FPR with CP intervals. The 1x reference is carried
    # in the FPR_1p0 columns of each paired row and is drawn explicitly so the
    # prespecified 1.5x-versus-1x contrast has both endpoints visible.
    for family, color, alpha, ls in [("Mphase", ORANGE, 1.0, "-"),
                                     ("Mpower", BLUE, 0.6, "--")]:
        sub = pc[pc.family == family]
        for run, off, marker, face in [("R1", -0.045, "o", color), ("R2", 0.045, "s", "white")]:
            pts = {}
            for _, r in sub[sub.run == run].iterrows():
                pts[density_value(r.density)] = (r.FPR_1p5, r.FPR_1p5_CP_low, r.FPR_1p5_CP_high)
                pts[1.0] = (r.FPR_1p0, r.FPR_1p0_CP_low, r.FPR_1p0_CP_high)
            xs = sorted(pts)
            ys = np.array([pts[x][0] for x in xs])
            los = np.array([pts[x][1] for x in xs])
            his = np.array([pts[x][2] for x in xs])
            xo = np.array(xs) + off
            errorbar_from_bounds(axa, xo, ys, los, his, fmt=marker,
                                 color=color, mfc=face, mec=color, alpha=alpha,
                                 capsize=2.5, lw=1.1, markersize=5.0, zorder=3)
            axa.plot(xo, ys, color=color, lw=1.0, ls=ls, alpha=alpha * 0.65, zorder=2)
    axa.axhline(0.05, color=GRAY, lw=0.9, ls="--")
    axa.text(0.30, 0.061, "nominal 0.05", ha="left", fontsize=6.8, color=GRAY,
             gid="data-label")
    axa.set_xticks([0.5, 1.0, 1.5], ["0.5×", "1×", "1.5×"])
    axa.set_xlim(0.28, 1.80)
    axa.set_ylim(0, 0.92)
    axa.set_xlabel("Schedule density relative to the original schedule")
    axa.set_ylabel("False-positive rate (95% CP interval)")
    axa.set_title("Mphase FPR declined with schedule density", loc="center", x=0.50,
                  fontsize=8.5, pad=10)
    route_legend(axa, INK, extra=[Line2D([], [], color=BLUE, ls="--", lw=1.2,
                                         label="Mpower")])

    # Panel B: paired 1.5x-versus-1x differences (Newcombe method-10 intervals);
    # row order mirrors Figure 1 (Mpower above, Mphase below); run encoded by
    # fill state (R1 filled, R2 white).
    rows = ["Mpower R1", "Mpower R2", "Mphase R1", "Mphase R2"]
    fam_run = [tuple(r.split(" ")) for r in rows]
    for i, (family, run) in enumerate(fam_run):
        r = pc[(pc.family == family) & (pc.run == run) & (pc.density == "1p5x")].iloc[0]
        y = len(rows) - 1 - i
        color = ORANGE if family == "Mphase" else BLUE
        axb.scatter(r.Delta, y, s=45, color=color, zorder=4)
        axb.hlines(y, r.CI_lower, r.CI_upper, color=color, lw=1.6, zorder=3)
        axb.vlines([r.CI_lower, r.CI_upper], y - 0.16, y + 0.16, color=color,
                   lw=1.0, zorder=3)
    axb.axvline(0, color=GRAY, lw=0.9)
    axb.axvline(-0.05, color=RED, lw=1.0, ls="--")
    axb.text(-0.058, -0.56, "attenuation\nthreshold", ha="right", va="top",
             fontsize=6.8, color=RED, gid="data-label", linespacing=1.35)
    axb.set_yticks(range(4), rows[::-1])
    axb.set_xlim(-0.48, 0.42)
    axb.set_ylim(-1.0, 3.8)
    axb.set_xlabel("Paired FPR difference (1.5× − 1×)")
    axb.set_title("Prespecified paired contrast", loc="center", x=0.50,
                  fontsize=8.5, pad=10)

    finish(fig, "Figure_S2", [[("A", axa), ("B", axb)]])


def figure8():
    """P0-7 aggregation remedy and injection-matched validation (Mphase).

    Panel A: pure-null FPR under the participant-median versus participant-mean
    reducer (confirmation arm, 120 fresh worlds), Clopper-Pearson 95% intervals.
    Panel B: mean-reducer FPR across the injection-matched arms — random-phase
    control P0 and phase-concentrated P4 (60 matched worlds) plus the
    gain-flattened P0F diagnostic on the same backgrounds.
    Panel C: world-mean observed versus refit-surrogate PPC for the two Stage B
    arms with the identity line; every P4 world sits below identity, showing
    that the 60/60 "detections" were negative-direction.
    """
    P07 = (ROOT / "processed" / "subject" /
           "group" / "ieeg_p0_7_aggregation_remedy")
    fa = pd.read_csv(P07 / "stage_a_confirmation" / "main_fpr.tsv", sep="\t")
    fb = pd.read_csv(P07 / "stage_b_injection_matched" / "main_fpr.tsv", sep="\t")
    fc = pd.read_csv(P07 / "stage_c_gain_flat" / "main_fpr.tsv", sep="\t")
    ppc = pd.read_csv(P07 / "stage_b_injection_matched" / "world_ppc_means.tsv", sep="\t")

    fig = plt.figure(figsize=(7.48, 2.85))
    gs = GridSpec(1, 3, figure=fig, wspace=0.52, width_ratios=[0.9, 1.25, 1.0])
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[0, 2])

    def fpr_points(ax, table, xpos, color, runs=("R1", "R2"), off=0.095):
        for run, sign, marker, face in [("R1", -1, "o", color), ("R2", 1, "s", "white")]:
            r = table[table.run == run]
            assert len(r) == 1
            r = r.iloc[0]
            errorbar_from_bounds(ax, np.array([xpos + sign * off]),
                                 np.array([r.fpr]), np.array([r.cp_low]),
                                 np.array([r.cp_high]), fmt=marker, color=color,
                                 mfc=face, mec=color, capsize=2.5, lw=1.1,
                                 markersize=5.0, zorder=3)

    # Panel A: pure-null confirmation, both reducers on identical pair values.
    for x, reducer, color in [(0.0, "median", GRAY), (1.0, "mean", ORANGE)]:
        fpr_points(axa, fa[fa.reducer == reducer], x, color)
    axa.axhline(0.05, color=GRAY, lw=0.9, ls="--")
    axa.text(-0.38, 0.066, "nominal 0.05", ha="left", fontsize=6.8, color=GRAY,
             gid="data-label")
    axa.set_xticks([0, 1], ["participant\nmedian", "participant\nmean"])
    axa.set_xlim(-0.42, 1.42)
    axa.set_ylim(0, 0.62)
    axa.set_ylabel("False-positive rate (95% CP)")
    axa.set_title("Pure null: mean aggregation\nrestored calibration", loc="center",
                  x=0.50, fontsize=8.5, pad=10)
    route_legend(axa, INK)

    # Panel B: injection-matched arms under the mean reducer.
    arms = [(0.0, fb[(fb.arm == "P0") & (fb.reducer == "mean")], GRAY,
             "P0\nrandom"),
            (1.0, fb[(fb.arm == "P4") & (fb.reducer == "mean")], ORANGE,
             "P4\nκ = 4"),
            (2.0, fc[fc.reducer == "mean"], BLUE, "P0F\nflat gain")]
    for x, table, color, _label in arms:
        fpr_points(axb, table, x, color)
    axb.axhline(0.05, color=GRAY, lw=0.9, ls="--")
    axb.set_xticks([a[0] for a in arms], [a[3] for a in arms])
    axb.set_xlim(-0.45, 2.45)
    axb.set_ylim(0, 1.08)
    axb.set_title("Injection-matched revalidation\nexposed residual miscalibration",
                  loc="center", x=0.50, fontsize=8.5, pad=10)

    # Panel C: world-mean observed vs refit-surrogate PPC, identity line.
    lim = (-0.015, 0.055)
    axc.plot(lim, lim, color=GRAY, lw=0.9, ls="--", zorder=1)
    for arm, color in [("P0", GRAY), ("P4", ORANGE)]:
        a = ppc[ppc.arm == arm]
        axc.scatter(a.surrogate_ppc_mean, a.observed_ppc_mean, s=9, color=color,
                    alpha=0.55, lw=0, zorder=2, rasterized=True)
    axc.text(-0.0125, 0.044, "P4: 60/60 worlds\nbelow identity", fontsize=6.8,
             color=ORANGE, gid="data-label", va="top", linespacing=1.35)
    axc.set_xlim(lim)
    axc.set_ylim(lim)
    axc.set_xlabel("Refit-surrogate PPC (world mean)")
    axc.set_ylabel("Observed PPC (world mean)")
    axc.set_title("Phase-concentrated arm rejected\nbelow its surrogates", loc="center",
                  x=0.50, fontsize=8.5, pad=10)

    finish(fig, "Figure_3", [[("A", axa), ("B", axb), ("C", axc)]])


def main():
    setup()
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    figure6()
    figure7()
    figure8()


if __name__ == "__main__":
    main()
