"""Plot the jet-radius metrics, one curve per model.

Reads JetRadius.pkl (from scan.py) in each model directory matching the
pattern; /store/... patterns are fetched from EOS. Each point is the mean over
truth jets.

  * <metric>.pdf: metric vs visible-jet radius R, at one truth radius R_DH
  * <score>_ptbins_<model>.pdf: each score vs R in bins of the truth jet's
    visible pT, and best_radius_vs_pt_<score>.pdf, the best R per pT bin
  * if the scans have several truth radii, for quality and iou:
    <metric>_2d_<model>.pdf, the mean over (R_DH, R) with the best R per R_DH
    marked, and best_radius_<metric>.pdf, the best-R stars of every model
    overlaid on the same axes

Also prints the best R per model and truth radius.

    python3 JetRadiusOptimization/plot.py -m "models/*" -x rinv
"""

import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
from magiconfig import ArgumentParser, ArgumentDefaultsRawHelpFormatter

import radius  # puts the repository root on sys.path
from common import resolve_models, set_plot_style

LABELS = {
    "quality": "quality = acceptance x purity",
    "iou": "$p_T$-weighted IoU with the truth jet",
    "acceptance": "acceptance ($p_T$ fraction in complete dark hadrons)",
    "purity": "purity ($f_{full}$)",
    "f_partial": "$p_T$ fraction from partly contained dark hadrons",
    "f_nodh": "$p_T$ fraction with no dark-hadron ancestor",
    "matched": "fraction of truth jets matched",
}


SCORES = ("quality", "iou", "acceptance", "purity")


def mean_grid(scan, metric, pt_range=(0, np.inf)):
    """Mean of a metric over truth jets in a visible-pT range: array[truth radius, radius].

    NaN values are skipped; a cell with no values is NaN.
    """
    m = scan["metrics"]
    n_r = len(scan["radii"])
    cell = np.searchsorted(scan["truth_radii"], m["truth_radius"]) * n_r + np.searchsorted(scan["radii"], m["radius"])
    values = m[metric]
    keep = (m["truth_pt"] >= pt_range[0]) & (m["truth_pt"] < pt_range[1]) & ~np.isnan(values)
    size = len(scan["truth_radii"]) * n_r
    total = np.bincount(cell[keep], weights=values[keep], minlength=size)
    count = np.bincount(cell[keep], minlength=size)
    with np.errstate(invalid="ignore"):
        return (total / count).reshape(-1, n_r)


def square_axes(ax, scan, line_color, pad=0.0):
    """Visible radius on x, truth radius on y: same range, ticks and scale, with R = R_DH.

    pad: extra margin on both axes, e.g. for points spread within their row
    """
    half = [np.min(np.diff(r)) / 2 if len(r) > 1 else 0.1 for r in (scan["radii"], scan["truth_radii"])]
    lo = max(scan["radii"][0] - half[0], scan["truth_radii"][0] - half[1])
    hi = min(scan["radii"][-1] + half[0], scan["truth_radii"][-1] + half[1])
    lo, hi = lo - pad, hi + pad
    ax.plot([lo, hi], [lo, hi], ":", color=line_color, lw=1)
    ticks = np.arange(np.ceil(lo * 5) / 5, hi + 1e-9, 0.2)
    ax.set(xlim=(lo, hi), ylim=(lo, hi), xticks=ticks, yticks=ticks, aspect="equal",
           xlabel="visible jet radius R", ylabel="truth radius $R_{DH}$")


def legend(ax):
    """Inside the axes for a few curves; below them, in columns, for many."""
    if len(ax.get_legend_handles_labels()[1]) <= 5:
        ax.legend(fontsize="small")
    else:
        ax.legend(fontsize="x-small", ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.18))


def save(fig, name):
    fig.savefig(os.path.join(args.dir, f"{name}.{args.format}"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsRawHelpFormatter)
    parser.add_argument("-m", "--models", required=True, help="pattern of model directories")
    parser.add_argument("-x", type=str, default="rinv", help="model parameter used to label the curves")
    parser.add_argument("-y", type=str, nargs="+", default=list(LABELS), choices=list(LABELS), help="metrics to plot vs R")
    parser.add_argument("--truth-radius", type=float, default=0.8, help="truth radius R_DH for the plots vs R")
    parser.add_argument("--pt-bins", type=float, nargs="+", default=[0, 50, 100, 200, np.inf], help="edges of the truth-jet visible pT bins [GeV]")
    parser.add_argument("--dir", type=str, default="JetRadius_plots", help="output directory")
    parser.add_argument("--format", type=str, default="pdf", help="plot file format")
    args = parser.parse_args()

    scans = []
    for model_dir in resolve_models(args.models):
        with open(os.path.join(model_dir, "JetRadius.pkl"), "rb") as f:
            scans.append(pickle.load(f))
    scans.sort(key=lambda s: s["model"][args.x])

    def label(scan):
        return f"{args.x} = {scan['model'][args.x]}"

    set_plot_style()
    os.makedirs(args.dir, exist_ok=True)

    for metric in args.y:
        fig, ax = plt.subplots()
        for scan in scans:
            t = scan["truth_radii"].index(args.truth_radius)
            ax.plot(scan["radii"], mean_grid(scan, metric)[t], marker="o", label=label(scan))
        ax.set_xlabel("visible jet radius R")
        ax.set_ylabel(f"mean {LABELS[metric]}")
        ax.set_title(f"truth radius $R_{{DH}}$ = {args.truth_radius}", fontsize=14)
        legend(ax)
        save(fig, metric)

    pt_bins = list(zip(args.pt_bins[:-1], args.pt_bins[1:]))
    pt_label = [f"{lo:g}-{hi:g} GeV" if np.isfinite(hi) else f"> {lo:g} GeV" for lo, hi in pt_bins]
    for metric in SCORES:
        fig_best, ax_best = plt.subplots()
        for scan in scans:
            t = scan["truth_radii"].index(args.truth_radius)
            fig, ax = plt.subplots()
            best = []
            for pt_range, name in zip(pt_bins, pt_label):
                curve = mean_grid(scan, metric, pt_range)[t]
                ax.plot(scan["radii"], curve, marker="o", label=f"truth $p_T$ {name}")
                best.append(scan["radii"][np.nanargmax(curve)] if np.isfinite(curve).any() else np.nan)
            ax.set_xlabel("visible jet radius R")
            ax.set_ylabel(f"mean {LABELS[metric]}")
            ax.set_title(f"{label(scan)}, $R_{{DH}}$ = {args.truth_radius}", fontsize=14)
            legend(ax)
            save(fig, f"{metric}_ptbins_{args.x}_{scan['model'][args.x]}")
            ax_best.plot(pt_label, best, marker="o", label=label(scan))
        ax_best.set_xlabel("truth jet visible $p_T$")
        ax_best.set_ylabel(f"best visible jet radius R ({metric})")
        ax_best.set_title(f"$R_{{DH}}$ = {args.truth_radius}", fontsize=14)
        legend(ax_best)
        save(fig_best, f"best_radius_vs_pt_{metric}")

    multi = [s for s in scans if len(s["truth_radii"]) > 1]
    for metric in ("quality", "iou") if multi else ():
        fig_best, ax_best = plt.subplots()
        colors = plt.cm.viridis(np.linspace(0, 0.9, len(multi)))
        # spread the models vertically within each R_DH row so equal best-R points stay visible
        step = np.min(np.diff(multi[0]["truth_radii"]))
        nudge = (np.arange(len(multi)) - (len(multi) - 1) / 2) * 0.6 * step / max(len(multi), 1)
        for scan, color, dy in zip(multi, colors, nudge):
            grid = mean_grid(scan, metric)
            best = np.array(scan["radii"])[np.nanargmax(grid, axis=1)]
            fig, ax = plt.subplots()
            image = ax.pcolormesh(scan["radii"], scan["truth_radii"], grid, shading="nearest")
            ax.plot(best, scan["truth_radii"], "w*", ms=10)
            square_axes(ax, scan, "w")
            ax.set_title(f"{label(scan)} (star = best R)", fontsize=14)
            fig.colorbar(image, ax=ax, label=f"mean {metric}")
            save(fig, f"{metric}_2d_{args.x}_{scan['model'][args.x]}")
            ax_best.plot(best, np.array(scan["truth_radii"]) + dy, "*-", color=color, ms=9, lw=1, label=label(scan))
            for t, r, q in zip(scan["truth_radii"], best, np.nanmax(grid, axis=1)):
                print(f"{label(scan)}, R_DH = {t}: best R by {metric} = {r} (mean {q:.3f})")
        square_axes(ax_best, multi[0], "k", pad=0.06)
        ax_best.set_title(f"best R by {metric}" + (" (rows spread for visibility)" if len(multi) > 1 else ""), fontsize=12)
        legend(ax_best)
        save(fig_best, f"best_radius_{metric}")
