"""Plot the jet-radius metrics, one curve per model.

Reads JetRadius.pkl (from scan.py) in each model directory matching the
pattern; /store/... patterns are fetched from EOS. Each point is the mean over
truth jets.

  * <metric>.pdf: metric vs visible-jet radius R, at one truth radius R_DH
  * <score>_ptbins_<model>.pdf: each score vs R in bins of the truth jet's
    visible pT, and best_radius_vs_pt_<score>.pdf, the best R per pT bin
  * <quantity>_dist_<model>.pdf: distributions at several radii R (--dist-radii),
    truth_radius_dist_<model>.pdf: the truth jets' 90/95/99% containment radii,
    and truth_r<pct>_by_rdh_<model>.pdf: each of them for several truth radii
    R_DH (--dist-truth-radii)
  * truth_radius_vs_rdh.pdf: radii containing 90/95/99% of the truth jets'
    visible pT, vs R_DH
  * if the scans have several truth radii, for quality and iou:
    <metric>_2d_<model>.pdf, the mean over (R_DH, R) with the best R per R_DH
    marked, and best_radius_<metric>.pdf, the best-R stars of every model
    overlaid on the same axes; both also show the truth jets' 90/95/99% radii
    (only the 95% one when several models are overlaid)
  * with --displays, for models with a local events.root (first
    --display-max-events events), at --truth-radius:
    event_<i>_<model>.pdf, one event per truth-jet pT bin: the visible
    particles around the leading truth jet's axis by origin, its truth
    90/95/99% radii, and the matched visible jet at each --display-radii R
    with its own 95% radius; stacked_<model>[_pt<bin>].pdf, all truth jets
    stacked around their axes, one image per origin (this truth jet's dark
    hadrons, other dark hadrons, no dark-hadron ancestor), and their pT vs
    DeltaR stacked by origin, inclusive and per pT bin

Also prints the best R per model and truth radius.

    python3 JetRadiusOptimization/plot.py -m "models/*" -x rinv
"""

import os
import pickle

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
from magiconfig import ArgumentParser, ArgumentDefaultsRawHelpFormatter
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle

import radius
from common import resolve_models, set_plot_style
from darkancestry import dark_hadron_ancestry

LABELS = {
    "quality": "quality = acceptance x purity",
    "iou": "$p_T$-weighted IoU with the truth jet",
    "acceptance": "acceptance ($p_T$ fraction in complete dark hadrons)",
    "purity": "purity ($f_{full}$)",
    "f_partial": "$p_T$ fraction from partly contained dark hadrons",
    "f_nodh": "$p_T$ fraction with no dark-hadron ancestor",
    "matched": "fraction of truth jets matched",
    "psi_truth": "$\\Psi_T(R)$: truth visible $p_T$ fraction within R",
    "jet_r95": "jet $\\Delta R_{95}$",
    "jet_girth": "jet girth",
    "jet_ptd": "jet $p_T^D$",
    "jet_major": "jet major axis",
    "jet_minor": "jet minor axis",
    "jet_nconst": "jet constituents",
}
DIST = ("quality", "iou", "acceptance", "purity", "psi_truth",
        "jet_r95", "jet_girth", "jet_ptd", "jet_major", "jet_minor", "jet_nconst")
# line style of each truth containment radius (percent of visible pT)
TRUTH_STYLES = {90: ":", 95: "--", 99: "-."}
# particle origin relative to one truth jet: name and colour (displays)
ORIGINS = (("this truth jet's dark hadrons", "tab:red", "Reds"),
           ("other dark hadrons", "tab:blue", "Blues"),
           ("no dark-hadron ancestor (ISR/UE)", "0.6", "Greys"))
DISPLAY_LIM = 2.5   # half-width of the displays in eta and phi
# smaller fonts than set_plot_style for the multi-panel displays
DISPLAY_STYLE = {"font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12, "legend.fontsize": 8,
                 "xtick.labelsize": 10, "ytick.labelsize": 10}


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
    ax.plot([lo, hi], [lo, hi], "-", color=line_color, lw=0.8, alpha=0.5, label="$R = R_{DH}$")
    ticks = np.arange(np.ceil(lo * 5) / 5, hi + 1e-9, 0.2)
    ax.set(xlim=(lo, hi), ylim=(lo, hi), xticks=ticks, yticks=ticks, aspect="equal",
           xlabel="visible jet radius R", ylabel="truth radius $R_{DH}$")


def legend(ax):
    """Inside the axes for a few curves; below them, in columns, for many."""
    if len(ax.get_legend_handles_labels()[1]) <= 5:
        ax.legend(fontsize="small")
    else:
        ax.legend(fontsize="x-small", ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.18))


def truth_radius(scan, pct):
    """Mean radius containing pct% of the truth jets' visible pT, per truth radius."""
    m = scan["metrics"]
    first = m["radius"] == scan["radii"][0]   # one row per truth jet
    return np.array([np.nanmean(m[f"truth_r{pct}"][first & (m["truth_radius"] == t)]) for t in scan["truth_radii"]])


def wrap(dphi):
    return (dphi + np.pi) % (2 * np.pi) - np.pi


def percentile_radius(dr, pt, pct):
    """Smallest DeltaR holding pct% of the pT (as radius_percentile in Histogram.py)."""
    order = np.argsort(dr)
    cumul = np.cumsum(pt[order])
    return dr[order][np.argmax(cumul >= pct / 100 * cumul[-1])]


def circle(ax, center, r, label=None, **kw):
    """Circle with a line (not a box) in the legend."""
    ax.add_patch(Circle(center, r, fill=False, **kw))
    ax.plot([], [], color=kw.get("color"), ls=kw.get("ls", "-"), lw=kw.get("lw", 1), label=label)


def truth_jet_views(events, masks, truth, i):
    """Per truth jet of event i: its particles' (deta, dphi) around the dark-hadron jet axis,
    origin (index into ORIGINS), the truth jet's visible pT and its 90/95/99% radii."""
    mask = ak.to_numpy(masks[i])
    descends = ((mask[None, :] >> np.arange(len(events.DarkHadronCandidate[i]), dtype=np.uint64)[:, None])
                & np.uint64(1)).astype(bool)
    gen, dark = events.GenCandidate[i], events.DarkHadronCandidate[i]
    pt = ak.to_numpy(gen.pt).astype(float)
    eta, phi = ak.to_numpy(gen.eta).astype(float), ak.to_numpy(gen.phi).astype(float)
    views = []
    for t, dh in enumerate(ak.to_numpy(j) for j in truth[i]):
        in_truth = descends[dh].any(axis=0)
        if not in_truth.any():
            continue
        px, py, pz = (ak.to_numpy(getattr(dark, c))[dh].sum() for c in ("px", "py", "pz"))
        axis = (np.arcsinh(pz / np.hypot(px, py)), np.arctan2(py, px))
        deta, dphi = eta - axis[0], wrap(phi - axis[1])
        dr = np.hypot(deta, dphi)
        views.append({"truth_jet": t, "dh": dh, "axis": axis, "deta": deta, "dphi": dphi, "pt": pt,
                      "origin": np.where(in_truth, 0, np.where(mask != 0, 1, 2)),
                      "truth_pt": pt[in_truth].sum(),
                      "truth_r": {p: percentile_radius(dr[in_truth], pt[in_truth], p) for p in radius.SHAPE_PCTS},
                      "descends": descends, "no_dh": mask == 0, "eta": eta, "phi": phi})
    return views


def event_display(events, masks, truth, visible, i, title, name):
    """Leading truth jet of event i: truth containment (left) and matched visible jets (right)."""
    views = truth_jet_views(events, masks, truth, i)
    v = max(views, key=lambda x: x["truth_pt"])
    deta, dphi, pt, origin = v["deta"], v["dphi"], v["pt"], v["origin"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.3))
    near = np.hypot(deta, dphi) < 1.5 * DISPLAY_LIM
    for ax in axes:
        for k in (2, 1, 0):
            sel = near & (origin == k)
            ax.scatter(deta[sel], dphi[sel], s=3 + 4 * pt[sel], c=ORIGINS[k][1], alpha=0.6, lw=0, label=ORIGINS[k][0])
        ax.plot(0, 0, "k+", ms=12, mew=2)
        ax.set(xlim=(-DISPLAY_LIM, DISPLAY_LIM), ylim=(-DISPLAY_LIM, DISPLAY_LIM), aspect="equal",
               xlabel="$\\Delta\\eta$ to dark-hadron jet axis", ylabel="$\\Delta\\phi$ to dark-hadron jet axis")

    ax = axes[0]
    dark = events.DarkHadronCandidate[i]
    ax.plot(ak.to_numpy(dark.eta)[v["dh"]] - v["axis"][0], wrap(ak.to_numpy(dark.phi)[v["dh"]] - v["axis"][1]),
            "kx", ms=7, mew=1.5, label=f"dark hadrons ({len(v['dh'])})")
    circle(ax, (0, 0), args.truth_radius, color="0.3", label=f"$R_{{DH}}$ = {args.truth_radius}")
    for p, style in TRUTH_STYLES.items():
        circle(ax, (0, 0), v["truth_r"][p], color="k", ls=style, lw=1.6,
               label=f"truth $\\Delta R_{{{p}}}$ = {v['truth_r'][p]:.2f}")
    ax.set_title(f"truth jet: visible $p_T$ = {v['truth_pt']:.0f} GeV", fontsize=12)
    ax.legend(loc="upper left", framealpha=0.8)

    ax = axes[1]
    circle(ax, (0, 0), v["truth_r"][95], color="k", ls="--", lw=1.2,
           label=f"truth $\\Delta R_{{95}}$ = {v['truth_r'][95]:.2f}")
    truth_dh = [ak.to_numpy(j) for j in truth[i]]
    gen = events.GenCandidate[i]
    for r, color in zip(args.display_radii, plt.cm.viridis(np.linspace(0, 0.85, len(args.display_radii)))):
        jets = [ak.to_numpy(j) for j in visible[r][i]]
        rows = radius.event_metrics(v["descends"], pt, v["no_dh"], truth_dh, jets, 0.95)
        row = next(x for x in rows if x["truth_jet"] == v["truth_jet"])
        if row["jet"] < 0:
            continue
        c = jets[row["jet"]]
        px, py, pz = (ak.to_numpy(getattr(gen, k))[c].sum() for k in ("px", "py", "pz"))
        j_eta, j_phi = np.arcsinh(pz / np.hypot(px, py)), np.arctan2(py, px)
        center = (j_eta - v["axis"][0], wrap(j_phi - v["axis"][1]))
        j_r95 = percentile_radius(np.hypot(v["eta"][c] - j_eta, wrap(v["phi"][c] - j_phi)), pt[c], 95)
        circle(ax, center, r, color=color, lw=2, label=f"R = {r}: quality {row['quality']:.2f}, IoU {row['iou']:.2f}")
        circle(ax, center, j_r95, color=color, ls="--", lw=1.2, label=f"   jet $\\Delta R_{{95}}$ = {j_r95:.2f}")
        ax.plot(*center, "o", color=color, ms=5)
    ax.set_title("matched visible jets: R (solid), jet $\\Delta R_{95}$ (dashed)", fontsize=12)
    ax.legend(loc="upper left", framealpha=0.8)
    fig.suptitle(f"{title}, event {i}, truth jet at $R_{{DH}}$ = {args.truth_radius}", fontsize=13)
    fig.tight_layout()
    save(fig, name)


def stacked_display(views, title, name):
    """All truth jets stacked around their axes, each normalised to its own visible pT."""
    if not views:
        return
    keep = [np.hypot(v["deta"], v["dphi"]) < 2 * DISPLAY_LIM for v in views]
    deta, dphi, origin = (np.concatenate([v[k][s] for v, s in zip(views, keep)]) for k in ("deta", "dphi", "origin"))
    w = np.concatenate([v["pt"][s] / v["truth_pt"] for v, s in zip(views, keep)]) / len(views)
    mean_r = {p: np.mean([v["truth_r"][p] for v in views]) for p in radius.SHAPE_PCTS}
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(args.display_radii)))

    fig = plt.figure(figsize=(18, 12))
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 0.85])
    bins = np.linspace(-DISPLAY_LIM, DISPLAY_LIM, 101)
    for k, (origin_label, _, cmap) in enumerate(ORIGINS):
        ax = fig.add_subplot(grid[0, k])
        sel = origin == k
        _, _, _, image = ax.hist2d(deta[sel], dphi[sel], bins=bins, weights=w[sel], norm=LogNorm(), cmap=cmap)
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03,
                     label=f"mean fraction of truth $p_T$ per bin (total {w[sel].sum():.2f})")
        for p, style in TRUTH_STYLES.items():
            circle(ax, (0, 0), mean_r[p], color="k", ls=style, lw=1.8, label=f"mean truth $\\Delta R_{{{p}}}$ = {mean_r[p]:.2f}")
        for r, color in zip(args.display_radii, colors):
            circle(ax, (0, 0), r, color=color, lw=1.5, label=f"R = {r}")
        ax.set(xlim=(-DISPLAY_LIM, DISPLAY_LIM), ylim=(-DISPLAY_LIM, DISPLAY_LIM), aspect="equal",
               xlabel="$\\Delta\\eta$ to dark-hadron jet axis", ylabel="$\\Delta\\phi$ to dark-hadron jet axis")
        ax.set_title(origin_label, fontsize=12)
        if k == 0:
            ax.legend(loc="upper left", framealpha=0.8)

    ax = fig.add_subplot(grid[1, :])
    dr = np.hypot(deta, dphi)
    rbins = np.linspace(0, DISPLAY_LIM, 51)
    ax.hist([dr[origin == k] for k in range(3)], bins=rbins, weights=[w[origin == k] for k in range(3)],
            stacked=True, color=[o[1] for o in ORIGINS], alpha=0.8, label=[o[0] for o in ORIGINS])
    for p, style in TRUTH_STYLES.items():
        ax.axvline(mean_r[p], color="k", ls=style, lw=1.8, label=f"mean truth $\\Delta R_{{{p}}}$ = {mean_r[p]:.2f}")
    for r, color in zip(args.display_radii, colors):
        ax.axvline(r, color=color, lw=1.5, label=f"R = {r}")
    ax_psi = ax.twinx()   # cumulative fraction of the truth jets' own pT (Psi)
    own = np.histogram(dr[origin == 0], bins=rbins, weights=w[origin == 0])[0]
    ax_psi.plot(rbins[1:], np.cumsum(own) / own.sum(), "r-", lw=1.2)
    ax_psi.set_ylim(0, 1.02)
    ax_psi.set_ylabel("cumulative fraction of truth visible $p_T$ (red line)", color="tab:red")
    ax.set(xlim=(0, DISPLAY_LIM), xlabel="$\\Delta R$ to dark-hadron jet axis",
           ylabel="mean $p_T$ per bin / truth visible $p_T$")
    ax.set_title("$p_T$ vs distance from the dark-hadron axis, stacked by origin", fontsize=12)
    ax.legend(loc="upper right", bbox_to_anchor=(0.97, 0.75), framealpha=0.9)
    fig.suptitle(f"{title}: {len(views)} truth jets at $R_{{DH}}$ = {args.truth_radius}", fontsize=13)
    fig.tight_layout()
    save(fig, name)


def displays(model_dir, title, name, pt_bins, pt_label):
    """Event displays and stacked truth jets for one model, from its events.root."""
    with plt.rc_context(DISPLAY_STYLE):
        _displays(model_dir, title, name, pt_bins, pt_label)


def _displays(model_dir, title, name, pt_bins, pt_label):
    events = radius.load(os.path.join(model_dir, "events.root"))[: args.display_max_events]
    masks = dark_hadron_ancestry(events)
    truth = radius.cluster(events.DarkHadronCandidate, args.truth_radius)
    visible = {r: radius.cluster(events.GenCandidate, r) for r in args.display_radii}
    views = [truth_jet_views(events, masks, truth, i) for i in range(len(events))]
    leading = np.array([max((v["truth_pt"] for v in vs), default=0) for vs in views])
    for (lo, hi), pt_name in zip(pt_bins, pt_label):
        in_bin = np.flatnonzero((leading >= lo) & (leading < hi))
        if len(in_bin):
            event_display(events, masks, truth, visible, int(in_bin[0]), f"{title}, leading truth $p_T$ {pt_name}",
                          f"event_{in_bin[0]}_{name}")
    flat = [v for vs in views for v in vs]
    stacked_display(flat, title, f"stacked_{name}")
    for (lo, hi), pt_name in zip(pt_bins, pt_label):
        stacked_display([v for v in flat if lo <= v["truth_pt"] < hi], f"{title}, truth $p_T$ {pt_name}",
                        f"stacked_{name}_pt{lo:g}-{hi:g}")


def save(fig, name):
    fig.savefig(os.path.join(args.dir, f"{name}.{args.format}"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsRawHelpFormatter)
    parser.add_argument("-m", "--models", required=True, help="pattern of model directories")
    parser.add_argument("-x", type=str, default="rinv", help="model parameter used to label the curves")
    parser.add_argument("-y", type=str, nargs="+", default=list(LABELS), choices=list(LABELS), help="metrics to plot vs R")
    parser.add_argument("--truth-radius", type=float, default=0.8, help="truth radius R_DH for the plots vs R")
    parser.add_argument("--dist-truth-radii", type=float, nargs="+", default=[0.4, 0.8, 1.2, 1.6, 2.0], help="truth radii R_DH for the truth-radius distributions")
    parser.add_argument("--dist-radii", type=float, nargs="+", default=[0.4, 0.8, 1.2, 1.6, 2.0], help="radii R for the distributions")
    parser.add_argument("--pt-bins", type=float, nargs="+", default=[0, 50, 100, 200, np.inf], help="edges of the truth-jet visible pT bins [GeV]")
    parser.add_argument("--displays", action="store_true", help="also make event displays and stacked truth jets (needs events.root)")
    parser.add_argument("--display-radii", type=float, nargs="+", default=[0.4, 0.8, 1.2, 1.6], help="radii R drawn in the displays")
    parser.add_argument("--display-max-events", type=int, default=1000, help="events read for the displays")
    parser.add_argument("--dir", type=str, default="JetRadius_plots", help="output directory")
    parser.add_argument("--format", type=str, default="pdf", help="plot file format")
    args = parser.parse_args()

    scans = []
    for model_dir in resolve_models(args.models):
        with open(os.path.join(model_dir, "JetRadius.pkl"), "rb") as f:
            scans.append(dict(pickle.load(f), dir=model_dir))
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

    for scan in scans:
        m, name = scan["metrics"], f"{args.x}_{scan['model'][args.x]}"
        at_truth = m["truth_radius"] == args.truth_radius
        for q in DIST:
            values = m[q][at_truth & np.isfinite(m[q])]
            if not len(values):
                continue
            bins = np.linspace(0, np.percentile(values, 99.5), 41) if q not in DIST[:5] else np.linspace(0, 1, 41)
            fig, ax = plt.subplots()
            for r in args.dist_radii:
                ax.hist(m[q][at_truth & (m["radius"] == r)], bins=bins, histtype="step", density=True, label=f"R = {r}")
            ax.set_xlabel(LABELS[q])
            ax.set_ylabel("truth jets (normalized)")
            ax.set_title(f"{label(scan)}, $R_{{DH}}$ = {args.truth_radius}", fontsize=14)
            legend(ax)
            save(fig, f"{q}_dist_{name}")
        fig, ax = plt.subplots()
        first = at_truth & (m["radius"] == scan["radii"][0])
        for pct in radius.SHAPE_PCTS:
            ax.hist(m[f"truth_r{pct}"][first], bins=np.linspace(0, 3, 61), histtype="step", density=True, label=f"$\\Delta R_{{{pct}}}$")
        ax.set_xlabel("radius around the dark-hadron jet axis containing x% of its visible $p_T$")
        ax.set_ylabel("truth jets (normalized)")
        ax.set_title(f"{label(scan)}, $R_{{DH}}$ = {args.truth_radius}", fontsize=14)
        legend(ax)
        save(fig, f"truth_radius_dist_{name}")
        first = m["radius"] == scan["radii"][0]   # one row per truth jet
        for pct in radius.SHAPE_PCTS:
            fig, ax = plt.subplots()
            for t in args.dist_truth_radii:
                if t in scan["truth_radii"]:
                    ax.hist(m[f"truth_r{pct}"][first & (m["truth_radius"] == t)], bins=np.linspace(0, 3, 61),
                            histtype="step", density=True, label=f"$R_{{DH}}$ = {t}")
            ax.set_xlabel(f"truth $\\Delta R_{{{pct}}}$: radius containing {pct}% of the visible $p_T$")
            ax.set_ylabel("truth jets (normalized)")
            ax.set_title(label(scan), fontsize=14)
            legend(ax)
            save(fig, f"truth_r{pct}_by_rdh_{name}")

    fig, ax = plt.subplots()
    for scan, color in zip(scans, plt.cm.viridis(np.linspace(0, 0.9, len(scans)))):
        for pct, style in TRUTH_STYLES.items():
            ax.plot(scan["truth_radii"], truth_radius(scan, pct), style, marker="o", ms=3, color=color,
                    label=f"{label(scan)}, {pct}%" if len(scans) == 1 else None)
        if len(scans) > 1:
            ax.plot([], [], "-", color=color, label=label(scan))
    if len(scans) > 1:   # line styles once, colours per model
        for pct, style in TRUTH_STYLES.items():
            ax.plot([], [], style, color="k", label=f"$\\Delta R_{{{pct}}}$")
    ax.set_xlabel("truth radius $R_{DH}$")
    ax.set_ylabel("mean truth $\\Delta R_{90/95/99}$")
    legend(ax)
    save(fig, "truth_radius_vs_rdh")

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
            for pct, style in TRUTH_STYLES.items():
                ax.plot(truth_radius(scan, pct), scan["truth_radii"], style, color="w", lw=1.5,
                        label=f"truth $\\Delta R_{{{pct}}}$")
            square_axes(ax, scan, "w")
            # upper left (small R, large R_DH) holds no stars or lines; dark frame for the white lines
            ax.legend(fontsize="x-small", loc="upper left", facecolor="0.2", labelcolor="w", framealpha=0.7)
            ax.set_title(f"{label(scan)} (star = best R)", fontsize=14)
            fig.colorbar(image, ax=ax, label=f"mean {metric}")
            save(fig, f"{metric}_2d_{args.x}_{scan['model'][args.x]}")
            ax_best.plot(best, np.array(scan["truth_radii"]) + dy, "*-", color=color, ms=9, lw=1, label=label(scan))
            for pct, style in TRUTH_STYLES.items() if len(multi) == 1 else [(95, TRUTH_STYLES[95])]:
                ax_best.plot(truth_radius(scan, pct), scan["truth_radii"], style, color=color, lw=1,
                             label=f"truth $\\Delta R_{{{pct}}}$" if len(multi) == 1 else None)
            for t, r, q in zip(scan["truth_radii"], best, np.nanmax(grid, axis=1)):
                print(f"{label(scan)}, R_DH = {t}: best R by {metric} = {r} (mean {q:.3f})")
        square_axes(ax_best, multi[0], "k", pad=0.06)
        truth_lines = "truth $\\Delta R_{90/95/99}$ (dotted/dashed/dash-dot)" if len(multi) == 1 else "truth $\\Delta R_{95}$ (dashed)"
        ax_best.set_title(f"best R by {metric} (stars), {truth_lines}" + (", rows spread" if len(multi) > 1 else ""), fontsize=11)
        if len(multi) == 1:
            ax_best.legend(fontsize="x-small", loc="upper left")
        else:
            legend(ax_best)
        save(fig_best, f"best_radius_{metric}")

    for scan in scans if args.displays else ():
        if os.path.exists(os.path.join(scan["dir"], "events.root")):
            displays(scan["dir"], label(scan), f"{args.x}_{scan['model'][args.x]}", pt_bins, pt_label)
        else:
            print(f"{scan['dir']}: no events.root, no displays")
