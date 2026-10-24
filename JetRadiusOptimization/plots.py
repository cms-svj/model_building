#!/usr/bin/env python3
"""Stage 4: turn the merged metric tables into the decision figures.

Reads every ``metrics_*.parquet`` from ``scan.py`` plus its ``.meta.json``
sidecar, and produces the plots that actually answer "which radius, and is one
radius defensible".

Figures
-------
``optimal_radius_map.png``      THE money plot. One panel, colour = argmax_R.
                                Flat means one radius works everywhere. A
                                gradient means the optimum tracks a physical
                                scale and you should be fitting a formula.

``regret_map.png``              Cost of the compromise at the chosen global R,
                                on the same axes. Near zero everywhere means
                                the recommendation holds; wherever it is not is
                                what has to be caveated.

``quality_vs_radius.png``       One curve per model point. This is where the
                                *plateau* lives -- whether the peak is sharp or
                                whether 0.8-1.0 are indistinguishable is
                                invisible in any heatmap, and it is the
                                difference between "R = 0.9" and "anything from
                                0.8 to 1.0, we suggest 0.8".

``decomposition_maps.png``      Why a corner fails: acceptance,
                                frac_partial_pt, frac_nodh_pt, and matched
                                fraction on the same axes at one radius.

``regret_by_radius.png``        The per-radius panel grid, plotting REGRET not
                                raw quality. Regret is normalised per model
                                point, so the colour scale is comparable across
                                panels and cells; raw quality is not, and its
                                scale gets consumed by "which point is easy"
                                rather than "which radius is good".

Two things to keep in mind when reading these:

  * ``mmed`` is a hidden third axis and it drives the boost harder than either
    of the plotted ones. Fix it per figure (``--fixed``) or the map is a
    projection averaging over it.
  * high-``rinv`` cells go sparse -- fewer visible constituents means jets fall
    below ``pt_min`` and matching drops -- so read the matched-fraction panel
    before trusting the top of that axis.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import awkward as ak  # noqa: E402
from optimize import choose_radius  # noqa: E402

MIN_JETS_PER_CELL = 50


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metrics", nargs="+", required=True,
                   help="metrics_*.parquet files, or globs.")
    p.add_argument("--outdir", type=Path, default=Path("radius_decision_plots"))
    p.add_argument("--x", default="rinv", help="Parameter on the x axis.")
    p.add_argument("--y", default="mpi", help="Parameter on the y axis.")
    p.add_argument("--metric", default="quality",
                   help="Column to optimise. Default 'quality' = A * P.")
    p.add_argument("--extra-curves", default="",
                   help="Comma-separated extra columns (e.g. "
                        "'frac_full_pt,iou_pt') to render as their own "
                        "one-curve-per-model-point plot, <column>_vs_radius.png, "
                        "against the same R* chosen from --metric. Any column "
                        "radius_metrics.scalars() emits works.")
    p.add_argument("--curves-only", action="store_true",
                   help="Skip the four 2D maps (optimal_radius_map, regret_map, "
                        "decomposition_maps, regret_by_radius) -- useful when "
                        "only one axis of the model-point grid actually varies, "
                        "e.g. a pure mass-point scan at fixed mmed/rinv.")
    p.add_argument("--fixed", nargs="*", default=[],
                   help="param=value filters, e.g. mmed=1000. Without these "
                        "the map averages over the unplotted axes.")
    p.add_argument("--quantile", type=float, default=1.0,
                   help="1.0 = strict minimax regret; 0.9 lets the worst ~10%% "
                        "of model points off the hook. Run both.")
    p.add_argument("--tolerance", type=float, default=0.01)
    return p.parse_args()


def load(patterns, fixed, metric):
    """Collect per-(model point, radius) means plus the axis parameters."""
    paths = sorted({p for pat in patterns for p in glob.glob(pat)})
    if not paths:
        raise SystemExit(f"no metrics files matched: {patterns}")

    filters = {}
    for item in fixed:
        key, _, value = item.partition("=")
        filters[key] = float(value)

    points, quality, matched = {}, {}, {}
    for path in paths:
        meta_path = Path(path).with_suffix(".meta.json")
        if not meta_path.is_file():
            print(f"[WARN] no sidecar for {path}, skipping", flush=True)
            continue
        meta = json.loads(meta_path.read_text())
        params = meta.get("params", {})
        if any(abs(params.get(k, np.nan) - v) > 1e-9 for k, v in filters.items()):
            continue

        table = ak.from_parquet(path)
        if metric not in ak.fields(table):
            raise SystemExit(f"{path} has no column '{metric}'")
        radius = np.asarray(table["radius"], dtype=float)
        values = np.asarray(table[metric], dtype=float)
        is_matched = np.asarray(table["matched"], dtype=float)

        label = meta.get("label", Path(path).stem)
        # Several chunks can share a label; merge rather than overwrite.
        row = quality.setdefault(label, {})
        mrow = matched.setdefault(label, {})
        for r in np.unique(radius):
            sel = radius == r
            if sel.sum() < MIN_JETS_PER_CELL:
                continue
            prev_n, prev_v = row.get(float(r), (0, 0.0))
            n = int(sel.sum())
            row[float(r)] = (prev_n + n,
                             prev_v + float(np.nansum(values[sel])))
            pn, pv = mrow.get(float(r), (0, 0.0))
            mrow[float(r)] = (pn + n, pv + float(np.nansum(is_matched[sel])))
        points[label] = params

    quality = {m: {r: v / n for r, (n, v) in row.items() if n}
               for m, row in quality.items() if row}
    matched = {m: {r: v / n for r, (n, v) in row.items() if n}
               for m, row in matched.items() if row}
    points = {m: points[m] for m in quality}
    if not quality:
        raise SystemExit("no model points survived the filters / stats cut")
    return points, quality, matched


def grid(points, quality, xkey, ykey, reduce_fn):
    """Bin model points onto the (x, y) axes, averaging any duplicates."""
    xs = sorted({points[m][xkey] for m in quality if xkey in points[m]})
    ys = sorted({points[m][ykey] for m in quality if ykey in points[m]})
    if not xs or not ys:
        raise SystemExit(f"model metadata has no '{xkey}'/'{ykey}' params")
    z = np.full((len(ys), len(xs)), np.nan)
    acc = {}
    for m, row in quality.items():
        p = points[m]
        if xkey not in p or ykey not in p:
            continue
        i, j = ys.index(p[ykey]), xs.index(p[xkey])
        acc.setdefault((i, j), []).append(reduce_fn(row))
    for (i, j), vals in acc.items():
        vals = [v for v in vals if v is not None and np.isfinite(v)]
        if vals:
            z[i, j] = float(np.mean(vals))
    return np.array(xs), np.array(ys), z


def heatmap(ax, xs, ys, z, xlabel, ylabel, title, cbar_label, cmap, vmin=None,
            vmax=None, fmt="{:.2f}"):
    mesh = ax.imshow(z, origin="lower", aspect="auto", cmap=cmap,
                     vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xticks(range(len(xs)), [f"{v:g}" for v in xs])
    ax.set_yticks(range(len(ys)), [f"{v:g}" for v in ys])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            if np.isfinite(z[i, j]):
                lo, hi = np.nanmin(z), np.nanmax(z)
                frac = 0.5 if hi == lo else (z[i, j] - lo) / (hi - lo)
                ax.text(j, i, fmt.format(z[i, j]), ha="center", va="center",
                        fontsize=7,
                        color="white" if frac > 0.6 else "black")
    cbar = ax.figure.colorbar(mesh, ax=ax)
    cbar.set_label(cbar_label, fontsize=9)
    return mesh


def curve_plot(outpath, points, rows, choice, xkey, metric, ylabel=None):
    """One curve per model point, coloured by ``points[m][xkey]``.

    Factored out of the original ``quality_vs_radius.png`` step so any column
    ``radius_metrics.scalars()`` emits -- not only ``quality`` -- can be
    plotted the same way, against the same R* line for reference.
    """
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    cvals = [points[m].get(xkey, np.nan) for m in rows]
    lo, hi = np.nanmin(cvals), np.nanmax(cvals)
    cmap = plt.get_cmap("coolwarm")
    for m, row in sorted(rows.items()):
        r = sorted(row)
        c = points[m].get(xkey, np.nan)
        shade = 0.5 if not np.isfinite(c) or hi == lo else (c - lo) / (hi - lo)
        ax.plot(r, [row[x] for x in r], marker="o", ms=2.5, lw=1.0,
                color=cmap(shade), alpha=0.75)
    ax.axvline(choice.radius, color="k", ls="--", lw=1.2,
               label=f"global R* = {choice.radius:.2f}")
    if choice.plateau:
        ax.axvspan(min(choice.plateau), max(choice.plateau), color="k",
                   alpha=0.07, label="plateau (within tolerance)")
    ax.set_xlabel("jet radius R")
    ax.set_ylabel(ylabel or metric)
    ax.set_title(f"{metric} vs radius, one curve per model point "
                 f"(colour = {xkey})", fontsize=11)
    ax.grid(alpha=0.2)
    ax.legend(fontsize=9)
    fig.colorbar(plt.cm.ScalarMappable(
        norm=plt.Normalize(lo, hi), cmap=cmap), ax=ax, label=xkey)
    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    points, quality, matched = load(args.metrics, args.fixed, args.metric)
    fixed_note = ", ".join(args.fixed) if args.fixed else "none (projected)"
    print(f"[INFO] {len(quality)} model points, fixed: {fixed_note}", flush=True)

    choice = choose_radius(quality, quantile=args.quantile,
                           tolerance=args.tolerance)
    print(choice.report(), flush=True)

    if args.curves_only:
        curve_plot(args.outdir / f"{args.metric}_vs_radius.png", points,
                   quality, choice, args.x, args.metric)
        for extra in (c.strip() for c in args.extra_curves.split(",")):
            if not extra:
                continue
            _, extra_rows, _ = load(args.metrics, args.fixed, extra)
            curve_plot(args.outdir / f"{extra}_vs_radius.png", points,
                       extra_rows, choice, args.x, extra)
        summary = {
            "metric": args.metric, "extra_curves": args.extra_curves,
            "fixed": args.fixed, "n_model_points": len(quality),
            "global_radius": choice.radius, "quantile": args.quantile,
            "worst_regret": choice.worst_regret,
            "plateau": list(choice.plateau),
            "bootstrap_16_84": list(choice.radius_uncertainty),
            "outlier_models": list(choice.outliers),
        }
        (args.outdir / "decision_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(f"[DONE] curves-only: {args.metric} plus "
              f"[{args.extra_curves}] in {args.outdir}", flush=True)
        return 0

    # --- 1. the money plot: which radius wins where -----------------------
    xs, ys, z_best = grid(points, quality, args.x, args.y,
                          lambda row: max(row, key=row.get))
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    heatmap(ax, xs, ys, z_best, args.x, args.y,
            f"Optimal radius per model point  (fixed: {fixed_note})",
            "argmax$_R$  " + args.metric, "viridis")
    fig.tight_layout()
    fig.savefig(args.outdir / "optimal_radius_map.png", dpi=160)
    plt.close(fig)

    # An argmax sitting on the first or last radius means the true optimum may
    # lie outside the scanned grid, so the map is censored there and the
    # scaling fit will be biased toward the boundary. Widen the grid rather
    # than quoting the edge value.
    all_radii = sorted({r for row in quality.values() for r in row})
    r_lo, r_hi = all_radii[0], all_radii[-1]
    at_edge = [m for m, row in quality.items()
               if max(row, key=row.get) in (r_lo, r_hi)]
    if at_edge:
        print(f"[WARN] {len(at_edge)}/{len(quality)} model points peak at the "
              f"edge of the radius grid ({r_lo:g} or {r_hi:g}). Their true "
              f"optimum may lie outside it -- extend --radii in scan.py before "
              f"trusting the map or the scaling fit there.", flush=True)
        for m in at_edge[:8]:
            print(f"         {m}  argmax = {max(quality[m], key=quality[m].get):g}",
                  flush=True)

    # --- 2. regret at the global choice -----------------------------------
    def regret_at(row, R=choice.radius):
        if not row:
            return None
        best = max(row.values())
        here = row.get(R, min(row.values()))
        return best - here

    _, _, z_reg = grid(points, quality, args.x, args.y, regret_at)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    heatmap(ax, xs, ys, z_reg, args.x, args.y,
            f"Regret at global R = {choice.radius:.2f}"
            f"  (q={args.quantile:g}, tol={args.tolerance:g})",
            "best $-$ achieved", "magma_r", vmin=0.0, fmt="{:.3f}")
    fig.tight_layout()
    fig.savefig(args.outdir / "regret_map.png", dpi=160)
    plt.close(fig)

    # --- 3. curves: where the plateau is ----------------------------------
    curve_plot(args.outdir / "quality_vs_radius.png", points, quality,
               choice, args.x, args.metric)
    for extra in (c.strip() for c in args.extra_curves.split(",")):
        if not extra:
            continue
        _, extra_rows, _ = load(args.metrics, args.fixed, extra)
        curve_plot(args.outdir / f"{extra}_vs_radius.png", points,
                   extra_rows, choice, args.x, extra)

    # --- 4. why a corner fails --------------------------------------------
    panels = [
        ("acceptance", "whole dark hadrons captured", "viridis"),
        ("frac_partial_pt", "pT from shredded dark hadrons", "magma_r"),
        ("frac_nodh_pt", "pT with no dark-hadron ancestor", "magma_r"),
    ]
    diag = {}
    for column, _, _ in panels:
        diag[column] = load(args.metrics, args.fixed, column)[1]
    fig, axes = plt.subplots(1, 4, figsize=(21, 4.8))
    for axis, (column, label, cmap_name) in zip(axes, panels):
        _, _, zc = grid(points, diag[column], args.x, args.y,
                        lambda row, R=choice.radius: row.get(R, np.nan))
        heatmap(axis, xs, ys, zc, args.x, args.y,
                f"{label}\nat R = {choice.radius:.2f}", column, cmap_name)
    _, _, zm = grid(points, matched, args.x, args.y,
                    lambda row, R=choice.radius: row.get(R, np.nan))
    heatmap(axes[3], xs, ys, zm, args.x, args.y,
            f"matched fraction\nat R = {choice.radius:.2f}",
            "matched", "cividis", vmin=0.0, vmax=1.0)
    fig.tight_layout()
    fig.savefig(args.outdir / "decomposition_maps.png", dpi=160)
    plt.close(fig)

    # --- 5. per-radius regret panels --------------------------------------
    radii = sorted({r for row in quality.values() for r in row})
    ncol = min(5, len(radii))
    nrow = math.ceil(len(radii) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow),
                             squeeze=False)
    vmax = np.nanmax([regret_at(row, R) or 0.0
                      for row in quality.values() for R in radii])
    for k, R in enumerate(radii):
        axis = axes[k // ncol][k % ncol]
        _, _, zr = grid(points, quality, args.x, args.y,
                        lambda row, RR=R: regret_at(row, RR))
        heatmap(axis, xs, ys, zr, args.x, args.y, f"R = {R:.2f}",
                "regret", "magma_r", vmin=0.0, vmax=vmax, fmt="{:.2f}")
    for k in range(len(radii), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle("Regret per radius (normalised per model point, so cells "
                 "are comparable across panels)", fontsize=12)
    fig.tight_layout()
    fig.savefig(args.outdir / "regret_by_radius.png", dpi=160)
    plt.close(fig)

    summary = {
        "metric": args.metric,
        "extra_curves": args.extra_curves,
        "axes": {"x": args.x, "y": args.y},
        "fixed": args.fixed,
        "n_model_points": len(quality),
        "global_radius": choice.radius,
        "quantile": args.quantile,
        "worst_regret": choice.worst_regret,
        "plateau": list(choice.plateau),
        "bootstrap_16_84": list(choice.radius_uncertainty),
        "outlier_models": list(choice.outliers),
        "models_peaking_at_grid_edge": at_edge,
    }
    (args.outdir / "decision_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"[DONE] figures + decision_summary.json in {args.outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
