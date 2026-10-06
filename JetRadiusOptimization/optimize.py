"""Pick one radius that works across the whole model grid.

The obvious move -- average the figure of merit over all model points and take
the argmax -- is wrong for this question.  Averaging lets a large score on
easy models (heavy mediator, low r_inv, few dark hadrons) mask a bad score on
hard ones, and the task force needs a radius that is *defensible everywhere*,
not one that is excellent on average.

This module instead works in **regret**: for each model point m,

    regret_m(R) = Q_m(R_m*) - Q_m(R)

is how much figure of merit you give up at m by using the global radius R
instead of m's own best radius.  Then

    R* = argmin_R  quantile_q [ regret_m(R) ]     over models m

with q = 1.0 recovering strict minimax and q ~ 0.9 giving a version that is
not hostage to a single pathological corner of the grid.

Two things get reported alongside R*:

  * the **plateau** -- every radius whose worst-case regret is within the
    bootstrap uncertainty of the best.  The honest deliverable is almost always
    "0.8 to 1.0 are indistinguishable, take 0.8" rather than a single number
    quoted to two decimals.
  * the **outliers** -- model points whose regret at R* exceeds a tolerance.
    If those cluster in a corner of parameter space (they usually do: low
    mediator mass, high r_inv), that is the physics result, and it is a much
    more useful answer than the number itself.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np

__all__ = ["RadiusChoice", "choose_radius", "fit_scaling", "regret_table"]


@dataclass
class RadiusChoice:
    radius: float
    quantile: float
    worst_regret: float
    plateau: tuple[float, ...]
    per_model_regret: dict[str, float]
    outliers: tuple[str, ...]
    bootstrap_radii: tuple[float, ...] = field(default=())

    @property
    def radius_uncertainty(self) -> tuple[float, float]:
        """16-84% interval on R* from the model-level bootstrap."""
        if not self.bootstrap_radii:
            return (np.nan, np.nan)
        arr = np.asarray(self.bootstrap_radii, dtype=float)
        return (float(np.percentile(arr, 16)), float(np.percentile(arr, 84)))

    def report(self) -> str:
        lo, hi = self.radius_uncertainty
        lines = [
            f"R* = {self.radius:.2f}",
            f"  selection rule      : minimise the {100*self.quantile:.0f}th-percentile "
            f"regret across {len(self.per_model_regret)} model points",
            f"  worst-case regret   : {self.worst_regret:.4f}",
            f"  plateau             : {', '.join(f'{r:.2f}' for r in self.plateau)}",
        ]
        if np.isfinite(lo):
            lines.append(f"  bootstrap 16-84%    : [{lo:.2f}, {hi:.2f}]")
        if self.outliers:
            lines.append(f"  models above tol    : {len(self.outliers)}")
            for name in self.outliers[:10]:
                lines.append(f"      {name}  regret={self.per_model_regret[name]:.4f}")
            if len(self.outliers) > 10:
                lines.append(f"      ... and {len(self.outliers)-10} more")
        else:
            lines.append("  models above tol    : none")
        return "\n".join(lines)


def regret_table(
    quality: Mapping[str, Mapping[float, float]],
    radii: Sequence[float] | None = None,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Turn {model: {radius: Q}} into a dense regret matrix.

    Returns ``(model_names, radii, regret)`` with ``regret`` shaped
    ``(n_models, n_radii)``.  NaN entries (a radius never evaluated for that
    model) propagate, so callers should check before reducing.
    """
    models = sorted(quality)
    if radii is None:
        seen: set[float] = set()
        for row in quality.values():
            seen.update(float(r) for r in row)
        radii_arr = np.array(sorted(seen), dtype=float)
    else:
        radii_arr = np.asarray(sorted(float(r) for r in radii), dtype=float)

    q = np.full((len(models), len(radii_arr)), np.nan)
    for i, model in enumerate(models):
        row = quality[model]
        for j, radius in enumerate(radii_arr):
            if radius in row:
                q[i, j] = float(row[radius])
            else:
                key = next((k for k in row if np.isclose(float(k), radius)), None)
                if key is not None:
                    q[i, j] = float(row[key])

    best = np.nanmax(q, axis=1, keepdims=True)
    return models, radii_arr, best - q


def choose_radius(
    quality: Mapping[str, Mapping[float, float]],
    radii: Sequence[float] | None = None,
    quantile: float = 0.90,
    tolerance: float = 0.01,
    n_bootstrap: int = 500,
    seed: int = 12345,
) -> RadiusChoice:
    """Select the single radius minimising high-quantile cross-model regret.

    Parameters
    ----------
    quality
        ``{model_label: {radius: figure_of_merit}}``.  Use the *mean over jets*
        of ``metrics.iou_pt`` for the headline result, with unmatched truth
        jets already entered as zero.
    quantile
        1.0 for strict minimax; 0.90 to let the worst ~10% of model points off
        the hook.  Report both -- if they disagree, the grid has a corner that
        the simplified model genuinely cannot cover with one radius, which is
        itself a finding worth writing down.
    tolerance
        Absolute figure-of-merit loss still considered acceptable, used only to
        flag outlier models.
    n_bootstrap
        Resamples *over model points* (not events) to get an uncertainty band
        on R*.  Event-level uncertainty belongs upstream, in Q_m itself.
    """
    models, radii_arr, regret = regret_table(quality, radii)
    if not models:
        raise ValueError("no model points supplied")

    def pick(rows: np.ndarray) -> tuple[int, np.ndarray]:
        sub = regret[rows]
        with np.errstate(invalid="ignore"):
            score = np.nanquantile(sub, quantile, axis=0)
        score = np.where(np.all(np.isnan(sub), axis=0), np.inf, score)
        return int(np.nanargmin(score)), score

    best_j, score = pick(np.arange(len(models)))
    best_radius = float(radii_arr[best_j])
    worst = float(score[best_j])

    plateau = tuple(
        float(r) for r, s in zip(radii_arr, score) if np.isfinite(s) and s <= worst + tolerance
    )

    per_model = {m: float(regret[i, best_j]) for i, m in enumerate(models)}
    outliers = tuple(
        m for m in models if np.isfinite(per_model[m]) and per_model[m] > tolerance
    )

    rng = np.random.default_rng(seed)
    boot: list[float] = []
    for _ in range(n_bootstrap):
        rows = rng.integers(0, len(models), size=len(models))
        try:
            j, _ = pick(rows)
            boot.append(float(radii_arr[j]))
        except ValueError:
            continue

    return RadiusChoice(
        radius=best_radius,
        quantile=quantile,
        worst_regret=worst,
        plateau=plateau,
        per_model_regret=per_model,
        outliers=outliers,
        bootstrap_radii=tuple(boot),
    )


def fit_scaling(
    optimal_radius: Mapping[str, float],
    predictor: Mapping[str, float],
    log: bool = True,
) -> dict[str, Any]:
    """Fit R*_m against a candidate dimensionless predictor.

    The most useful outcome of the whole study is not a number but a *rule*.
    If R*_m collapses onto one curve against, say, 2*m_med/pT_jet, then the
    recommendation to experiments becomes a formula that stays valid outside
    the scanned grid, which is far stronger than a fixed 0.8.

    Candidate predictors worth trying, in rough order of expected explanatory
    power:

        2 * m_med / pT_jet_visible     opening angle of a two-prong decay
        m_dark / pT_jet_visible        angular scale of the dark hadron itself
        Lambda / m_med                 shower hardness
        m_pi / Lambda                  2-body vs 3-body regime
        r_inv                          shifts pT_visible, so partly degenerate
                                       with the first two -- check the residual
                                       after regressing on those

    Returns slope, intercept, and R^2 of a straight-line fit in either linear
    or log-log space.  A low R^2 is informative: it says the radius does *not*
    follow that scaling and a fixed value is the right recommendation.
    """
    shared = sorted(set(optimal_radius) & set(predictor))
    if len(shared) < 3:
        raise ValueError("need at least 3 shared model points to fit")

    x = np.array([float(predictor[m]) for m in shared])
    y = np.array([float(optimal_radius[m]) for m in shared])
    good = np.isfinite(x) & np.isfinite(y)
    if log:
        good &= (x > 0) & (y > 0)
    x, y = x[good], y[good]
    if len(x) < 3:
        raise ValueError("too few finite points after filtering")

    xf, yf = (np.log(x), np.log(y)) if log else (x, y)
    slope, intercept = np.polyfit(xf, yf, 1)
    pred = slope * xf + intercept
    ss_res = float(np.sum((yf - pred) ** 2))
    ss_tot = float(np.sum((yf - yf.mean()) ** 2))

    return {
        "n_points": int(len(x)),
        "space": "log-log" if log else "linear",
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan,
        "residual_rms": float(np.sqrt(ss_res / len(x))),
        "models": shared,
    }
