"""Figures of merit for choosing a jet radius.

Thin by design. ``core.matched_jet_metrics`` and
``core.dark_hadron_group_metrics`` compute the pieces; this module only
combines them into scores and defines what an unmatched truth jet is worth.
Nothing here recomputes anything ``core`` already provides.

The decomposition
-----------------
Every constituent of a matched jet is classified by the **containment status
of the dark hadron it came from**, not by which truth partition owns that dark
hadron:

    FULL      from a dark hadron captured whole by this jet
    PARTIAL   from a dark hadron only fragments of which are in this jet
    NODH      no dark-hadron ancestor at all (ISR, UE, mediator-side SM)

    pt_jet = pt_full + pt_partial + pt_nodh

Ownership is the wrong axis. A dark hadron belonging to *this* truth jet that
is only half captured is just as damaging as a neighbour's leaking in, while a
neighbour's dark hadron arriving fully contained still carries its own correct
mass. ``core``'s ``cross_truth_dark_hadron_contamination_pt`` splits on
ownership and so mixes those two cases; it is still written out for continuity
but is not what the scores below use.

PARTIAL is the worst of the three cases, and worse than simply missing a dark
hadron. Fully contained, a dark hadron's decay products sum to its mass.
Entirely absent, it costs acceptance but tells no lies. Half in, the jet gains
a random fraction of its momentum with no mass meaning attached -- heavier and
less correct at once.

The two scores
--------------
    A = n_fully_contained_dark_hadrons / n_eligible_dark_hadrons
    P = pt_full / pt_jet = 1 - frac_partial_pt - frac_nodh_pt

A is a **count**, deliberately. A dark hadron is either whole or it is not;
pT-weighting a binary means nothing.

The two line up with the physics without a tuning constant, because a
partially captured dark hadron is penalised **twice**: it fails to count in
A's numerator, and its fragments land in ``frac_partial_pt``, which drags down
P. A dark hadron entirely outside the jet costs only A. Junk with no
dark-hadron ancestor costs only P.

Small R leaves dark hadrons shredded at the jet edge and both terms suffer.
Large R lets whole neighbours in (P barely moves, they are fully contained)
along with UE and ISR (P falls through ``frac_nodh_pt``).

``iou_pt`` is kept as a cross-check. It scores set overlap, which is a
different object from "how many complete dark hadrons survived", and it weights
each constituent by pT linearly while jet mass responds to pT x dR^2 -- so it
under-charges the soft wide junk that a too-large radius admits. Do not use it
as the headline.
"""

from __future__ import annotations

import math
from typing import Any, Mapping

__all__ = [
    "acceptance",
    "purity",
    "quality",
    "purity_dh",
    "quality_dh",
    "iou_pt",
    "iou_pt_dh",
    "orphan_dh_pt_fraction",
    "scalars",
    "unmatched_scalars",
    "SCALAR_NAMES",
]

SCALAR_NAMES = (
    "acceptance",
    "purity",
    "quality",
    "purity_dh",
    "quality_dh",
    "frac_partial_pt",
    "frac_nodh_pt",
    "frac_full_pt",
    "frac_dh_pt",
    "frac_orphan_dh_pt",
    "n_dark_hadrons_full",
    "n_dark_hadrons_eligible",
    "n_dark_hadrons_partial",
    "iou_pt",
    "iou_pt_dh",
    "containment_pt",
    "matched",
)


def _finite(value: Any) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def acceptance(metrics: Mapping[str, Any]) -> float:
    """Fraction of this truth jet's dark hadrons captured whole.

    Denominator is dark hadrons with at least one resolved visible descendant;
    one with none cannot be called contained and is outside the question.
    """
    full = _finite(metrics.get("n_fully_contained_dark_hadrons"))
    eligible = _finite(metrics.get("n_dark_hadrons_with_visible_descendants"))
    if math.isnan(full) or math.isnan(eligible) or eligible <= 0.0:
        return math.nan
    return full / eligible


def purity(metrics: Mapping[str, Any]) -> float:
    """Fraction of jet pT coming from dark hadrons captured whole.

    Caveat: ``1 - frac_partial_pt - frac_nodh_pt`` silently folds in
    ``orphan_dh_pt_fraction`` (a dark hadron with a real ancestor that was
    just never assigned to any truth partition) as if it were clean pT --
    see that function's docstring. Not corrected here; existing plots and
    closure gates consume this value as defined.
    """
    partial = _finite(metrics.get("partially_contained_dark_hadron_pt_fraction"))
    nodh = _finite(metrics.get("non_dark_hadron_constituent_pt_fraction"))
    if math.isnan(partial) or math.isnan(nodh):
        return math.nan
    return max(0.0, 1.0 - partial - nodh)


def quality(metrics: Mapping[str, Any]) -> float:
    """``A * P`` -- the headline scalar fed to ``optimize.choose_radius``.

    Bounded in [0, 1] and zero if either term collapses, so a radius cannot buy
    a good score by capturing whole dark hadrons into a jet full of junk, nor
    by producing a spotless jet that shredded most of the dark hadrons.

    Report the two terms alongside it. The product says *how good*; only the
    split says *why*, and the writeup needs the why.
    """
    a, p = acceptance(metrics), purity(metrics)
    if math.isnan(a) or math.isnan(p):
        return math.nan
    return a * p


def purity_dh(metrics: Mapping[str, Any]) -> float:
    """Purity of the jet's dark-hadron content specifically.

    ``frac_full_pt / frac_dh_pt`` -- of the pT in this jet that traces to ANY
    dark hadron at all (whole, partial, or orphan), what fraction came from
    ones captured whole. Unlike ``purity``, background (``frac_nodh_pt``)
    never enters this at all: a jet that is 90% UE/ISR but whose small
    dark-hadron slice is entirely clean scores 1.0 here, where ``purity``
    would be dragged down by the background it never claimed to measure.

    Also closes the ``purity`` caveat directly: ``frac_orphan_dh_pt`` sits in
    the denominator here, penalised the same way a partial capture is,
    instead of riding along inside the numerator as ``purity`` does.

    Requested and defined by the user, 2026-09-28; added alongside ``purity``
    rather than replacing it, following the rule used for ``contamination_pt``:
    add a new field, do not redefine a validated one.

    NaN if the jet has no dark-hadron-sourced pT at all (``frac_dh_pt`` ~ 0):
    "how clean is the dark-hadron part" has no answer when there isn't one.
    """
    full = _finite(metrics.get("fully_contained_dark_hadron_pt_fraction"))
    dh = 1.0 - _finite(metrics.get("non_dark_hadron_constituent_pt_fraction"))
    if math.isnan(full) or math.isnan(dh) or dh <= 1e-9:
        return math.nan
    return full / dh


def quality_dh(metrics: Mapping[str, Any]) -> float:
    """``A * purity_dh`` -- quality using the dark-hadron-scoped purity.

    Same shape as ``quality``, with ``purity`` swapped for ``purity_dh``. Not
    the headline metric ``optimize.choose_radius`` consumes; a companion for
    comparing the two purity definitions' effect on the radius choice.
    """
    a, p = acceptance(metrics), purity_dh(metrics)
    if math.isnan(a) or math.isnan(p):
        return math.nan
    return a * p


def orphan_dh_pt_fraction(metrics: Mapping[str, Any]) -> float:
    """pT from a dark hadron that exists but was never assigned to any
    truth ``DarkHadronJet`` partition.

    ``core.py`` classifies constituents into three named buckets --
    ``fully_contained_dark_hadron_pt_fraction``,
    ``partially_contained_dark_hadron_pt_fraction``,
    ``non_dark_hadron_constituent_pt_fraction`` -- but the first two only sum
    over dark hadrons owned by *some* truth jet (``target_records +
    contaminating_records`` in ``core.dark_hadron_group_metrics``), while the
    third tests ancestry against *every* ``DarkHadronCandidate`` in the event
    regardless of ownership. A dark hadron that exists but was never assigned
    to a truth partition therefore counts toward neither "full"/"partial" nor
    "nodh" -- it is missing from all three, which is why they do not sum to
    1.0. This is that missing residual, made explicit rather than silently
    absorbed into ``purity`` (see ``purity``'s docstring caveat).

    Measured non-zero on both the DRAGON mmed=2000 mass grid and the local
    mmed=1000 validation sample: mean ~0.014, up to 0.85 on individual jets.
    """
    full = _finite(metrics.get("fully_contained_dark_hadron_pt_fraction"))
    partial = _finite(metrics.get("partially_contained_dark_hadron_pt_fraction"))
    nodh = _finite(metrics.get("non_dark_hadron_constituent_pt_fraction"))
    return 1.0 - full - partial - nodh


def iou_pt(metrics: Mapping[str, Any]) -> float:
    """pT-weighted Jaccard index. Cross-check only -- see the module docstring."""
    a = _finite(metrics.get("containment_pt"))
    p = purity(metrics)
    if math.isnan(a) or math.isnan(p):
        return math.nan
    if a <= 0.0 or p <= 0.0:
        return 0.0
    return 1.0 / (1.0 / a + 1.0 / p - 1.0)


def iou_pt_dh(metrics: Mapping[str, Any]) -> float:
    """Same construction as ``iou_pt``, with ``purity`` swapped for
    ``purity_dh``: ``1 / (1/containment_pt + 1/purity_dh - 1)``.

    Cross-check only, same as ``iou_pt`` -- and inherits the same caveat that
    ``containment_pt`` (recall over the target jet's own visible pT) and
    ``purity_dh`` (precision over the jet's dark-hadron-sourced pT, own +
    contaminating) don't share one intersection set, so this is an IoU-shaped
    blend, not a literal Jaccard index.
    """
    a = _finite(metrics.get("containment_pt"))
    p = purity_dh(metrics)
    if math.isnan(a) or math.isnan(p):
        return math.nan
    if a <= 0.0 or p <= 0.0:
        return 0.0
    return 1.0 / (1.0 / a + 1.0 / p - 1.0)


def scalars(metrics: Mapping[str, Any]) -> dict[str, float]:
    """Derived scores for one matched truth jet, from core's metric dict."""
    return {
        "acceptance": acceptance(metrics),
        "purity": purity(metrics),
        "quality": quality(metrics),
        "purity_dh": purity_dh(metrics),
        "quality_dh": quality_dh(metrics),
        "frac_partial_pt": _finite(
            metrics.get("partially_contained_dark_hadron_pt_fraction")),
        "frac_nodh_pt": _finite(
            metrics.get("non_dark_hadron_constituent_pt_fraction")),
        "frac_full_pt": _finite(
            metrics.get("fully_contained_dark_hadron_pt_fraction")),
        "frac_dh_pt": 1.0 - _finite(
            metrics.get("non_dark_hadron_constituent_pt_fraction")),
        "frac_orphan_dh_pt": orphan_dh_pt_fraction(metrics),
        "n_dark_hadrons_full": _finite(
            metrics.get("n_fully_contained_dark_hadrons")),
        "n_dark_hadrons_eligible": _finite(
            metrics.get("n_dark_hadrons_with_visible_descendants")),
        "n_dark_hadrons_partial": _finite(
            metrics.get("n_partially_contained_dark_hadrons_in_jet")),
        "iou_pt": iou_pt(metrics),
        "iou_pt_dh": iou_pt_dh(metrics),
        "containment_pt": _finite(metrics.get("containment_pt")),
        "matched": 1.0,
    }


def unmatched_scalars(n_eligible: float = math.nan) -> dict[str, float]:
    """Scores for a truth jet with no match at this radius.

    Zero, not NaN, and not omitted. A truth jet that fails to match still
    exists; excluding it rewards radii that fail to reconstruct the object at
    all, which is the easiest way to argue yourself into too small a radius.

    Purity is NaN rather than 0 -- there is no jet, so "what fraction of the
    jet is clean" has no answer. ``quality`` is 0 because no dark hadron was
    captured whole, which is the statement that matters.
    """
    return {
        "acceptance": 0.0,
        "purity": math.nan,
        "quality": 0.0,
        "purity_dh": math.nan,
        "quality_dh": 0.0,
        "frac_partial_pt": math.nan,
        "frac_nodh_pt": math.nan,
        "frac_full_pt": math.nan,
        "frac_dh_pt": math.nan,
        "frac_orphan_dh_pt": math.nan,
        "n_dark_hadrons_full": 0.0,
        "n_dark_hadrons_eligible": _finite(n_eligible),
        "n_dark_hadrons_partial": 0.0,
        "iou_pt": 0.0,
        "iou_pt_dh": 0.0,
        "containment_pt": 0.0,
        "matched": 0.0,
    }
