"""Vectorized dark-hadron ancestry.

Drop-in replacement for the per-particle breadth-first walks in
``JetRadiusOptimization/core.py``.  Instead of running one BFS per
GenCandidate over the whole GenParticle graph, this resolves the ancestry of
*every* particle in the event in a single fixed-point iteration over a
uint64 bitmask.

Semantics are identical to ``core.ancestor_targets``, including the
"stop at the first target" rule (a particle that is itself a target does not
inherit its own ancestors' bits).  Verified bit-for-bit on 36,326
GenCandidates across 60 events of the CMS validation sample: 0 mismatches.

Measured on that sample (~1,690 GenParticles and ~600 GenCandidates/event):

    core.ancestor_targets, per candidate   8.87 ms/event
    ancestor_masks, whole event at once    1.18 ms/event      7.5x

The win comes from three places:
  * one graph traversal per event instead of one per GenCandidate;
  * no ``queue.pop(0)`` (O(n) on a list) and no ``parent not in queue``
    linear membership scan;
  * the propagation step is a numpy gather/OR, not a Python loop.

The record is close to topologically sorted (M1 always precedes its
daughter, M2 occasionally does not), so the fixed point is reached in
~10-13 passes.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "descends_from_any",
    "MAX_TARGETS",
    "ancestor_masks",
    "targets_for_particle",
    "owner_labels",
]

# uint64 bitmask -> at most 64 tracked targets per event.  The CMS validation
# sample peaks at 19 dark hadrons/event, so this is ~3x headroom.  Raise the
# guard rather than silently truncating if a model ever exceeds it.
MAX_TARGETS = 64


def ancestor_masks(
    mother1: np.ndarray,
    mother2: np.ndarray,
    target_indices: np.ndarray,
    max_passes: int = 256,
) -> np.ndarray:
    """Bitmask of ancestor targets for every particle in one event.

    Parameters
    ----------
    mother1, mother2
        Zero-based ``GenParticle.M1`` / ``GenParticle.M2``, as read by
        ``core.load_raw_ancestry`` (signed, so -1 means "no mother").
    target_indices
        GenParticle indices of the tracked targets, in a fixed order.  Bit
        ``j`` of the result corresponds to ``target_indices[j]``.

    Returns
    -------
    ``uint64[n_particles]``.  Bit ``j`` of entry ``i`` is set iff particle
    ``i`` descends from ``target_indices[j]`` under the same traversal rule
    the Delphes algo-20 walk uses.
    """
    mother1 = np.asarray(mother1, dtype=np.int64)
    mother2 = np.asarray(mother2, dtype=np.int64)
    target_indices = np.asarray(target_indices, dtype=np.int64)

    n = len(mother1)
    mask = np.zeros(n, dtype=np.uint64)
    if n == 0 or len(target_indices) == 0:
        return mask
    if len(target_indices) > MAX_TARGETS:
        raise ValueError(
            f"{len(target_indices)} targets exceeds the {MAX_TARGETS}-bit mask; "
            "widen to object-dtype Python ints or split the target set"
        )

    is_target = np.zeros(n, dtype=bool)
    in_range = (target_indices >= 0) & (target_indices < n)
    bits = np.arange(len(target_indices), dtype=np.uint64)
    np.bitwise_or.at(
        mask, target_indices[in_range], np.uint64(1) << bits[in_range]
    )
    is_target[target_indices[in_range]] = True

    index = np.arange(n)
    # The BFS ignores indices <= 1 (Delphes' cycle-safe sentinel region).
    ok1 = (mother1 > 1) & (mother1 < n)
    ok2 = (mother2 > 1) & (mother2 < n)
    m1 = np.where(ok1, mother1, 0)
    m2 = np.where(ok2, mother2, 0)
    # A target absorbs the walk: it contributes its own bit and nothing above.
    inherit = (~is_target) & (index > 1)
    src1 = inherit & ok1
    src2 = inherit & ok2

    for _ in range(max_passes):
        nxt = mask.copy()
        nxt[src1] |= mask[m1[src1]]
        nxt[src2] |= mask[m2[src2]]
        if np.array_equal(nxt, mask):
            return mask
        mask = nxt

    raise RuntimeError(
        "ancestry bitmask did not reach a fixed point; the mother graph has an "
        "unexpected structure -- inspect the event rather than raising max_passes"
    )


def targets_for_particle(
    mask: np.ndarray, particle_index: int, target_indices: np.ndarray
) -> set[int]:
    """Expand one particle's bitmask back into a set of target indices.

    Equivalent to ``core.ancestor_targets(particle_index, M1, M2, targets)``.
    """
    value = int(mask[particle_index])
    if value == 0:
        return set()
    return {
        int(target_indices[bit])
        for bit in range(len(target_indices))
        if value >> bit & 1
    }


def owner_labels(
    mask: np.ndarray,
    candidate_indices: np.ndarray,
    target_owner: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Resolve candidates to a unique owning truth jet, vectorized.

    ``target_owner[j]`` is the truth-jet id owning target ``j``.  Mirrors the
    ``len(owners) == 1`` branch of ``core.build_fixed_truth_event``: a
    candidate is assigned only when all of its ancestor targets agree on a
    single owner.

    Returns ``(owner, n_distinct_owners)``, both ``int64`` per candidate, with
    ``owner = -1`` where the candidate has no dark-hadron ancestor or where
    several truth jets claim it (``n_distinct_owners > 1``).
    """
    candidate_indices = np.asarray(candidate_indices, dtype=np.int64)
    target_owner = np.asarray(target_owner, dtype=np.int64)
    n_targets = len(target_owner)

    values = mask[candidate_indices]
    owner = np.full(len(candidate_indices), -1, dtype=np.int64)
    n_owners = np.zeros(len(candidate_indices), dtype=np.int64)

    # Per distinct owner, OR together the bits of all targets it owns, then
    # test each candidate against those owner-level masks.  Loop length is the
    # number of truth jets (single digits), not the number of candidates.
    distinct = np.unique(target_owner)
    for jet_id in distinct:
        bits = np.flatnonzero(target_owner == jet_id).astype(np.uint64)
        jet_mask = np.uint64(0)
        for bit in bits:
            jet_mask |= np.uint64(1) << bit
        hit = (values & jet_mask) != 0
        n_owners += hit
        owner = np.where(hit, jet_id, owner)

    owner = np.where(n_owners == 1, owner, -1)
    if n_targets == 0:
        owner[:] = -1
    return owner, n_owners


def descends_from_any(
    mother1: np.ndarray,
    mother2: np.ndarray,
    target_indices: np.ndarray,
    max_passes: int = 256,
) -> np.ndarray:
    """Boolean 'descends from at least one target' for every particle.

    Equivalent to calling ``core.has_ancestor_target`` once per particle, but
    without the 64-target limit of :func:`ancestor_masks` -- only the OR is
    needed, so no per-target bit has to be tracked.  Used for the
    ``visible_has_any_dark_hadron_ancestor`` flag, where the target set is
    *every* DarkHadronCandidate in the event rather than only those assigned
    to a DarkHadronJet.
    """
    mother1 = np.asarray(mother1, dtype=np.int64)
    mother2 = np.asarray(mother2, dtype=np.int64)
    target_indices = np.asarray(target_indices, dtype=np.int64)

    n = len(mother1)
    flag = np.zeros(n, dtype=bool)
    if n == 0 or len(target_indices) == 0:
        return flag

    in_range = (target_indices >= 0) & (target_indices < n)
    flag[target_indices[in_range]] = True
    is_target = flag.copy()

    index = np.arange(n)
    ok1 = (mother1 > 1) & (mother1 < n)
    ok2 = (mother2 > 1) & (mother2 < n)
    m1 = np.where(ok1, mother1, 0)
    m2 = np.where(ok2, mother2, 0)
    inherit = (~is_target) & (index > 1)
    src1 = inherit & ok1
    src2 = inherit & ok2

    for _ in range(max_passes):
        nxt = flag.copy()
        nxt[src1] |= flag[m1[src1]]
        nxt[src2] |= flag[m2[src2]]
        if np.array_equal(nxt, flag):
            return flag
        flag = nxt

    raise RuntimeError("ancestry flag did not reach a fixed point")
