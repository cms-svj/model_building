"""Dark-hadron ancestry for Delphes events.

Labels generator-level particles with the dark hadrons they descend from, by
following the GenParticle mother links (M1, M2) upward. For each particle the
result is a bitmask: bit j is set if the particle descends from target j.

The walk uses the same rules as the DarkHadronVisibleMatch mode
(JetAlgorithm 20) of FastJetFinder in the cms-svj Delphes fork, which fills
DarkHadronVisibleJet:
  * GenParticle indices 0 and 1 (the incoming beams) end the walk;
  * a target ends the walk: a particle descended from target j is not also
    credited to the targets above j.
One difference: Delphes stops at the first dark hadron it reaches, while this
module records every dark hadron reachable through M1 or M2. A particle with
more than one dark-hadron ancestor therefore has more than one bit set.

Usage, with events from common.load_events:
    masks = dark_hadron_ancestry(events)   # one uint64 per GenCandidate
    from_dh_0 = (masks & 1) != 0            # descends from DarkHadronCandidate 0

Check the labels against the stored DarkHadronVisibleJet with
    python3 darkancestry.py -i events.root
"""

from collections import Counter

import awkward as ak
import numpy as np

from common import get_constituents_kernel, load_events

# bits in a uint64 mask; the CMS benchmark model has at most 22 dark hadrons per event
MAX_TARGETS = 64


def genparticle_index(uids, gen_uids):
    """Index into GenParticle of each particle with unique ID uids (per event)."""
    return ak.Array([
        get_constituents_kernel(
            ak.to_numpy(u).astype(np.int32), ak.to_numpy(g).astype(np.uint32)
        )
        for u, g in zip(uids, gen_uids)
    ])


def ancestor_masks(mother1, mother2, targets):
    """Bitmask of target ancestors for every GenParticle.

    mother1, mother2: GenParticle.M1 and M2 (per event, -1 = no mother)
    targets: GenParticle indices of the targets (per event); bit j of the
             result corresponds to targets[j]
    """
    n_targets = ak.num(targets)
    if len(n_targets) and ak.max(n_targets) > MAX_TARGETS:
        raise ValueError(f"more than {MAX_TARGETS} targets in one event")

    # flatten all events into one array; event boundaries via offsets
    counts = ak.to_numpy(ak.num(mother1))
    offsets = np.concatenate([[0], np.cumsum(counts)[:-1]])
    n_in_event = np.repeat(counts, counts)
    first = np.repeat(offsets, counts)
    local = np.arange(len(n_in_event)) - first

    mask = np.zeros(len(n_in_event), dtype=np.uint64)
    target_global = ak.to_numpy(ak.flatten(targets + offsets)).astype(np.int64)
    target_bit = ak.to_numpy(ak.flatten(ak.local_index(targets))).astype(np.uint64)
    np.bitwise_or.at(mask, target_global, np.uint64(1) << target_bit)
    is_target = np.zeros(len(mask), dtype=bool)
    is_target[target_global] = True

    # each particle inherits the bits of its mothers, except targets and beams
    inherits = ~is_target & (local > 1)
    links = []
    for mother in (mother1, mother2):
        m = ak.to_numpy(ak.flatten(mother)).astype(np.int64)
        ok = inherits & (m > 1) & (m < n_in_event)
        links.append((np.flatnonzero(ok), (m + first)[ok]))

    # iterate until no mask changes; converges in O(depth of the decay chain)
    while True:
        new = mask.copy()
        for daughter, mother in links:
            new[daughter] |= mask[mother]
        if np.array_equal(new, mask):
            return ak.unflatten(mask, counts)
        mask = new


def dark_hadron_ancestry(events, dark_hadrons=None):
    """Bitmask per GenCandidate: bit j set if it descends from dark hadron j.

    dark_hadrons: indices into DarkHadronCandidate (per event) to use as
                  targets, in bit order; default is all of them
    """
    gen_uids = events.GenParticle.fUniqueID
    dh_index = genparticle_index(events.DarkHadronCandidate.fUniqueID, gen_uids)
    if dark_hadrons is not None:
        dh_index = dh_index[dark_hadrons]
    # item access on purpose: attribute access .M2 is vector's mass squared
    masks = ancestor_masks(events.GenParticle["M1"], events.GenParticle["M2"], dh_index)
    return masks[genparticle_index(events.GenCandidate.fUniqueID, gen_uids)]


def check_visible_jets(events):
    """Compare ancestry with the stored DarkHadronVisibleJet constituents.

    For each DarkHadronJet, the GenCandidates descended from its dark hadrons
    should be exactly the constituents of one DarkHadronVisibleJet. Visible jets
    are pT-ordered separately, so jets are compared as sets within each event.
    """
    # targets: the dark hadrons in DarkHadronJets, jet by jet (as Delphes uses)
    dh_uids = events.DarkHadronCandidate.fUniqueID
    dh_per_jet = ak.Array([
        [get_constituents_kernel(ak.to_numpy(r), ak.to_numpy(u)) for r in refs]
        for refs, u in zip(events.DarkHadronJet.Constituents.refs, dh_uids)
    ])
    # jet number of each target bit
    jet_of_bit = ak.flatten(
        ak.broadcast_arrays(ak.local_index(dh_per_jet, axis=1), dh_per_jet)[0], axis=2
    )
    masks = dark_hadron_ancestry(events, ak.flatten(dh_per_jet, axis=2))

    stats = {"events": len(events), "jets": 0, "mismatched_jets": 0, "ambiguous_candidates": 0}
    for mask, jet, refs, uids in zip(
        masks, jet_of_bit, events.DarkHadronVisibleJet.Constituents.refs,
        events.GenCandidate.fUniqueID,
    ):
        mask, jet, uids = ak.to_numpy(mask), ak.to_numpy(jet), ak.to_numpy(uids)
        n_jets = int(jet.max()) + 1 if len(jet) else 0
        # descends[j, k]: GenCandidate k descends from a dark hadron of jet j
        descends = np.zeros((n_jets, len(mask)), dtype=bool)
        for bit, j in enumerate(jet):
            descends[j] |= (mask >> np.uint64(bit)) & np.uint64(1) == 1
        stats["ambiguous_candidates"] += int((descends.sum(axis=0) > 1).sum())

        ours = Counter(tuple(np.flatnonzero(d)) for d in descends)
        stored = Counter(
            tuple(sorted(get_constituents_kernel(ak.to_numpy(r), uids))) for r in refs
        )
        matched = sum((ours & stored).values())
        stats["jets"] += n_jets
        stats["mismatched_jets"] += max(n_jets, len(refs)) - matched
    stats["passed"] = stats["mismatched_jets"] == 0 and stats["ambiguous_candidates"] == 0
    return stats


if __name__ == "__main__":
    from magiconfig import ArgumentParser, ArgumentDefaultsRawHelpFormatter

    parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsRawHelpFormatter)
    parser.add_argument("-i", "--input", required=True, help="Delphes events.root")
    parser.add_argument("-n", "--max-events", type=int, default=None, help="number of events to check")
    args = parser.parse_args()

    events = load_events(args.input)
    if args.max_events is not None:
        events = events[: args.max_events]
    print(check_visible_jets(events))
