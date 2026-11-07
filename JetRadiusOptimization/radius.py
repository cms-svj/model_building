"""Containment and contamination of dark-hadron jets as a function of jet radius,
at two stages of the shower.

Truth jets are anti-kT jets of the dark hadrons (DarkHadronCandidate) at radius
R_DH; visible jets are anti-kT jets of the visible particles (GenCandidate) at
radius R. Both use pT > 15 GeV, so at R_DH = R = 0.8 they are the Delphes
DarkHadronJet and GenFatJet. For each (R_DH, R), each truth jet is matched
one-to-one to the visible jet that shares the most visible pT with it, and the
matched jet's pT is split by where it came from, using darkancestry (complete
and partial as defined below):

    f_full     from complete dark hadrons
    f_partial  from partial dark hadrons (some, but not enough, visible pT in the jet)
    f_nodh     from particles with no dark-hadron ancestor (ISR, UE, ...)

These add up to 1 (a particle descended from two dark hadrons would count for
both; there are none in the CMS benchmark model).

All scores are fractions of pT. For a truth jet T matched to a visible jet J,
with pT(S) the scalar pT sum over a set of particles S, V_d the visible
descendants of dark hadron d and V_T the union of V_d over the dark hadrons of T:

    d is complete in J  if  pT(V_d & J) >= min_containment * pT(V_d)
                        (default 0.95; 1 requires every visible descendant)

    acceptance = sum over complete d in T  of pT(V_d & J) / pT(V_T)
                 fraction of the truth jet's visible pT that J holds through
                 complete dark hadrons
    purity     = sum over complete d in J  of pT(V_d & J) / pT(J)  (= f_full)
                 fraction of J's pT from complete dark hadrons, of any truth jet
    quality    = acceptance * purity
    iou        = pT(V_T & J) / pT(V_T | J)
                 overlap of J with the truth jet's visible particles; does not
                 depend on min_containment

n_dh and n_full count the dark hadrons of T with visible descendants and the
complete ones. Truth jets with n_dh = 0 (only invisible decays) are skipped. An
unmatched truth jet has acceptance = quality = iou = 0 and no purity.

Jet shapes, with the Histogram.py definitions (DeltaR in eta-phi):
    truth_r90/95/99  radius around the dark-hadron jet axis containing 90/95/99%
                     of pT(V_T), i.e. the inverse of the integrated jet shape
    psi_truth        Psi_T(R) = pT(V_T within R of that axis) / pT(V_T)
    jet_r90/95/99, jet_girth, jet_ptd, jet_major, jet_minor, jet_nconst
                     shapes of the matched visible jet J over its own constituents

Check the clustering against the stored GenFatJet (and DarkHadronJet with
--dark-hadron-jets) with
    python3 JetRadiusOptimization/radius.py -i events.root
"""

from collections import Counter

import awkward as ak
import fastjet
import numpy as np
import vector
from scipy.optimize import linear_sum_assignment

from common import DelphesSchema2, get_constituents_kernel, load_events
from darkancestry import dark_hadron_ancestry
from Histogram import calc_axis1_axis2, calculate_girth, calculate_ptD, radius_percentile

vector.register_awkward()

# same as GenFatJet and DarkHadronJet in cards/delphes_card_CMS.tcl
PT_MIN = 15.0

# one row per (event, truth jet, truth radius, radius); truth_pt is the truth
# jet's visible pT, jet_pt the scalar pT sum of the matched jet's constituents
COLUMNS = (
    "event", "truth_jet", "truth_radius", "radius", "matched", "n_dh", "n_full",
    "f_full", "f_partial", "f_nodh", "acceptance", "purity", "quality", "iou",
    "truth_pt", "jet_pt",
    # jet shapes (Histogram.py): radii containing 90/95/99% of the pT, girth, pT^D, axes
    "truth_r90", "truth_r95", "truth_r99", "psi_truth",
    "jet_r90", "jet_r95", "jet_r99", "jet_girth", "jet_ptd", "jet_major", "jet_minor", "jet_nconst",
)

SHAPE_PCTS = (90, 95, 99)


def load(path):
    """Delphes events; DelphesSchema2 drops the fBits branches, which fastjet would otherwise read."""
    return load_events(path, schema=DelphesSchema2)


def cluster(particles, radius, pt_min=PT_MIN):
    """Anti-kT jets of the particles: constituent indices per event and jet."""
    jet_def = fastjet.JetDefinition(fastjet.antikt_algorithm, radius)
    # to_packed: fastjet clusters the whole underlying buffer, ignoring slices
    return fastjet.ClusterSequence(ak.to_packed(particles), jet_def).constituent_index(min_pt=pt_min)


def nudged(particles, nudge, rng):
    """The particles with each kinematic field multiplied by an independent factor 1 + nudge * N(0, 1)."""
    counts = ak.num(particles)
    n = int(ak.sum(counts))
    for field in ("PT", "Eta", "Phi", "Mass", "E", "Px", "Py", "Pz"):
        particles = ak.with_field(particles, particles[field] * ak.unflatten(1 + nudge * rng.standard_normal(n), counts), field)
    return particles


def event_metrics(descends, pt, no_dh, truth, jets, min_containment):
    """Rows of metrics for one event, one truth radius and one radius.

    descends: bool[n_dark_hadrons, n_particles], particle descends from dark hadron
    pt: particle pT; no_dh: particle has no dark-hadron ancestor
    truth: dark-hadron indices of each truth jet
    jets: constituent indices of each visible jet
    """
    if not truth:
        return []
    in_jet = np.zeros((len(jets), len(pt)), dtype=bool)
    for j, constituents in enumerate(jets):
        in_jet[j, constituents] = True
    jet_pt = in_jet @ pt

    # per (dark hadron, jet): how much of each dark hadron each jet captured
    captured_n = descends.astype(int) @ in_jet.T.astype(int)
    captured_pt = descends @ (in_jet * pt).T
    n_desc = descends.sum(axis=1)[:, None]
    if min_containment >= 1:
        full = captured_n == n_desc
    else:
        full = captured_pt >= min_containment * (descends @ pt)[:, None]
    full &= captured_n > 0
    partial = (captured_n > 0) & ~full
    with np.errstate(invalid="ignore", divide="ignore"):
        f_full = (full * captured_pt).sum(axis=0) / jet_pt
        f_partial = (partial * captured_pt).sum(axis=0) / jet_pt
        f_nodh = (in_jet * pt * no_dh).sum(axis=1) / jet_pt

    # one-to-one matching of truth jets to visible jets by shared visible pT
    truth_visible = np.array([descends[dh].any(axis=0) for dh in truth]).reshape(len(truth), -1)
    overlap = (truth_visible * pt) @ in_jet.T
    match = dict(zip(*linear_sum_assignment(-overlap))) if overlap.size else {}

    rows = []
    for t, dh in enumerate(truth):
        n_dh = int((n_desc[dh, 0] > 0).sum())
        if n_dh == 0:
            continue
        j = match.get(t)
        row = {"truth_jet": t, "jet": -1, "n_dh": n_dh, "truth_pt": float(truth_visible[t] @ pt)}
        if j is None or overlap[t, j] <= 0:
            row.update(matched=0, n_full=0, f_full=np.nan, f_partial=np.nan, f_nodh=np.nan,
                       acceptance=0.0, purity=np.nan, quality=0.0, iou=0.0, jet_pt=0.0)
        else:
            n_full = int(full[dh, j].sum())
            acceptance = (full[dh, j] * captured_pt[dh, j]).sum() / row["truth_pt"]
            iou = overlap[t, j] / (row["truth_pt"] + jet_pt[j] - overlap[t, j])
            row.update(jet=int(j), matched=1, n_full=n_full, f_full=f_full[j], f_partial=f_partial[j],
                       f_nodh=f_nodh[j], acceptance=acceptance, purity=f_full[j],
                       quality=acceptance * f_full[j], iou=iou, jet_pt=jet_pt[j])
        rows.append(row)
    return rows


def gather(particles, idx):
    """particles[event][idx[event][jet][k]] as an array[event][jet][k]."""
    counts = ak.to_numpy(ak.num(particles))
    offsets = np.concatenate([[0], np.cumsum(counts)[:-1]])
    flat = ak.flatten(particles)[ak.to_numpy(ak.flatten(idx + offsets, axis=None))]
    return ak.unflatten(ak.unflatten(flat, ak.flatten(ak.num(idx, axis=2))), ak.num(idx, axis=1))


def four_vectors(particles):
    return ak.zip({"px": particles.px, "py": particles.py, "pz": particles.pz, "E": particles.energy},
                  with_name="Momentum4D")


def jet_records(constituents, axis=None):
    """Jets with a Constituents field, as Histogram.py's shape functions expect.

    axis: four-vectors giving each jet's axis; default is the constituent sum
    """
    if axis is None:
        axis = ak.sum(constituents, axis=-1)
    return ak.zip({"px": axis.px, "py": axis.py, "pz": axis.pz, "E": axis.E, "Constituents": constituents},
                  depth_limit=2, with_name="Momentum4D")


def jet_shapes(jets, radii=()):
    """Per-jet shape variables (lists per event), using the Histogram.py definitions.

    DeltaR is in (eta, phi), as in Histogram.py. Jets without constituents give NaN.
    """
    # empty jets (0/0) and rounding in the axes (sqrt of a tiny negative) only give NaN
    with np.errstate(invalid="ignore", divide="ignore"):
        shapes = {f"r{pct}": radius_percentile(jets, pct) for pct in SHAPE_PCTS}
        shapes["girth"] = calculate_girth(jets)
        shapes["ptd"] = calculate_ptD(jets)
        shapes["major"], shapes["minor"] = calc_axis1_axis2(jets)
        shapes["nconst"] = ak.num(jets.Constituents, axis=-1)
        # integrated jet shape Psi(r): fraction of the constituent pT within r of the axis
        dr, pt = jets.deltaR(jets.Constituents), jets.Constituents.pt
        for r in radii:
            shapes[f"psi_{r}"] = ak.sum(pt * (dr < r), axis=-1) / ak.sum(pt, axis=-1)
    return {k: ak.to_list(ak.fill_none(v, np.nan)) for k, v in shapes.items()}


def radius_scan(events, radii, truth_radii=(0.8,), pt_min=PT_MIN, min_containment=0.95):
    """Metrics for every truth jet at every (truth radius, radius), as numpy columns."""
    # bit j of each GenCandidate's mask: descends from DarkHadronCandidate j
    masks = dark_hadron_ancestry(events)
    n_dh = ak.to_numpy(ak.num(events.DarkHadronCandidate))
    pt = events.GenCandidate.pt
    candidates, dark_hadrons = four_vectors(events.GenCandidate), four_vectors(events.DarkHadronCandidate)
    truth = {r: cluster(events.DarkHadronCandidate, r, pt_min) for r in truth_radii}
    visible = {r: cluster(events.GenCandidate, r, pt_min) for r in radii}

    def descendants(i):
        mask = ak.to_numpy(masks[i])
        bits = np.arange(n_dh[i], dtype=np.uint64)[:, None]
        return mask, ((mask[None, :] >> bits) & np.uint64(1)).astype(bool)

    # visible jets: their own shapes (axis = constituent sum)
    jet_shape = {r: jet_shapes(jet_records(gather(candidates, visible[r]))) for r in radii}
    # truth jets: shapes of their visible descendants around the dark-hadron jet axis
    truth_shape = {}
    for r_dh in truth_radii:
        visible_of_truth = ak.Array([
            [np.flatnonzero(d[dh].any(axis=0)) for dh in (ak.to_numpy(t) for t in truth[r_dh][i])]
            for i, d in ((i, descendants(i)[1]) for i in range(len(events)))
        ])
        axis = ak.sum(gather(dark_hadrons, truth[r_dh]), axis=-1)
        truth_shape[r_dh] = jet_shapes(jet_records(gather(candidates, visible_of_truth), axis), radii)

    rows = []
    for i in range(len(events)):
        mask, descends = descendants(i)
        event_pt = ak.to_numpy(pt[i]).astype(np.float64)
        for r_dh in truth_radii:
            event_truth = [ak.to_numpy(t) for t in truth[r_dh][i]]
            ts = truth_shape[r_dh]
            for r in radii:
                jets = [ak.to_numpy(j) for j in visible[r][i]]
                js = jet_shape[r]
                for row in event_metrics(descends, event_pt, mask == 0, event_truth, jets, min_containment):
                    t, j = row["truth_jet"], row.pop("jet")
                    row.update({f"truth_r{p}": ts[f"r{p}"][i][t] for p in SHAPE_PCTS})
                    row["psi_truth"] = ts[f"psi_{r}"][i][t]
                    for k in ("r90", "r95", "r99", "girth", "ptd", "major", "minor", "nconst"):
                        row[f"jet_{k}"] = js[k][i][j] if j >= 0 else np.nan
                    rows.append(dict(row, event=i, truth_radius=r_dh, radius=r))
    return {c: np.array([row[c] for row in rows], dtype=float) for c in COLUMNS}


def check_clustering(events, candidates, jets, radius=0.8, nudge=1e-6, trials=10):
    """Compare our clustering of events[candidates] with the stored events[jets].

    Delphes clustered the unrounded momenta; we cluster the float32 values
    stored in the file. In rare events anti-kT meets a near-tie that this
    rounding can tip. A mismatched event counts as "unstable" instead of a
    failure if nudging its inputs by a relative 1e-6 changes our own clustering.
    """
    def jet_sets(particles):
        return [Counter(tuple(sorted(ak.to_numpy(j))) for j in e) for e in cluster(particles, radius)]

    stats = {"jets": 0, "mismatched_jets": 0, "unstable_jets": 0}
    mismatched = {}
    for i, (mine, refs, uids) in enumerate(zip(jet_sets(events[candidates]), events[jets].Constituents.refs, events[candidates].fUniqueID)):
        stored = Counter(tuple(sorted(get_constituents_kernel(ak.to_numpy(r), ak.to_numpy(uids)))) for r in refs)
        stats["jets"] += len(refs)
        n_bad = max(sum(mine.values()), len(refs)) - sum((mine & stored).values())
        if n_bad:
            mismatched[i] = (n_bad, mine)
    rng = np.random.default_rng(1)
    for i, (n_bad, mine) in mismatched.items():
        unstable = any(jet_sets(nudged(events[candidates][i : i + 1], nudge, rng))[0] != mine for _ in range(trials))
        stats["unstable_jets" if unstable else "mismatched_jets"] += n_bad
    stats["passed"] = stats["mismatched_jets"] == 0
    return stats


if __name__ == "__main__":
    from magiconfig import ArgumentParser, ArgumentDefaultsRawHelpFormatter

    parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsRawHelpFormatter)
    parser.add_argument("-i", "--input", required=True, help="Delphes events.root")
    parser.add_argument("-n", "--max-events", type=int, default=None, help="number of events to check")
    parser.add_argument("--dark-hadron-jets", action="store_true", help="also check DarkHadronJet")
    args = parser.parse_args()

    events = load(args.input)
    if args.max_events is not None:
        events = events[: args.max_events]
    checks = [("GenCandidate", "GenFatJet")] + [("DarkHadronCandidate", "DarkHadronJet")] * args.dark_hadron_jets
    for candidates, jets in checks:
        print(jets, check_clustering(events, candidates, jets))
