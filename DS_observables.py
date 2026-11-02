import numpy as np
import awkward as ak
import fastjet
import uproot
import vector
from common import load_events

vector.register_awkward()

# value stored for jets that cannot be assigned an observable (too few constituents, etc.)
SENTINEL = -99.0

# jet collections for which the observables are calculated: only the two leading jets of each collection are used
COLLECTIONS = ["FatJet", "DarkPartonJet", "DarkHadronJet", "DarkHadronVisibleJet", "DarkHadronStableJet"]

# settings of the substructure observables
KT_CUTS = [1, 2, 5, 10]  # kt cuts in GeV of the primary Lund multiplicities
DR_PCTS = [90, 95, 99]  # percentages of the pT of the jet shape radii
# arguments of FastJet's energy_correlator
ecfs = {
    "e2": dict(npoint=2),
    "e3": dict(npoint=3),
    "C2": dict(func="c2"),
    "D2": dict(func="d2"),
}

# observables of each jet (shape: event -> jet), written into one record per collection
JET_SUBSTRUCTURE = (
    ["N2", "N3", "girth", "ptD", "nconst", "majoraxis", "minoraxis"]
    + [f"lundMult{k}" for k in KT_CUTS]
    + list(ecfs)
    + [f"radius{pct}" for pct in DR_PCTS]
)
# primary Lund-plane points of each collection (shape: event -> Lund point), written as separate branches
LUND_OBSERVABLES = ["LundX", "LundY"]

# branches of the friend tree that get merged into the Delphes tree (used by run_model and passed to merge_trees.C)
FRIEND_BRANCHES = [f"{coll}_{name}" for coll in COLLECTIONS for name in JET_SUBSTRUCTURE + LUND_OBSERVABLES]


def to_f32_jagged(arr, fill_value=0.0):
    # Fill missing whole events with empty lists: None -> []
    arr = ak.fill_none(arr, [], axis=0)

    # Fill missing jet-level values with fill_value: [None] -> [fill_value]
    arr = ak.fill_none(arr, fill_value, axis=-1)

    arr = ak.values_astype(arr, np.float64)
    arr = ak.nan_to_num(arr, nan=fill_value, posinf=fill_value, neginf=fill_value)
    return ak.values_astype(arr, np.float32)


# Force the number of entries per event in an observable to match the number of jets
def force_jet_layout(arr, jets, fill_value=0.0):
    arr = to_f32_jagged(arr, fill_value)

    # position of each jet within its event, e.g. [[0, 1], [], [0], ...]
    jet_index = ak.local_index(jets, axis=1)

    # pad every event to the longest list, so that every jet position exists in every event
    longest = int(max(ak.max(ak.num(arr, axis=1), initial=0), ak.max(ak.num(jets, axis=1), initial=0)))
    padded = ak.pad_none(arr, longest, axis=1)

    # keep exactly one entry per jet (extra values are dropped), and fill the padding with fill_value
    fixed = ak.fill_none(padded[jet_index], fill_value)
    return ak.to_packed(ak.values_astype(fixed, np.float32))


def count_constituents(jets):
    # number of constituents of each jet in a flat jet array
    return ak.to_numpy(ak.fill_none(ak.num(jets["Constituents"], axis=1), 0))


def cluster_jets(jets, R):
    # one FastJet cluster sequence per jet, made directly from the (vector) constituents
    jet_def = fastjet.JetDefinition(fastjet.cambridge_algorithm, R)
    return fastjet.ClusterSequence(ak.to_packed(jets["Constituents"]), jet_def)


def builtin_ecf(cs, func, beta):
    # FastJet's built-in ECF variation (e.g. "n2", "n3") of each jet when clustered exclusively into one jet
    values = cs.exclusive_jets_energy_correlator(njets=1, beta=beta, func=func)
    values = ak.fill_none(values, np.nan)
    return ak.to_numpy(ak.flatten(values, axis=None)).astype(np.float64)


def lund_points_per_event(values, n_per_jet, n_per_event):
    """
    Flat list of Lund points -> shape: event -> Lund point (float32).
    n_per_jet: number of points of each jet, n_per_event: number of jets of each event.
    All the points of the jets of an event are put together.
    """
    per_jet = ak.unflatten(values, n_per_jet)  # jet -> point
    per_event = ak.unflatten(per_jet, n_per_event)  # event -> jet -> point
    return to_f32_jagged(ak.flatten(per_event, axis=2))  # event -> point


# jet substructure observables

def deltaR(jet):
    return jet.deltaR(jet.Constituents)

def calculate_girth(jet):
    particle_dR = deltaR(jet)
    girth = ak.sum(jet.Constituents.pt * particle_dR, axis=-1)
    girth = np.divide(girth,jet.pt) #normalize wrt jet pt
    return girth

def calculate_ptD(jet):
    sum_pt = ak.sum(jet.Constituents.pt,axis=-1)
    sum_pt2 = ak.sum(jet.Constituents.pt ** 2,axis=-1)
    ptD = np.sqrt(sum_pt2) / sum_pt

    # Return the result (ptD)
    return ptD

def calc_axis1_axis2(jet):
    jet_constpt = jet.Constituents.pt
    deta_particle = np.abs(jet.eta-jet.Constituents.eta)
    dphi_particle = np.abs(jet.deltaphi(jet.Constituents))

    # Calculate weights (pt^2) for each constituent
    weights_pt = jet_constpt**2

    # Calculate weighted sums for each event
    sum_weight = ak.sum(weights_pt, axis=-1)  # Sum of weights (pt^2) for each event

    sum_deta = ak.sum(deta_particle * weights_pt, axis=-1)
    sum_dphi = ak.sum(dphi_particle * weights_pt, axis=-1)
    sum_deta2 = ak.sum(deta_particle**2 * weights_pt, axis=-1)
    sum_dphi2 = ak.sum(dphi_particle**2 * weights_pt, axis=-1)
    sum_detadphi = ak.sum(deta_particle * dphi_particle * weights_pt, axis=-1)

    # Calculate averages
    ave_deta = sum_deta / sum_weight
    ave_dphi = sum_dphi / sum_weight
    ave_deta2 = sum_deta2 / sum_weight
    ave_dphi2 = sum_dphi2 / sum_weight

    # Calculate covariance matrix components
    a = ave_deta2 - ave_deta**2
    b = ave_dphi2 - ave_dphi**2
    c = -(sum_detadphi / sum_weight - ave_deta * ave_dphi)

    # Calculate the discriminant (delta) for each jet
    delta = np.sqrt(np.abs((a - b)**2 + 4 * c**2))

    # Calculate axis1 (major) and axis2 (minor) for each jet
    axis1 = np.sqrt(0.5 * (a + b + delta))
    axis2 = np.sqrt(0.5 * (a + b - delta))

    return axis1, axis2

def fj_cluster_sequence(jets):
    # Re-cluster the constituents with Cambridge/Aachen.
    # Using a very large R guarantees a single C/A jet containing all constituents
    jet_def = fastjet.JetDefinition(fastjet.cambridge_algorithm, 1000.0)
    # Jets with zero constituents (e.g. DHV jets whose dark hadrons were all stable) cause issues
    # to_packed is required: fastjet clusters the whole underlying buffer, ignoring masks/slices
    # so the empty jets must be removed from it.
    # Masking via indices (not ak.mask on jets directly) guarantees an
    # IndexedOptionArray, which to_packed compacts even when the input has no Nones yet
    nonempty = ak.fill_none(ak.num(jets.Constituents, axis=1) > 0, False)
    jets = jets[ak.mask(ak.local_index(jets, axis=0), nonempty)]
    cs = fastjet.ClusterSequence(ak.to_packed(jets.Constituents), jet_def)
    return cs

def getLundMultiplicity(cluster_seq, kt = 1):
    # Retrieve the primary Lund-plane declusterings for the single jet and apply kt cut
    lund = cluster_seq.exclusive_jets_lund_declusterings(njets=1)
    # firsts, not flatten: both strip the njets=1 axis, but flatten drops padded
    # None jets entirely, which breaks the unflatten back to (event, jet)
    kt_values = ak.firsts(lund)["kt"]
    mult = ak.sum(kt_values > kt, axis=-1)
    return mult

def getECF(cluster_seq, **kwargs):
    ecf = cluster_seq.exclusive_jets_energy_correlator(njets=1, beta=1, **kwargs)
    return ecf

def jet_const_cumsum(array):
    counts = ak.num(array, axis=-1)
    flat_counts = ak.flatten(counts, axis=None)
    flat_array = ak.flatten(array, axis=None)
    global_cumsum = np.cumsum(flat_array)
    offsets = np.zeros(len(flat_counts)+1, dtype=int)
    offsets[1:] = np.cumsum(flat_counts)
    start_sums = np.zeros_like(global_cumsum)
    prev_tot = np.zeros(len(flat_counts))
    prev_tot[1:] = global_cumsum[offsets[1:-1]-1]
    subtractions = np.repeat(prev_tot, flat_counts)
    # unflatten in two stages
    per_jet_flat = global_cumsum - subtractions
    jets_unflat = ak.unflatten(per_jet_flat, flat_counts)
    return ak.unflatten(jets_unflat, ak.num(counts, axis=1))


def calculate_radius(jet12, pct):
    """
    Jet shape radius of each jet: the pT-weighted percentile of the DeltaR between the jet and its constituents,
    i.e. the DeltaR within which pct % of the scalar sum of the constituent pT is contained.
    jet12: the two leading jets of each event, padded to exactly two with None (shape: event -> jet)
    """
    none_mask = ak.is_none(jet12, axis=1)
    const_prop = ak.zip({
        'dr': deltaR(jet12),
        'pt': jet12.Constituents.pt,
    })[~none_mask]
    sort_indices = ak.argsort(const_prop.dr, axis=-1)
    sorted_prop = const_prop[sort_indices]
    cumul_pt = jet_const_cumsum(sorted_prop.pt)
    # last element of sum is total from all constituents
    total_pt = ak.fill_none(ak.pad_none(cumul_pt, 1, axis=-1)[:, :, -1], 0)
    target_pt = pct/100 * total_pt
    mask_pt = cumul_pt >= target_pt
    return ak.pad_none(ak.firsts(sorted_prop[mask_pt].dr, axis=-1), target=2, axis=1)


def calculate_shape_substructure(jet12):
    """
    Substructure observables of the two leading jets of each event, calculated from the jet constituents:
    girth, ptD, number of constituents, major and minor axis, primary Lund multiplicities, ECFs and jet shape radii.
    jet12: the two leading jets of each event, padded to exactly two with None (shape: event -> jet)
    Returns a dictionary with an array (shape: event -> jet, padded jets are None) for each name in
    JET_SUBSTRUCTURE, except for N2 and N3, which are calculated in calculate_n2_n3_lund.
    """
    obs = {}
    obs["girth"] = calculate_girth(jet12)
    obs["ptD"] = calculate_ptD(jet12)
    # axis=2, not -1: reco constituents have a ragged Particles field, so -1 overshoots
    obs["nconst"] = ak.num(jet12.Constituents, axis=2)
    obs["majoraxis"], obs["minoraxis"] = calc_axis1_axis2(jet12)

    jet12_shape = ak.num(jet12, axis=1)
    cs = fj_cluster_sequence(ak.flatten(jet12, axis=1))
    for k in KT_CUTS:
        obs[f"lundMult{k}"] = ak.unflatten(getLundMultiplicity(cs, kt = k), jet12_shape)
    for name, kwargs in ecfs.items():
        obs[name] = ak.unflatten(getECF(cs, **kwargs), jet12_shape)

    for pct in DR_PCTS:
        obs[f"radius{pct}"] = calculate_radius(jet12, pct)
    return obs


def calculate_n2_n3_lund(event_jets, beta=1.0, fill_value=SENTINEL, lund_R=0.8):
    """
    N2, N3 and the primary Lund-plane points of each jet, using FastJet.
    All jets are clustered once (Cambridge/Aachen), and all observables are taken from that clustering.

        N2: FastJet's built-in "n2", needs at least 3 constituents
        N3: FastJet's built-in "n3", needs at least 4 constituents
        Lund points: LundX = log(lund_R/Delta), LundY = log(kt), need at least 2 constituents

    Jets that cannot be assigned an N2 or N3 value get fill_value.
    Returns two dictionaries and the number of jets used for the Lund points:
        jet_sub: "N2", "N3" (shape: event -> jet), written into the record of the collection
        lund_obs: "LundX", "LundY" (shape: event -> Lund point), written as separate branches
    """
    jets = ak.flatten(event_jets, axis=1)  # all jets of all events in one list
    n_const = count_constituents(jets)

    n2 = np.full(len(jets), fill_value, dtype=np.float64)
    n3 = np.full(len(jets), fill_value, dtype=np.float64)
    n_points = np.zeros(len(jets), dtype=np.int64)  # number of Lund points of each jet
    x_flat = np.zeros(0, dtype=np.float64)
    y_flat = np.zeros(0, dtype=np.float64)
    n_jets = 0

    # Cluster once
    clustered = np.flatnonzero(n_const >= 2)
    if len(clustered) > 0:
        cs = cluster_jets(jets[clustered], R=1000.0)
        nc = n_const[clustered]


        n2_all = builtin_ecf(cs, "n2", beta)
        n3_all = builtin_ecf(cs, "n3", beta)
        ok2 = (nc >= 3) & np.isfinite(n2_all)
        ok3 = (nc >= 4) & np.isfinite(n3_all)
        n2[clustered[ok2]] = n2_all[ok2]
        n3[clustered[ok3]] = n3_all[ok3]

        # Primary Lund declustering
        lund = cs.exclusive_jets_lund_declusterings(njets=1)
        delta = lund["Delta"]
        kt = lund["kt"]

        # make sure there is exactly one list of declusterings per jet
        while delta.ndim > 2:
            delta = ak.flatten(delta, axis=1)
            kt = ak.flatten(kt, axis=1)

        good = np.isfinite(delta) & np.isfinite(kt) & (delta > 0.0) & (kt > 0.0)
        x = np.log(lund_R / delta[good])
        y = np.log(kt[good])

        n_jets = len(clustered)
        n_points[clustered] = ak.to_numpy(ak.num(x, axis=1))
        x_flat = ak.to_numpy(ak.flatten(x, axis=1))
        y_flat = ak.to_numpy(ak.flatten(y, axis=1))

    n_per_event = ak.to_numpy(ak.fill_none(ak.num(event_jets, axis=1), 0))

    # back to event -> jet
    N2 = ak.to_packed(ak.unflatten(n2.astype(np.float32), n_per_event))
    N3 = ak.to_packed(ak.unflatten(n3.astype(np.float32), n_per_event))

    # shape: event -> Lund point
    LundX = lund_points_per_event(x_flat, n_points, n_per_event)
    LundY = lund_points_per_event(y_flat, n_points, n_per_event)

    # per-jet values (event -> jet) go into the record of the collection
    jet_sub = {"N2": N2, "N3": N3}
    # separate branches: the number of Lund points differs from the number of jets
    lund_obs = {"LundX": LundX, "LundY": LundY}

    return jet_sub, lund_obs, n_jets


def calculate_jet_substructure(jets, beta=1.0, fill_value=SENTINEL, lund_R=0.8):
    """
    All the observables (see JET_SUBSTRUCTURE and LUND_OBSERVABLES) of the jets of one jet collection.
    jets: the two leading jets of each event (shape: event -> jet, at most two jets per event)
    Returns the same as calculate_n2_n3_lund, with the substructure observables added to jet_sub.
    """
    jet_sub, lund_obs, n_jets = calculate_n2_n3_lund(jets, beta=beta, fill_value=fill_value, lund_R=lund_R)
    # the substructure code works on exactly two jets per event (padded with None)
    jet_sub.update(calculate_shape_substructure(ak.pad_none(jets, target=2, axis=1)))
    return jet_sub, lund_obs, n_jets


def make_friend_tree(root_fname, friend_fname, with_constituents=True, verbose=False):
    """
    Compute the additional dark-sector observables of the two leading jets of each jet collection in COLLECTIONS
    (N2, N3, jet substructure and the primary Lund-plane points) in root_fname and write them to a "friend" tree
    in friend_fname, with exactly one entry per entry of the Delphes tree.
    """
    if not with_constituents:
        raise ValueError("the dark-sector observables require jet constituents")

    events = load_events(root_fname, with_constituents=True)
    n_evt = len(events)

    # branches to write: one record with the jet-level observables per collection,
    # plus one separate branch per Lund observable
    branches = {}
    for coll in COLLECTIONS:
        # only the two leading jets of each event are used
        jets = events[coll][:, 0:2]

        jet_sub, lund_obs, n_jets = calculate_jet_substructure(jets, beta=1.0, fill_value=SENTINEL, lund_R=0.8)

        # make sure the friend tree has the same number of entries as the Delphes tree
        for name, arr in {**jet_sub, **lund_obs}.items():
            if len(arr) != n_evt:
                raise RuntimeError(
                    f"{coll}_{name} has {len(arr)} events but Delphes has {n_evt}: "
                    "do not drop or reorder events when computing friend branches."
                )

        # force all jet-level observables to use the exact jet layout (None and NaN become the sentinel)
        jet_sub = {name: force_jet_layout(jet_sub[name], jets, fill_value=SENTINEL) for name in JET_SUBSTRUCTURE}

        # Zip the jet-level variables into one record per collection, so there is only one counter
        # per group structure (event -> jet, and event -> Lund point)
        branches[coll] = ak.zip(jet_sub, depth_limit=2)
        for name in LUND_OBSERVABLES:
            branches[f"{coll}_{name}"] = lund_obs[name]

        if verbose:
            print(f"{coll}: {n_jets} jets used for the Lund points")

    with uproot.recreate(friend_fname) as fout:
        fout.mktree("DelphesFriend", {name: ak.type(arr) for name, arr in branches.items()})
        fout["DelphesFriend"].extend(branches)

    if verbose:
        print(f"Wrote friend tree with {n_evt} events: {friend_fname}")
