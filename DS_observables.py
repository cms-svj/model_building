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
# FastJet's built-in ECF ratios, which need at least 3 (N2) or 4 (N3) constituents to be defined
ratio_ecfs = {
    "N2": dict(func="n2"),
    "N3": dict(func="n3"),
}
MIN_CONSTITUENTS = {"N2": 3, "N3": 4}

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

def getLundMultiplicity(lund, kt = 1):
    # Number of primary Lund-plane declusterings with kt above the cut
    # lund: the primary Lund-plane declusterings of each jet (see calculate_jet_substructure)
    kt_values = lund["kt"]
    mult = ak.sum(kt_values > kt, axis=-1)
    return mult

def getECF(cluster_seq, beta=1.0, **kwargs):
    # FastJet's energy correlation function of each jet when clustered exclusively into one jet
    # kwargs select the function (npoint=2 or 3 for e2 or e3, func="c2", "d2", "n2" or "n3")
    ecf = cluster_seq.exclusive_jets_energy_correlator(njets=1, beta=beta, **kwargs)
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


def calculate_shape_substructure(jet12, lund):
    """
    Substructure observables of the two leading jets of each event that are calculated from the jet constituents
    and the primary Lund-plane declusterings (no FastJet call here):
    girth, ptD, number of constituents, major and minor axis, primary Lund multiplicities and jet shape radii.
    jet12: the two leading jets of each event, padded to exactly two with None (shape: event -> jet)
    lund: the primary Lund-plane declusterings of each jet of jet12, flattened (see calculate_jet_substructure)
    Returns a dictionary with an array (shape: event -> jet, padded jets are None) for each name
    in JET_SUBSTRUCTURE that is not an ECF (see calculate_ecfs).
    """
    obs = {}
    obs["girth"] = calculate_girth(jet12)
    obs["ptD"] = calculate_ptD(jet12)
    # axis=2, not -1: reco constituents have a ragged Particles field, so -1 overshoots
    obs["nconst"] = ak.num(jet12.Constituents, axis=2)
    obs["majoraxis"], obs["minoraxis"] = calc_axis1_axis2(jet12)

    jet12_shape = ak.num(jet12, axis=1)
    for k in KT_CUTS:
        obs[f"lundMult{k}"] = ak.unflatten(getLundMultiplicity(lund, kt = k), jet12_shape)

    for pct in DR_PCTS:
        obs[f"radius{pct}"] = calculate_radius(jet12, pct)
    return obs


def calculate_ecfs(cs, jet12_shape, n_const, beta=1.0):
    """
    All the energy correlation functions (the ECFs in ecfs and N2, N3 in ratio_ecfs) of each jet, with getECF.
    cs: the FastJet cluster sequence of the jets (see calculate_jet_substructure)
    jet12_shape: number of jets of each event, including the padded ones
    n_const: number of constituents of each jet (flat)
    N2 and N3 are only defined for jets with at least MIN_CONSTITUENTS constituents, the others are NaN.
    Returns a dictionary with an array (shape: event -> jet) for each ECF.
    """
    obs = {}
    for name, kwargs in {**ratio_ecfs, **ecfs}.items():
        values = getECF(cs, beta, **kwargs)
        if name in MIN_CONSTITUENTS:
            values = ak.to_numpy(ak.flatten(ak.fill_none(values, np.nan), axis=None)).astype(np.float64)
            values[n_const < MIN_CONSTITUENTS[name]] = np.nan
        obs[name] = ak.unflatten(values, jet12_shape)
    return obs


def calculate_lund_points(lund, n_per_event, lund_R=0.8):
    """
    Primary Lund-plane points of each jet: LundX = log(lund_R/Delta), LundY = log(kt).
    lund: the primary Lund-plane declusterings of each jet (see calculate_jet_substructure)
    n_per_event: number of jets of each event, including the padded ones
    Returns a dictionary with "LundX", "LundY" (shape: event -> Lund point).
    """
    # jets without declusterings (padded jets, jets with a single constituent) have no points
    delta = ak.fill_none(lund["Delta"], [], axis=0)
    kt = ak.fill_none(lund["kt"], [], axis=0)

    good = np.isfinite(delta) & np.isfinite(kt) & (delta > 0.0) & (kt > 0.0)
    x = np.log(lund_R / delta[good])
    y = np.log(kt[good])

    n_points = ak.to_numpy(ak.num(x, axis=1))  # number of Lund points of each jet
    x_flat = ak.to_numpy(ak.flatten(x, axis=1))
    y_flat = ak.to_numpy(ak.flatten(y, axis=1))
    return {
        "LundX": lund_points_per_event(x_flat, n_points, n_per_event),
        "LundY": lund_points_per_event(y_flat, n_points, n_per_event),
    }


def calculate_jet_substructure(jets, beta=1.0, lund_R=0.8):
    """
    All the observables (see JET_SUBSTRUCTURE and LUND_OBSERVABLES) of the jets of one jet collection.
    jets: the two leading jets of each event (shape: event -> jet, at most two jets per event)

    The FastJet calls are all made here, once: the jets are clustered once (Cambridge/Aachen, one jet containing
    all the constituents), and the primary Lund-plane declusterings are calculated once. The results are passed to
    the functions that calculate the observables.

    Returns two dictionaries and the number of jets used for the Lund points (jets with at least 2 constituents):
        jet_sub: all the names in JET_SUBSTRUCTURE (shape: event -> jet, padded jets and undefined values are None
                 or NaN), written into the record of the collection
        lund_obs: "LundX", "LundY" (shape: event -> Lund point), written as separate branches
    """
    # the substructure code works on exactly two jets per event (padded with None)
    jet12 = ak.pad_none(jets, target=2, axis=1)
    jet12_shape = ak.num(jet12, axis=1)
    jets_flat = ak.flatten(jet12, axis=1)
    n_const = count_constituents(jets_flat)

    cs = fj_cluster_sequence(jets_flat)
    # firsts, not flatten: both strip the njets=1 axis, but flatten drops padded
    # None jets entirely, which breaks the unflatten back to (event, jet)
    lund = ak.firsts(cs.exclusive_jets_lund_declusterings(njets=1))

    jet_sub = calculate_shape_substructure(jet12, lund)
    jet_sub.update(calculate_ecfs(cs, jet12_shape, n_const, beta))
    lund_obs = calculate_lund_points(lund, ak.to_numpy(jet12_shape), lund_R)
    n_jets = int(np.count_nonzero(n_const >= 2))
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

        jet_sub, lund_obs, n_jets = calculate_jet_substructure(jets, beta=1.0, lund_R=0.8)

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
