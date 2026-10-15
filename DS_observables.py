import numpy as np
import awkward as ak
import fastjet
import uproot
import vector
from common import load_events

vector.register_awkward()

# value stored for jets that cannot be assigned an observable (too few constituents, etc.)
SENTINEL = -99.0

# branches of the friend tree that get merged into the Delphes tree (used by run_model and passed to merge_trees.C)
FRIEND_BRANCHES = ["FatJet_N2", "FatJet_N3", "FatJet_LundX", "FatJet_LundY"]


def to_f32_jagged(arr):
    # Fill missing whole events with empty lists: None -> []
    arr = ak.fill_none(arr, [], axis=0)

    # Fill missing jet-level values with 0.0: [None] -> [0.0]
    arr = ak.fill_none(arr, 0.0, axis=-1)

    arr = ak.values_astype(arr, np.float64)
    arr = ak.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    return ak.values_astype(arr, np.float32)


# Force the number of entries per event in an observable to match the number of FatJets
def force_fatjet_layout(arr, FJ):
    arr = to_f32_jagged(arr)

    # position of each FatJet within its event, e.g. [[0, 1], [], [0], ...]
    jet_index = ak.local_index(FJ, axis=1)

    # pad every event to the longest list, so that every FatJet position exists in every event
    longest = int(max(ak.max(ak.num(arr, axis=1), initial=0), ak.max(ak.num(FJ, axis=1), initial=0)))
    padded = ak.pad_none(arr, longest, axis=1)

    # keep exactly one entry per FatJet (extra values are dropped), and fill the padding with zeros
    fixed = ak.fill_none(padded[jet_index], 0.0)
    return ak.to_packed(ak.values_astype(fixed, np.float32))

def extract_ecf_value(ecf_array):
    flat = ak.to_numpy(ak.flatten(ecf_array, axis=None))
    if len(flat) == 0:
        return np.nan
    return float(flat[0])


def constituents_to_fastjet_array(constituents):
    particles = []

    for constituent in constituents:
        pt = float(constituent.PT)
        eta = float(constituent.Eta)
        phi = float(constituent.Phi)

        if hasattr(constituent, "Mass") and constituent.Mass is not None:
            mass = float(constituent.Mass)
        else:
            mass = 0.0

        px = pt * np.cos(phi)
        py = pt * np.sin(phi)
        pz = pt * np.sinh(eta)
        energy = np.sqrt(px**2 + py**2 + pz**2 + mass**2)

        particles.append({"px": px, "py": py, "pz": pz, "E": energy})

    # Outer list means one FastJet event.
    return ak.Array([particles])


def calculate_n2_n3_fastjet(jets, beta=1.0, fill_value=SENTINEL):
    """
    Generalized ECF observables of each FatJet:

        N2 = (2e3) / (1e2)^2
        N3 = (2e4) / (1e3)^2

    Output shape: event -> FatJet
    """
    all_n2 = []
    all_n3 = []

    # The constituents already belong to one FatJet. The large radius ensures that
    # FastJet combines all of them into one exclusive jet.
    jet_def = fastjet.JetDefinition(fastjet.cambridge_algorithm, 1000.0)

    for jet_group in jets:
        if jet_group is None or len(jet_group) == 0:
            all_n2.append([])
            all_n3.append([])
            continue

        event_n2 = []
        event_n3 = []

        for jet in jet_group:
            if jet is None:
                event_n2.append(fill_value)
                event_n3.append(fill_value)
                continue

            constituents = jet["Constituents"]

            if constituents is None:
                event_n2.append(fill_value)
                event_n3.append(fill_value)
                continue

            n_constituents = len(constituents)

            if n_constituents < 3:
                event_n2.append(fill_value)
                event_n3.append(fill_value)
                continue

            particles = constituents_to_fastjet_array(constituents)
            cluster = fastjet.ClusterSequence(particles, jet_def)

            # 1e2
            e2_1 = extract_ecf_value(
                cluster.exclusive_jets_energy_correlator(
                    njets=1, beta=beta, npoint=2, angles=1, func="generalized"
                )
            )

            # 2e3
            e3_2 = extract_ecf_value(
                cluster.exclusive_jets_energy_correlator(
                    njets=1, beta=beta, npoint=3, angles=2, func="generalized"
                )
            )

            if np.isfinite(e2_1) and np.isfinite(e3_2) and e2_1 > 0.0:
                n2 = e3_2 / (e2_1**2)
            else:
                n2 = fill_value

            if not np.isfinite(n2):
                n2 = fill_value

            event_n2.append(n2)

            if n_constituents < 4:
                event_n3.append(fill_value)
                continue

            # 1e3
            e3_1 = extract_ecf_value(
                cluster.exclusive_jets_energy_correlator(
                    njets=1, beta=beta, npoint=3, angles=1, func="generalized"
                )
            )

            # 2e4
            e4_2 = extract_ecf_value(
                cluster.exclusive_jets_energy_correlator(
                    njets=1, beta=beta, npoint=4, angles=2, func="generalized"
                )
            )

            if np.isfinite(e3_1) and np.isfinite(e4_2) and e3_1 > 0.0:
                n3 = e4_2 / (e3_1**2)
            else:
                n3 = fill_value

            if not np.isfinite(n3):
                n3 = fill_value

            event_n3.append(n3)

        all_n2.append(event_n2)
        all_n3.append(event_n3)

    N2 = ak.values_astype(ak.Array(all_n2), np.float32)
    N3 = ak.values_astype(ak.Array(all_n3), np.float32)

    return ak.to_packed(N2), ak.to_packed(N3)


def calculate_primary_lund_xy_fastjet(FJ, R=0.8, min_jet_pt=20.0):
    """
    Primary Lund-plane points of each FatJet, using FastJet.

    Output shape: event -> Lund point (all jets of an event are put together)
    """
    all_lund_x = []
    all_lund_y = []

    n_jets = 0

    jet_def = fastjet.JetDefinition(fastjet.cambridge_algorithm, R)

    for event_jets in FJ:
        event_lund_x = []
        event_lund_y = []

        if event_jets is None:
            all_lund_x.append(event_lund_x)
            all_lund_y.append(event_lund_y)
            continue

        for jet in event_jets:

            if float(jet.PT) < min_jet_pt:
                continue

            constituents = jet.Constituents

            if constituents is None or len(constituents) < 2:
                continue

            # This jet is included in the Lund calculation.
            n_jets += 1

            particle_records = []

            for constituent in constituents:
                mass = 0.0

                if hasattr(constituent, "Mass"):
                    mass_value = constituent.Mass

                    if mass_value is not None:
                        mass = float(mass_value)

                particle_records.append(
                    {
                        "pt": float(constituent.PT),
                        "eta": float(constituent.Eta),
                        "phi": float(constituent.Phi),
                        "M": mass,
                    }
                )

            # Outer list means one FastJet event.
            particles = ak.Array([particle_records], with_name="Momentum4D")

            cluster = fastjet.ClusterSequence(particles, jet_def)

            # Recluster the constituents into one exclusive C/A jet,
            # then obtain its primary Lund declustering sequence.
            lund = cluster.exclusive_jets_lund_declusterings(njets=1)

            if "Delta" not in lund.fields or "kt" not in lund.fields:
                raise RuntimeError(
                    "FastJet Lund output does not contain the expected "
                    f"'Delta' and 'kt' fields. Found: {lund.fields}"
                )

            delta = ak.to_numpy(ak.flatten(lund["Delta"], axis=None))
            kt = ak.to_numpy(ak.flatten(lund["kt"], axis=None))

            good = np.isfinite(delta) & np.isfinite(kt) & (delta > 0.0) & (kt > 0.0)

            delta = delta[good]
            kt = kt[good]

            event_lund_x.extend(np.log(R / delta).tolist())
            event_lund_y.extend(np.log(kt).tolist())

        all_lund_x.append(event_lund_x)
        all_lund_y.append(event_lund_y)

    LundX = to_f32_jagged(ak.Array(all_lund_x))
    LundY = to_f32_jagged(ak.Array(all_lund_y))

    return LundX, LundY, n_jets


def make_friend_tree(root_fname, friend_fname, with_constituents=True, verbose=False):
    """
    Compute the additional dark-sector observables of the FatJets in root_fname
    (N2, N3 and the primary Lund-plane points) and write them to a "friend" tree
    in friend_fname, with exactly one entry per entry of the Delphes tree.
    """
    if not with_constituents:
        raise ValueError("the dark-sector observables require jet constituents")

    events = load_events(root_fname, with_constituents=True)
    n_evt = len(events)
    FJ = events.FatJet

    # N2, N3
    N2, N3 = calculate_n2_n3_fastjet(FJ, beta=1.0, fill_value=SENTINEL)
    jets_like = {"N2": N2, "N3": N3}

    # primary Lund plane points. Shape: event -> Lund point
    FatJet_LundX, FatJet_LundY, n_jets = calculate_primary_lund_xy_fastjet(
        FJ, R=0.8, min_jet_pt=20.0
    )

    # make sure the friend tree has the same number of entries as the Delphes tree
    all_arrays = dict(jets_like, LundX=FatJet_LundX, LundY=FatJet_LundY)
    for name, arr in all_arrays.items():
        if len(arr) != n_evt:
            raise RuntimeError(
                f"{name} has {len(arr)} events but Delphes has {n_evt}: "
                "do not drop or reorder events when computing friend branches."
            )

    # force all jet-level observables to use the exact FatJet layout
    jets_like = {name: force_fatjet_layout(arr, FJ) for name, arr in jets_like.items()}

    # Zip the jet-level variables into one FatJet record, so there is only one counter
    # per group structure (event -> jet, and event -> Lund point)
    FatJet_record = ak.zip(
        {"N2": jets_like["N2"], "N3": jets_like["N3"]},
        depth_limit=2,
    )

    with uproot.recreate(friend_fname) as fout:
        fout.mktree(
            "DelphesFriend",
            {
                "FatJet": ak.type(FatJet_record),
                "FatJet_LundX": ak.type(FatJet_LundX),
                "FatJet_LundY": ak.type(FatJet_LundY),
            },
        )

        fout["DelphesFriend"].extend(
            {
                "FatJet": FatJet_record,
                "FatJet_LundX": FatJet_LundX,
                "FatJet_LundY": FatJet_LundY,
            }
        )

    if verbose:
        print(f"Wrote friend tree with {n_evt} events and {n_jets} Lund jets: {friend_fname}")
