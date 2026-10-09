#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from coffea.nanoevents import DelphesSchema

# import matplotlib.ticker as ticker
# import pandas as pd
import gc
import time

import fastjet
import vector



DelphesSchema.mixins.update({
    "GenParticle": "Particle",
    "GenCandidate": "Particle",
    "ParticleFlowCandidate": "Particle",
    "DarkHadronCandidate": "Particle",
    "FatJet": "Jet",
    "GenFatJet": "Jet",
    "DarkHadronJet": "Jet",
})


# In[ ]:


import os, sys

# Notebook is run from scripts/; set up the same layout the Condor job has
if os.path.basename(os.getcwd()) == "scripts":
    SCRIPTS_DIR = os.getcwd()
    REPO_ROOT = os.path.abspath(os.path.join(SCRIPTS_DIR, "..", ".."))
    sys.path.insert(0, REPO_ROOT)     # so svjHelper and other root modules are found
    sys.path.insert(0, SCRIPTS_DIR)   # so scripts/common.py wins over the root common.py
    os.chdir(REPO_ROOT)               # common.py uses relative "models/..." paths


# In[ ]:


# Have this cell read the models from the models file instead of harcoding them 
# like in the cell below 

import glob, os

def pick_models(glob_pattern):
    matches = sorted(glob.glob(glob_pattern), key=os.path.getmtime)
    if not matches:
        raise FileNotFoundError(f"No models matched: {glob_pattern}")
    return matches

# e.g. analyze all snowmass cmslike models, any rinv
model_dirs = pick_models("models/s-channel_*_rinv-*")

samples = [{"name": os.path.basename(d), "model": os.path.basename(d)} for d in model_dirs]


# In[ ]:


# Hardcode the names of the models in here (mostly for internal testing)

# samples = [
#     {"name": "CMS", "model": "s-channel_mmed-1000_Nc-2_Nf-2_scale-35.1539_mq-10_mpi-20_mrho-20_pvector-0.75_spectrum-cms_rinv-0.3"},
#      {"name": "Snowmass (CMS-like)", "model": "s-channel_mmed-1000_Nc-3_Nf-3_scale-33.3333_mq-33.73_mpi-20_mrho-83.666_pvector-0.5_spectrum-snowmass_cmslike_rinv-0.3"}
# ]

# samples = [{"name": "CMS", "model": "s-channel_mmed-1000_Nc-2_Nf-2_scale-35.1539_mq-10_mpi-20_mrho-20_pvector-0.75_spectrum-cms_rinv-0.3"}]


# samples = [{"name": "Snowmass (CMS-like)", "model": "s-channel_mmed-1000_Nc-3_Nf-3_scale-33.3333_mq-33.73_mpi-20_mrho-83.666_pvector-0.5_spectrum-snowmass_cmslike_rinv-0.3"}]



# samples = [{"name": "Master Model", "model": "s-channel_mmed-1000_Nc-3_Nf-3_scale-33.3333_mq-33.73_mpi-20_mrho-83.666_pvector-0.5_spectrum-INDEPENDENTmodel_rinv-0.333333"}]


# In[ ]:


# Sentinel for values in arrays that are "None" or zero (Used mostly for plotting)
 
SENTINEL = -99.0


# In[ ]:





# In[ ]:


# Start the counter to see how long it takes to load the constituents

print()
print()

start = time.perf_counter()
print("Start loading events")

print()
print()


# In[ ]:


# Load constituents

from common import load_sample 
# from svjHelper import svjHelper

for sample in samples:
    load_sample(sample, with_constituents=True) 
#     sample["num_events"] = ak.num(sample["events"], axis=0)


import common
print("[check] common loaded from:", common.__file__)


# In[ ]:


print()
print()

end = time.perf_counter()
total_seconds = end - start
print(f"Finished loading events in {total_seconds:.3f} seconds ~{total_seconds/60:.3f} minutes")

print()
print()


# In[ ]:





# In[ ]:


# Pretty printer for each event and jet


# Kinematics and girth calculation 

def normalize_angle(angle):
    angle = np.mod(angle, 2 * np.pi)
    angle = np.where(angle >= np.pi, angle - 2 * np.pi, angle)
    return angle

def deltaR(jet):
    deta_particle = np.abs(jet["Eta"] - jet["Constituents"]["Eta"])
    dphi_particle = np.abs(normalize_angle(jet["Phi"] - jet["Constituents"]["Phi"]))
    return np.sqrt(deta_particle**2 + dphi_particle**2)

def calculate_girth(jet):
    particle_dR = deltaR(jet)
    girth = ak.sum(jet["Constituents"]["PT"] * particle_dR, axis=-1)
    return girth / jet["PT"]



def print_jets_info(jets, girth, sentinel=-99.0):

    for ievt, (evt_jets, evt_girth) in enumerate(zip(jets, girth)):

        # --- Handle missing jets ---
        if evt_jets is None:

            print(f"\n=== Event {ievt} : 1 jet (EMPTY → sentinel) ===")
            print(
                f"  Jet  0: "
                f"eta = {sentinel: .3f}, "
                f"phi = {sentinel: .3f}, "
                f"pT = {sentinel: .3f} GeV, "
                f"girth = {sentinel: .4f}"
            )
            continue


        # --- Handle missing girth ---
        if evt_girth is None:
            evt_girth = [sentinel] * len(evt_jets)


        # --- Normal case ---
        print(f"\n=== Event {ievt} : {len(evt_jets)} jets ===")

        for ijet, (eta, phi, pt, g) in enumerate(
            zip(evt_jets.Eta, evt_jets.Phi, evt_jets.PT, evt_girth)
        ):
            print(
                f"  Jet {ijet:2d}: "
                f"eta = {eta: .3f}, "
                f"phi = {phi: .3f}, "
                f"pT = {pt: .3f} GeV, "
                f"girth = {g: .4f}"
            )



for sample in samples:
    events = sample["events"]
    jets = events.FatJet  
    girth = calculate_girth(jets)
    
    
    # print pretty jet info
    print_jets_info(jets, girth)  


# In[ ]:





# In[ ]:





# In[ ]:





# <br>
# 
# **Calculate Energy Correlation functions ECF**
# 
# <br>

# <br>
# 
# $$
# e_2^{(\beta)} = \sum_{1 \leq i < j \leq n_J} z_i z_j \, \theta_{ij}^{\beta}
# $$
# 
# $$
# e_3^{(\beta)} = \sum_{1 \leq i < j < k \leq n_J} z_i z_j z_k \, \theta_{ij}^{\beta} \theta_{ik}^{\beta} \theta_{jk}^{\beta}
# $$
# 
# $$
# e_4^{(\beta)} = \sum_{1 \leq i < j < k < \ell \leq n_J} z_i z_j z_k z_\ell \, \theta_{ij}^{\beta} \theta_{ik}^{\beta} \theta_{jk}^{\beta} \theta_{i\ell}^{\beta} \theta_{j\ell}^{\beta} \theta_{k\ell}^{\beta}
# $$
# 
# 
# $$
# z_i \equiv \frac{p_{Ti}}{\sum_{j \in \text{jet}} p_{Tj}},
# \quad
# \theta_{ij}^2 \equiv R_{ij}^2 = (\phi_i - \phi_j)^2 + (y_i - y_j)^2
# $$
# 
# 
# <br>
# 
# 

# **CALCULATION OF ECF WITH FASTJET**
# <br>
# <br>
# <br>

# In[ ]:


print()
print()

start = time.perf_counter()
print("Starting ECFs calculation")

print()
print()


# In[ ]:


def extract_ecf_value(ecf_array):
    flat = ak.to_numpy(
        ak.flatten(ecf_array, axis=None)
    )

    if len(flat) == 0:
        return np.nan

    return float(flat[0])


# In[ ]:


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

        energy = np.sqrt(
            px**2
            + py**2
            + pz**2
            + mass**2
        )

        particles.append(
            {
                "px": px,
                "py": py,
                "pz": pz,
                "E": energy,
            }
        )

    # Outer list means one FastJet event.
    return ak.Array([particles])


# In[ ]:


def calculate_n2_n3_fastjet(
    jets,
    beta=1.0,
    fill_value=SENTINEL,
):
    """
    Calculate generalized ECF observables:

        N2 = (_2 e3) / (_1 e2)^2
        N3 = (_2 e4) / (_1 e3)^2

    Output shape:
        event -> FatJet
    """

    all_n2 = []
    all_n3 = []

    # The constituents already belong to one FatJet.
    # The large radius ensures that FastJet combines all of them
    # into one exclusive jet.
    jet_def = fastjet.JetDefinition(
        fastjet.cambridge_algorithm,
        1000.0,
    )

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

            particles = constituents_to_fastjet_array(
                constituents
            )

            cluster = fastjet.ClusterSequence(
                particles,
                jet_def,
            )

            # _1 e2
            e2_1 = extract_ecf_value(
                cluster.exclusive_jets_energy_correlator(
                    njets=1,
                    beta=beta,
                    npoint=2,
                    angles=1,
                    func="generalized",
                )
            )

            # _2 e3
            e3_2 = extract_ecf_value(
                cluster.exclusive_jets_energy_correlator(
                    njets=1,
                    beta=beta,
                    npoint=3,
                    angles=2,
                    func="generalized",
                )
            )

            if (
                np.isfinite(e2_1)
                and np.isfinite(e3_2)
                and e2_1 > 0.0
            ):
                n2 = e3_2 / (e2_1**2)
            else:
                n2 = fill_value

            if not np.isfinite(n2):
                n2 = fill_value

            event_n2.append(n2)

            if n_constituents < 4:
                event_n3.append(fill_value)
                continue

            # _1 e3
            e3_1 = extract_ecf_value(
                cluster.exclusive_jets_energy_correlator(
                    njets=1,
                    beta=beta,
                    npoint=3,
                    angles=1,
                    func="generalized",
                )
            )

            # _2 e4
            e4_2 = extract_ecf_value(
                cluster.exclusive_jets_energy_correlator(
                    njets=1,
                    beta=beta,
                    npoint=4,
                    angles=2,
                    func="generalized",
                )
            )

            if (
                np.isfinite(e3_1)
                and np.isfinite(e4_2)
                and e3_1 > 0.0
            ):
                n3 = e4_2 / (e3_1**2)
            else:
                n3 = fill_value

            if not np.isfinite(n3):
                n3 = fill_value

            event_n3.append(n3)

        all_n2.append(event_n2)
        all_n3.append(event_n3)

    N2 = ak.values_astype(
        ak.Array(all_n2),
        np.float32,
    )

    N3 = ak.values_astype(
        ak.Array(all_n3),
        np.float32,
    )

    return ak.to_packed(N2), ak.to_packed(N3)


# In[ ]:


beta = 1.0

for sample in samples:
    events = sample["events"]
    jets = events.FatJet

    N2, N3 = calculate_n2_n3_fastjet(
        jets,
        beta=beta,
        fill_value=SENTINEL,
    )

    events["N2"] = N2
    events["N3"] = N3

    sample["events"] = events


# In[ ]:


plt.figure(figsize=(14, 5))

colors = ["blue", "red"]

for sample, color in zip(samples, colors):
    n2_values = ak.to_numpy(
        ak.flatten(
            sample["events"]["N2"],
            axis=None,
        )
    )

    valid = (
        np.isfinite(n2_values)
        & (n2_values != SENTINEL)
    )

    n2_values = n2_values[valid]

    plt.hist(
        n2_values,
        bins=50,
        histtype="step",
        edgecolor=color,
        linewidth=2,
        density=True,
        label=sample["name"],
    )

plt.xlabel(
    r"$N_2^{(1)}="
    r"\frac{{}_2e_3^{(1)}}"
    r"{\left({}_1e_2^{(1)}\right)^2}$",
    fontsize=14,
)

plt.ylabel("Normalized number of jets")
plt.title(r"$N_2^{(1)}$ Distribution")
plt.yscale("log")
plt.grid(False)
plt.legend()
plt.show()


# In[ ]:


plt.figure(figsize=(14, 5))

colors = ["blue", "red"]

for sample, color in zip(samples, colors):
    n3_values = ak.to_numpy(
        ak.flatten(
            sample["events"]["N3"],
            axis=None,
        )
    )

    valid = (
        np.isfinite(n3_values)
        & (n3_values != SENTINEL)
    )

    n3_values = n3_values[valid]

    plt.hist(
        n3_values,
        bins=50,
        range=(0, 4),
        histtype="step",
        edgecolor=color,
        linewidth=2,
        density=True,
        label=sample["name"],
    )

plt.xlabel(
    r"$N_3^{(1)}="
    r"\frac{{}_2e_4^{(1)}}"
    r"{\left({}_1e_3^{(1)}\right)^2}$",
    fontsize=14,
)

plt.ylabel("Normalized number of jets")
plt.title(r"$N_3^{(1)}$ Distribution")
plt.yscale("log")
plt.grid(False)
plt.legend()
plt.show()


# In[ ]:


print()
print()

end = time.perf_counter()
total_seconds = end - start
print(f"Finished ECF calculation and plotting in {total_seconds:.3f} seconds  ~{total_seconds/60:.3f} minutes")

print()
print()


# In[ ]:


del N2, N3, n2_values, n3_values
gc.collect()


# **Lund Plane**
# <br>
# <br>
# <br>

# In[ ]:


# Prepare arrays so they are safe for ROOT
# and also for the friend tree later in this code

def to_f32_jagged(arr):
    # Fill missing whole events with empty lists: None -> []
    arr = ak.fill_none(arr, [], axis=0)

    # Fill missing jet-level values with 0.0: [None] -> [0.0]
    arr = ak.fill_none(arr, 0.0, axis=-1)

    arr = ak.values_astype(arr, np.float64)
    arr = ak.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    return ak.values_astype(arr, np.float32)


# In[ ]:


def plot_primary_lund_image(
    H,
    x_edges,
    y_edges,
    title="Primary Lund Plane",
    text_label=None,
    save=False,
    filename="primary_lund_plane.png",
    vmax=None
):
    fig, ax = plt.subplots(figsize=(7, 6))

    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        H.T,
        shading="auto",
        cmap="RdBu_r",
        vmin=0,
        vmax=vmax
    )

    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Average declusterings per jet")

    ax.set_xlabel(r"$\ln\left(\frac{R}{\Delta R}\right)$", fontsize=14)
    ax.set_ylabel(r"$\ln\left(k_t/\mathrm{GeV}\right)$", fontsize=14)
    ax.set_title(title, fontsize=14)

    ax.grid(True, linestyle="--", alpha=0.45)

    if text_label is not None:
        ax.text(
            0.96,
            0.96,
            text_label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=13,
            color="white",
            fontweight="bold"
        )

    plt.tight_layout()

    if save:
        plt.savefig(filename, dpi=200, bbox_inches="tight")

    plt.show()


# In[ ]:


def make_primary_lund_image_from_xy(
    LundX,
    LundY,
    x_range=(0, 7),
    y_range=(-1, 7),
    bins=(50, 50),
):
    # Flatten event -> Lund point
    flat_x = ak.to_numpy(ak.flatten(LundX, axis=None))
    flat_y = ak.to_numpy(ak.flatten(LundY, axis=None))

    # Keep only finite values inside the plot range
    mask = (
        np.isfinite(flat_x) &
        np.isfinite(flat_y) &
        (flat_x >= x_range[0]) & (flat_x <= x_range[1]) &
        (flat_y >= y_range[0]) & (flat_y <= y_range[1])
    )

    flat_x = flat_x[mask]
    flat_y = flat_y[mask]

    H, x_edges, y_edges = np.histogram2d(
        flat_x,
        flat_y,
        bins=bins,
        range=[x_range, y_range]
    )

    return H, x_edges, y_edges


# In[ ]:





# **Lund Calculation with fast jet algorithm**
# <br>

# In[ ]:


#  Calculate primary Lund-plane points using FastJet.

vector.register_awkward()


def calculate_primary_lund_xy_fastjet(
    FJ,
    R=0.8,
    min_jet_pt=20.0,
    debug=False,
):

    all_lund_x = []
    all_lund_y = []

    n_jets = 0
    debug_printed = False

    jet_def = fastjet.JetDefinition(
        fastjet.cambridge_algorithm,
        R,
    )

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
            particles = ak.Array(
                [particle_records],
                with_name="Momentum4D",
            )

            cluster = fastjet.ClusterSequence(
                particles,
                jet_def,
            )

            # Reclustering the constituents into one exclusive C/A jet,
            # then obtaining its primary Lund declustering sequence.
            lund = cluster.exclusive_jets_lund_declusterings(
                njets=1
            )

            if "Delta" not in lund.fields or "kt" not in lund.fields:
                raise RuntimeError(
                    "FastJet Lund output does not contain the expected "
                    f"'Delta' and 'kt' fields. Found: {lund.fields}"
                )

            delta = ak.to_numpy(
                ak.flatten(lund["Delta"], axis=None)
            )

            kt = ak.to_numpy(
                ak.flatten(lund["kt"], axis=None)
            )

            good = (
                np.isfinite(delta)
                & np.isfinite(kt)
                & (delta > 0.0)
                & (kt > 0.0)
            )

            delta = delta[good]
            kt = kt[good]

            x = np.log(R / delta)
            y = np.log(kt)

            event_lund_x.extend(x.tolist())
            event_lund_y.extend(y.tolist())

        all_lund_x.append(event_lund_x)
        all_lund_y.append(event_lund_y)

    LundX = to_f32_jagged(
        ak.Array(all_lund_x)
    )

    LundY = to_f32_jagged(
        ak.Array(all_lund_y)
    )

    return LundX, LundY, n_jets


# In[ ]:


for isample, sample in enumerate(samples):
    events = sample["events"]
    FJ = events.FatJet

    # Calculate the Lund variables once using FastJet.
    FatJet_LundX, FatJet_LundY, n_jets = (
        calculate_primary_lund_xy_fastjet(
            FJ,
            R=0.8,
            min_jet_pt=20.0,
        )
    )

    # Store under the same names expected later by the friend-tree code.
    events["LundX"] = FatJet_LundX
    events["LundY"] = FatJet_LundY
    sample["events"] = events

    # Build the histogram from the already-computed Lund arrays.
    H, x_edges, y_edges = make_primary_lund_image_from_xy(
        FatJet_LundX,
        FatJet_LundY,
        x_range=(0, 7),
        y_range=(-1, 7),
        bins=(50, 50),
    )

    # Normalize by the number of jets actually processed by FastJet.
    H = H / max(n_jets, 1)

    n_lund_points = int(
        ak.sum(ak.num(FatJet_LundX, axis=1))
    )

    print(
        f'{sample["name"]}: '
        f'used {n_jets} jets and produced '
        f'{n_lund_points} Lund points'
    )

    plot_primary_lund_image(
        H,
        x_edges,
        y_edges,
        title="Primary Lund Plane",

        text_label=r"$r_{\mathrm{inv}} = 30\%$",

        save=False,
        filename=(
            f'{sample["name"]}_primary_lund_plane_fastjet.png'
        ),
        vmax=0.016,
    );


# In[ ]:


print()
print()

end = time.perf_counter()
total_seconds = end - start
print(f"Finished Lund Plane calculation and plotting in {total_seconds:.3f} seconds  ~{total_seconds/60:.3f} minutes")

print()
print()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# Helper functions to create friend tree


# Force size_of_FatJets == size_observables
def force_fatjet_layout(arr, FJ, name):


    # Clean missing values and bad numerical values first
    arr = ak.fill_none(arr, [], axis=0)
    arr = ak.fill_none(arr, 0.0, axis=-1)

    arr = ak.values_astype(arr, np.float64)
    arr = ak.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    arr = ak.values_astype(arr, np.float32)

    # FatJet counts per event. Missing FatJet collections count as 0.
    n_fatjet = ak.fill_none(ak.num(FJ, axis=1), 0)
    n_fatjet = ak.to_numpy(ak.values_astype(n_fatjet, np.int64))

    # Convert observable to Python
    arr_list = ak.to_list(arr)

    fixed = []

    for iev, njet in enumerate(n_fatjet):
        vals = arr_list[iev]

        if vals is None:
            vals = []

        # Make sure vals is a normal list
        vals = list(vals)

        # If there are too many values, keep only the real FatJet entries
        if len(vals) > njet:
            vals = vals[:njet]

        # If there are too few values, pad with zeros
        if len(vals) < njet:
            vals = vals + [0.0] * (njet - len(vals))

        fixed.append(vals)

    fixed = ak.Array(fixed)
    fixed = ak.values_astype(fixed, np.float32)

    return ak.to_packed(fixed)




# In[ ]:


# Create friend tree with new variables created

import os, re, uproot

BASE_MODELS = os.environ.get(
    "MB_MODELS_BASE",
    os.path.join(os.getcwd(), "models")
)

print("Friend trees saved under:", BASE_MODELS)


# --- Main loop ----

for sample in samples:
    events = sample["events"]
    FJ = events.FatJet


    # --- N2, N3 ---
    if {"N2", "N3"}.issubset(set(events.fields)):
        FatJet_N2 = to_f32_jagged(events["N2"])
        FatJet_N3 = to_f32_jagged(events["N3"])
    else:
        FatJet_N2 = to_f32_jagged(calculate_n2_loop(FJ, beta=1.0))
        FatJet_N3 = to_f32_jagged(calculate_n3_loop(FJ, beta=1.0))


    # --- Primary Lund plane points ---
    # Shape: event -> Lund point
    if {"LundX", "LundY"}.issubset(set(events.fields)):
        FatJet_LundX = to_f32_jagged(events["LundX"])
        FatJet_LundY = to_f32_jagged(events["LundY"])
    else:
        raise RuntimeError(
            "LundX/LundY are not stored in events. Run the Lund calculation/plotting cell first."
        )

    # --- event-level alignment guard ---
    n_evt = len(events)

    jets_like = {
        "N2":         FatJet_N2,
        "N3":         FatJet_N3,
    }

    
    # Make sure friend tree branch has the same number of events as the
    # original tree
    for name, arr in jets_like.items():
        if len(arr) != n_evt:
            raise RuntimeError(
                f"{name} has {len(arr)} events but Delphes has {n_evt} — "
                "don’t drop/reorder events when computing friend branches."
            )

  
            
    # Force all jet-level observables to use the exact FatJet layout
    jets_like = {
        name: force_fatjet_layout(arr, FJ, name)
        for name, arr in jets_like.items()
    }

    # Unique path construction 
    model_dir = os.path.join(BASE_MODELS, sample["model"])
    
    # These lines check that the file "events.root" exist in the directory to basically
    # make sure that the friend tree created below is created inside the same directory
    events_path = os.path.join(model_dir, "events.root")
    if not os.path.isfile(events_path):
        raise FileNotFoundError(f"Missing {events_path}")
     

    cluster = os.environ.get("ClusterId")
    proc = os.environ.get("ProcId")

    
    # Check if the code finds a local model directory or it is running in the cluster
    if not cluster or not proc:
        jobtag = os.environ.get("MB_JOBTAG", "jlocal.0")
        m = re.match(r"^j(\d+)\.(\d+)$", jobtag)

        if m:
            cluster, proc = m.group(1), m.group(2)
        else:
            cluster, proc = "local", "0"

    # Get the string for the model used and make sure the characters are safe to save,
    # it substitutes any "unsafe" charachters with a "-"
    model_tag = os.path.basename(model_dir)
    model_tag = re.sub(r"[^A-Za-z0-9._-]+", "-", model_tag)[:120]

    out_basename = f"events_friend_{model_tag}_j{cluster}.{proc}.root"
    friend_path = os.path.join(model_dir, out_basename)

    
    
    
    # The zipping in this next section is to avoid having a "nFatJet..." counter for each
    # variable, so there is only a counter per group structure (event->jet, and event->Lund point)
    
    # Zip jet-level variables into one FatJet record
    FatJet_record = ak.zip(
        {
            "N2": jets_like["N2"],
            "N3": jets_like["N3"],
        },
        depth_limit=2,
    )

    # Zip Lund variables into one Lund record 
    Lund_record = ak.zip(
        {
            "X": FatJet_LundX,
            "Y": FatJet_LundY,
        },
        depth_limit=2,
    )

    
    
    
    with uproot.recreate(friend_path) as fout:
        fout.mktree(
            "DelphesFriend",
            {
                "FatJet": ak.type(FatJet_record),
                "FatJet_LundX": ak.type(FatJet_LundX),
                "FatJet_LundY": ak.type(FatJet_LundY),
            }
        )

        fout["DelphesFriend"].extend(
            {
                "FatJet": FatJet_record,
                "FatJet_LundX": FatJet_LundX,
                "FatJet_LundY": FatJet_LundY,
            }
        )

    print("Wrote friend with zipped FatJet and Lund records:", friend_path)

    



# In[ ]:





# In[ ]:




