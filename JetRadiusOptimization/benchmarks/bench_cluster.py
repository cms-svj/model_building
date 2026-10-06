"""Benchmark: per-event Python PseudoJet loop vs awkward batched clustering.

Reproduces core.recluster_event and compares against the vectorized
fastjet._pyjet / awkward interface that clusters all events in one C++ call.
"""
import sys
import time
import numpy as np
import awkward as ak
import uproot
import fastjet as fj
import vector
vector.register_awkward()

if len(sys.argv) != 2:
    sys.exit("usage: python3 JetRadiusOptimization/benchmarks/bench_cluster.py <events.root>")
P = sys.argv[1]
t = uproot.open(P, handler=uproot.source.file.MemmapSource)['Delphes']

# GenCandidate four-vectors, exactly the pool core.py clusters for GenFatJet
br = t.arrays(['GenCandidate.PT', 'GenCandidate.Eta', 'GenCandidate.Phi',
               'GenCandidate.Mass'], entry_stop=200)
pt = br['GenCandidate.PT']
eta = br['GenCandidate.Eta']
phi = br['GenCandidate.Phi']
mass = br['GenCandidate.Mass']

px = pt * np.cos(phi)
py = pt * np.sin(phi)
pz = pt * np.sinh(eta)
E = np.sqrt(px**2 + py**2 + pz**2 + mass**2)

N_EVENTS = 60
RADII = [round(0.2 + 0.1 * i, 2) for i in range(15)]


def recluster_event(px_e, py_e, pz_e, e_e, radius, pt_min=15.0):
    """core.py's approach: one Python PseudoJet per particle, per event."""
    pjs = []
    for idx, (a_, b_, c_, d_) in enumerate(zip(px_e, py_e, pz_e, e_e)):
        pj = fj.PseudoJet(float(a_), float(b_), float(c_), float(d_))
        pj.set_user_index(idx)
        pjs.append(pj)
    if not pjs:
        return []
    jd = fj.JetDefinition(fj.antikt_algorithm, float(radius))
    cs = fj.ClusterSequence(pjs, jd)
    out = []
    for jet in fj.sorted_by_pt(cs.inclusive_jets(float(pt_min))):
        idxs = tuple(sorted(c.user_index() for c in jet.constituents()))
        out.append((float(jet.pt()), float(jet.eta()), float(jet.phi_std()),
                    float(jet.m()), idxs))
    return out


# ---- build the awkward record once (shared across all radii) ----
evts = ak.zip({
    "px": px[:N_EVENTS], "py": py[:N_EVENTS],
    "pz": pz[:N_EVENTS], "E": E[:N_EVENTS],
}, with_name="Momentum4D")

print(f"events={N_EVENTS}  radii={len(RADII)}  "
      f"mean constituents/event={float(np.mean(ak.num(pt[:N_EVENTS]))):.0f}")
print()

# ---------- A: current approach ----------
s = time.perf_counter()
ref = {}
for r in RADII:
    for ev in range(N_EVENTS):
        ref[(r, ev)] = recluster_event(px[ev], py[ev], pz[ev], E[ev], r)
t_loop = time.perf_counter() - s
print(f"A  per-event PseudoJet loop : {t_loop:7.2f} s "
      f"({1e3*t_loop/N_EVENTS/len(RADII):6.2f} ms/event/radius)")

# ---------- B: awkward batched, all events per radius in one C++ call ----------
s = time.perf_counter()
got = {}
for r in RADII:
    jd = fj.JetDefinition(fj.antikt_algorithm, float(r))
    cs = fj.ClusterSequence(evts, jd)
    jets = cs.inclusive_jets(min_pt=15.0)
    consts = cs.constituent_index(min_pt=15.0)
    jets_l = ak.to_list(jets)
    consts_l = ak.to_list(consts)
    for ev in range(N_EVENTS):
        got[(r, ev)] = (jets_l[ev], consts_l[ev])
t_awk = time.perf_counter() - s
print(f"B  awkward batched          : {t_awk:7.2f} s "
      f"({1e3*t_awk/N_EVENTS/len(RADII):6.2f} ms/event/radius)")
print(f"   speedup                  : x{t_loop/t_awk:.1f}")
print()

# ---------- correctness: identical constituent partitions ----------
mismatch = 0
checked = 0
for r in RADII:
    for ev in range(N_EVENTS):
        a_sets = sorted(tuple(sorted(j[4])) for j in ref[(r, ev)])
        b_sets = sorted(tuple(sorted(int(i) for i in c)) for c in got[(r, ev)][1])
        checked += 1
        if a_sets != b_sets:
            mismatch += 1
print(f"constituent-partition checks: {checked}, mismatches: {mismatch}")
print("CLUSTERING IDENTICAL" if mismatch == 0 else ">>> partitions differ <<<")
