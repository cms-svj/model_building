"""Does clustering at R_max once, then reclustering only within each seed jet,
reproduce direct anti-kT clustering at every smaller R?

If yes, a 15-point radius scan costs 1 large + 15 small clusterings on much
smaller particle sets instead of 15 full-event clusterings.
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
    sys.exit("usage: python3 JetRadiusOptimization/benchmarks/bench_seeded.py <events.root>")
P = sys.argv[1]
t = uproot.open(P, handler=uproot.source.file.MemmapSource)['Delphes']
br = t.arrays(['GenCandidate.PT', 'GenCandidate.Eta', 'GenCandidate.Phi',
               'GenCandidate.Mass'], entry_stop=200)
pt, eta, phi, m = (br['GenCandidate.PT'], br['GenCandidate.Eta'],
                   br['GenCandidate.Phi'], br['GenCandidate.Mass'])
px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
E = np.sqrt(px**2 + py**2 + pz**2 + m**2)

N = 50
RADII = [round(0.2 + 0.1 * i, 2) for i in range(15)]  # 0.2 .. 1.6
R_MAX = 1.6
PT_MIN = 15.0


def direct(radius):
    evts = ak.zip({"px": px[:N], "py": py[:N], "pz": pz[:N], "E": E[:N]},
                  with_name="Momentum4D")
    cs = fj.ClusterSequence(evts, fj.JetDefinition(fj.antikt_algorithm, radius))
    return ak.to_list(cs.constituent_index(min_pt=PT_MIN))


# ---- direct reference for every radius ----
s = time.perf_counter()
ref = {r: direct(r) for r in RADII}
t_direct = time.perf_counter() - s

# ---- seeded: one R_MAX clustering with NO pt cut, then recluster inside ----
s = time.perf_counter()
evts = ak.zip({"px": px[:N], "py": py[:N], "pz": pz[:N], "E": E[:N]},
              with_name="Momentum4D")
seed_cs = fj.ClusterSequence(evts, fj.JetDefinition(fj.antikt_algorithm, R_MAX))
seed_idx = ak.to_list(seed_cs.constituent_index(min_pt=0.0))

seeded = {r: [] for r in RADII}
for ev in range(N):
    pxe, pye, pze, Ee = (np.asarray(px[ev]), np.asarray(py[ev]),
                         np.asarray(pz[ev]), np.asarray(E[ev]))
    for r in RADII:
        out = []
        for grp in seed_idx[ev]:
            g = np.asarray([int(i) for i in grp], dtype=np.int64)
            if len(g) == 0:
                continue
            sub = ak.zip({"px": [pxe[g]], "py": [pye[g]],
                          "pz": [pze[g]], "E": [Ee[g]]}, with_name="Momentum4D")
            cs = fj.ClusterSequence(sub, fj.JetDefinition(fj.antikt_algorithm, r))
            for jc in ak.to_list(cs.constituent_index(min_pt=PT_MIN))[0]:
                out.append(sorted(int(g[int(i)]) for i in jc))
        seeded[r].append(out)
t_seeded = time.perf_counter() - s

print(f"events={N}  radii={len(RADII)}  R_max={R_MAX}  pt_min={PT_MIN}")
print(f"direct  (15 full-event clusterings) : {t_direct:6.2f} s")
print(f"seeded  (1 seed + in-jet reclusters): {t_seeded:6.2f} s")
print()

bad = {}
for r in RADII:
    nmis = 0
    for ev in range(N):
        a = sorted(tuple(sorted(int(i) for i in c)) for c in ref[r][ev])
        b = sorted(tuple(x) for x in seeded[r][ev])
        if a != b:
            nmis += 1
    if nmis:
        bad[r] = nmis

if not bad:
    print("EXACT for every radius -> seeding is safe")
else:
    print("MISMATCHES (radius: n events differing out of %d)" % N)
    for r, n in sorted(bad.items()):
        print(f"   R={r}: {n}")
