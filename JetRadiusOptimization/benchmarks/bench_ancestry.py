"""Corrected benchmark: bitmask ancestry iterated to a true fixed point."""
import sys
import time
import numpy as np
import uproot

if len(sys.argv) != 2:
    sys.exit("usage: python3 JetRadiusOptimization/benchmarks/bench_ancestry.py <events.root>")
P = sys.argv[1]
t = uproot.open(P, handler=uproot.source.file.MemmapSource)['Delphes']
a = t.arrays(['GenParticle.fUniqueID', 'GenParticle.M1', 'GenParticle.M2',
              'GenCandidate.fUniqueID', 'DarkHadronCandidate.fUniqueID'], entry_stop=200)


def ancestor_targets(start_index, mother1, mother2, targets):
    """Verbatim from core.py."""
    found = set()
    queue = [int(start_index)]
    visited = set()
    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        visited.add(current)
        if current in targets:
            found.add(current)
            continue
        if current <= 1 or current >= len(mother1):
            continue
        for parent in (int(mother1[current]), int(mother2[current])):
            if 1 < parent < len(mother1) and parent not in visited and parent not in queue:
                queue.append(parent)
    return found


def ancestor_masks(mother1, mother2, target_indices, max_passes=200):
    """Vectorized: uint64 bitmask per GenParticle, iterated to a fixed point.

    Bit j of mask[i] is set iff particle i descends from target j, using the
    same 'stop at the first target' rule as the BFS (targets do not inherit
    their own ancestors' bits).
    """
    n = len(mother1)
    mask = np.zeros(n, dtype=np.uint64)
    is_target = np.zeros(n, dtype=bool)
    for bit, tgt in enumerate(target_indices):
        if 0 <= tgt < n:
            mask[tgt] |= np.uint64(1) << np.uint64(bit)
            is_target[tgt] = True
    if not len(target_indices):
        return mask, 0

    idx = np.arange(n)
    ok1 = (mother1 > 1) & (mother1 < n)
    ok2 = (mother2 > 1) & (mother2 < n)
    m1 = np.where(ok1, mother1, 0).astype(np.int64)
    m2 = np.where(ok2, mother2, 0).astype(np.int64)
    inherit = (~is_target) & (idx > 1)
    src1 = inherit & ok1
    src2 = inherit & ok2

    for npass in range(1, max_passes + 1):
        acc = mask.copy()
        acc[src1] |= mask[m1[src1]]
        acc[src2] |= mask[m2[src2]]
        if np.array_equal(acc, mask):
            return mask, npass
        mask = acc
    raise RuntimeError("ancestry bitmask did not converge")


N = 60
t_bfs = t_vec = 0.0
mismatch = 0
tot_cand = 0
passes = []

for ev in range(N):
    uid = np.asarray(a['GenParticle.fUniqueID'][ev], dtype=np.int64)
    m1 = np.asarray(a['GenParticle.M1'][ev], dtype=np.int64)
    m2 = np.asarray(a['GenParticle.M2'][ev], dtype=np.int64)
    gen_index = {int(u): i for i, u in enumerate(uid)}
    targets = sorted({gen_index[int(u)] for u in a['DarkHadronCandidate.fUniqueID'][ev]
                      if int(u) in gen_index})
    starts = [gen_index[int(u)] for u in a['GenCandidate.fUniqueID'][ev] if int(u) in gen_index]
    tot_cand += len(starts)
    tset = set(targets)

    s = time.perf_counter()
    ref = [ancestor_targets(c, m1, m2, tset) for c in starts]
    t_bfs += time.perf_counter() - s

    s = time.perf_counter()
    mask, np_ = ancestor_masks(m1, m2, targets)
    got = []
    for c in starts:
        mv = int(mask[c])
        got.append({targets[b] for b in range(len(targets)) if mv >> b & 1})
    t_vec += time.perf_counter() - s
    passes.append(np_)

    mismatch += sum(x != y for x, y in zip(ref, got))

print(f"events                     : {N}")
print(f"GenCandidates              : {tot_cand:,}")
print(f"fixed-point passes         : mean {np.mean(passes):.1f}, max {max(passes)}")
print()
print(f"current per-particle BFS   : {1e3*t_bfs/N:7.2f} ms/event")
print(f"vectorized bitmask         : {1e3*t_vec/N:7.2f} ms/event")
print(f"speedup                    : x{t_bfs/t_vec:.1f}")
print()
print(f"result mismatches          : {mismatch} / {tot_cand}")
print("PHYSICS IDENTICAL" if mismatch == 0 else ">>> NOT identical - investigate <<<")
