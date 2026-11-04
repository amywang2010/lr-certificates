"""B03: exhaustive small-instance check of the flow ends and the attainment criteria.

The three sections exercise the construction on real data but cannot check it against the
definition of the uncertainty set on instances small enough to enumerate: the signature
contraction, the shared budgets and the joint-attainment criterion. This script does that. It
draws random tiny
instances (a handful of cells, a handful of movable transcripts, two genes with possibly
equal, overlapping or disjoint counted sets), enumerates every feasible assignment in the
declared set, and compares:

  * max_A N_g(A) - N0_g      against  the per-gene max flow rho_g
  * N0_g - min_A N_g(A)      against  the closed form D_g = sum_{s in C_g} min(a_g(s), L_s)
  * existence of one A attaining both per-gene maxima
                             against  rho_LR == rho_L + rho_R
  * existence of one A attaining both per-gene minima
                             against  D_LR == D_L + D_R

The assignment model is the one Definition 2 states, with nothing added: a movable transcript
stays, drops to the extracellular compartment, or moves to any cell of its donor set; departures
from a cell are capped by its loss budget and arrivals into a cell by its gain budget, and both
budgets are properties of the cell, so the two genes compete for them. Counted sets are per
gene. A transcript whose source already lies in its own counted set is excluded from the flow,
as Proposition 3 states, and the enumeration is what shows the exclusion is without loss.

The flows here are built from the instance description, independently of the release's
vectorised construction in b03_exact_end.py, so agreement is evidence about the construction
and not about a shared bug.

Usage:  python src/b03_enumerate_check.py [n_instances]
Writes: results/enumeration_check.json
"""
import itertools
import json
import os
import random
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXTRA = -1


class Dinic:
    """Minimal integer max flow, so the checker does not share code with the release."""

    def __init__(self, n):
        self.n = n
        self.to = []
        self.cap = []
        self.head = [[] for _ in range(n)]

    def add(self, u, v, c):
        self.head[u].append(len(self.to))
        self.to.append(v)
        self.cap.append(int(c))
        self.head[v].append(len(self.to))
        self.to.append(u)
        self.cap.append(0)

    def maxflow(self, s, t):
        flow = 0
        while True:
            level = [-1] * self.n
            level[s] = 0
            q = [s]
            while q:
                u = q.pop(0)
                for e in self.head[u]:
                    if self.cap[e] > 0 and level[self.to[e]] < 0:
                        level[self.to[e]] = level[u] + 1
                        q.append(self.to[e])
            if level[t] < 0:
                return flow
            it = [0] * self.n
            while True:
                pushed = self._dfs(s, t, 10 ** 9, level, it)
                if pushed == 0:
                    break
                flow += pushed

    def _dfs(self, u, t, f, level, it):
        if u == t:
            return f
        while it[u] < len(self.head[u]):
            e = self.head[u][it[u]]
            v = self.to[e]
            if self.cap[e] > 0 and level[v] == level[u] + 1:
                d = self._dfs(v, t, min(f, self.cap[e]), level, it)
                if d > 0:
                    self.cap[e] -= d
                    self.cap[e ^ 1] += d
                    return d
            it[u] += 1
        return 0


def enumerate_instance(inst):
    """All feasible assignments, as (arrivals, departures) per gene."""
    outcomes = set()
    tx = inst["tx"]
    opts = []
    for t in tx:
        o = ["stay", "extra"] + list(t["donors"])
        opts.append(o)
    for choice in itertools.product(*opts):
        loss = {}
        gain = {}
        ok = True
        for t, c in zip(tx, choice):
            if c != "stay" and t["src"] >= 0:
                loss[t["src"]] = loss.get(t["src"], 0) + 1
            if c not in ("stay", "extra"):
                gain[c] = gain.get(c, 0) + 1
        for s, k in loss.items():
            if k > inst["L_cap"][s]:
                ok = False
                break
        if ok:
            for d, k in gain.items():
                if k > inst["G_cap"][d]:
                    ok = False
                    break
        if not ok:
            continue
        arr = {"L": 0, "R": 0}
        dep = {"L": 0, "R": 0}
        for t, c in zip(tx, choice):
            C = inst["C"][t["gene"]]
            if c not in ("stay", "extra") and c in C:
                arr[t["gene"]] += 1
            if t["src"] >= 0 and t["src"] in C and c != "stay":
                dep[t["gene"]] += 1
        outcomes.add((arr["L"] - dep["L"], arr["R"] - dep["R"]))
    return outcomes


def flow_value(inst, genes):
    """Max total inflow on the signature-contracted network for the given genes."""
    sig = {}
    for t in inst["tx"]:
        if t["gene"] not in genes:
            continue
        C = inst["C"][t["gene"]]
        if t["src"] in C:
            continue                      # Proposition 3's exclusion, imposed, not assumed
        dest = tuple(sorted(d for d in t["donors"] if d in C))
        if not dest:
            continue
        sig[(t["gene"], t["src"], dest)] = sig.get((t["gene"], t["src"], dest), 0) + 1
    srcs = sorted({k[1] for k in sig})
    dests = sorted({d for k in sig for d in k[2]})
    S = 0
    idx_src = {s: 2 + i for i, s in enumerate(srcs)}
    off_grp = 2 + len(srcs)
    idx_grp = {k: off_grp + i for i, k in enumerate(sig)}
    off_dst = off_grp + len(sig)
    idx_dst = {d: off_dst + i for i, d in enumerate(dests)}
    n = off_dst + len(dests) + 1
    T = n - 1
    g = Dinic(n)
    big = sum(sig.values()) + 1
    for s in srcs:
        g.add(S, idx_src[s], inst["L_cap"][s] if s >= 0 else big)
    for k, cnt in sig.items():
        g.add(idx_src[k[1]], idx_grp[k], cnt)
        for d in k[2]:
            g.add(idx_grp[k], idx_dst[d], cnt)
    for d in dests:
        g.add(idx_dst[d], T, inst["G_cap"][d])
    return g.maxflow(S, T)


def closed_form_lower(inst, gene):
    C = inst["C"][gene]
    a = {}
    for t in inst["tx"]:
        if t["gene"] == gene and t["src"] in C:
            a[t["src"]] = a.get(t["src"], 0) + 1
    return sum(min(a.get(s, 0), inst["L_cap"][s]) for s in C)


def joint_lower(inst):
    cells = inst["C"]["L"] | inst["C"]["R"]
    tot = 0
    for s in cells:
        a = 0
        for t in inst["tx"]:
            if t["src"] == s and s in inst["C"][t["gene"]]:
                a += 1
        tot += min(a, inst["L_cap"][s])
    return tot


def draw(rng, n_cells, n_tx):
    cells = list(range(n_cells))
    mode = rng.choice(["same", "overlap", "disjoint"])
    if mode == "same":
        C = {g: set(rng.sample(cells, max(2, n_cells - 1))) for g in "LR"}
        C["R"] = set(C["L"])
    elif mode == "overlap":
        base = set(rng.sample(cells, max(2, n_cells - 1)))
        C = {"L": base, "R": set(rng.sample(sorted(base), max(1, len(base) - 1)))
             | set(rng.sample([c for c in cells if c not in base], 1))}
    else:
        half = max(1, n_cells // 2)
        left = rng.sample(cells, half)
        right = [c for c in cells if c not in left]
        C = {"L": set(left), "R": set(right)}
    L_cap = {c: rng.randint(0, 3) for c in cells}
    G_cap = {c: rng.randint(0, 3) for c in cells}
    tx = []
    for _ in range(n_tx):
        gene = rng.choice("LR")
        pool = [c for c in cells if c in C[gene]]
        if not pool:
            continue
        src = rng.choice(cells + [EXTRA, EXTRA])
        donors = rng.sample(pool, rng.randint(1, min(2, len(pool))))
        tx.append({"gene": gene, "src": src, "donors": donors})
    return {"C": C, "L_cap": L_cap, "G_cap": G_cap, "tx": tx, "mode": mode}


def main():
    n_inst = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    rng = random.Random(20261024)
    stats = dict(instances=0, empty=0, upper_attained=0, lower_attained=0,
                 upper_criterion_true=0, lower_criterion_true=0,
                 flow_mismatch=0, lower_mismatch=0, upper_criterion_wrong=0,
                 lower_criterion_wrong=0, same_class=0, disjoint=0)
    first_bad = []
    for _ in range(n_inst):
        inst = draw(rng, rng.randint(3, 5), rng.randint(2, 5))
        if not inst["tx"]:
            continue
        stats["instances"] += 1
        stats["same_class"] += int(inst["mode"] == "same")
        stats["disjoint"] += int(inst["mode"] == "disjoint")
        out = enumerate_instance(inst)
        if not out:
            stats["empty"] += 1
            continue
        dL = max(o[0] for o in out)
        dR = max(o[1] for o in out)
        lL = min(o[0] for o in out)
        lR = min(o[1] for o in out)
        up = (dL, dR) in out
        lo = (lL, lR) in out
        stats["upper_attained"] += int(up)
        stats["lower_attained"] += int(lo)

        rhoL = flow_value(inst, {"L"})
        rhoR = flow_value(inst, {"R"})
        rhoLR = flow_value(inst, {"L", "R"})
        if (rhoL, rhoR) != (dL, dR):
            stats["flow_mismatch"] += 1
            if len(first_bad) < 3:
                first_bad.append(("flow", inst, (rhoL, rhoR), (dL, dR)))
        DL, DR = closed_form_lower(inst, "L"), closed_form_lower(inst, "R")
        if (DL, DR) != (-lL, -lR):
            stats["lower_mismatch"] += 1
            if len(first_bad) < 3:
                first_bad.append(("lower", inst, (DL, DR), (-lL, -lR)))
        crit_up = (rhoLR == rhoL + rhoR)
        crit_lo = (joint_lower(inst) == DL + DR)
        stats["upper_criterion_true"] += int(crit_up)
        stats["lower_criterion_true"] += int(crit_lo)
        if crit_up != up:
            stats["upper_criterion_wrong"] += 1
            if len(first_bad) < 3:
                first_bad.append(("upper-criterion", inst, crit_up, up))
        if crit_lo != lo:
            stats["lower_criterion_wrong"] += 1
            if len(first_bad) < 3:
                first_bad.append(("lower-criterion", inst, crit_lo, lo))

    stats["ok"] = (stats["flow_mismatch"] == 0 and stats["lower_mismatch"] == 0
                   and stats["upper_criterion_wrong"] == 0
                   and stats["lower_criterion_wrong"] == 0)
    stats["seed"] = 20261024
    stats["max_cells"] = 5
    stats["max_movable_transcripts"] = 5
    with open(os.path.join(ROOT, "results", "enumeration_check.json"), "w",
              newline="\n") as fh:
        json.dump(stats, fh, indent=1, sort_keys=True)
    print(json.dumps(stats, indent=1, sort_keys=True))
    for b in first_bad:
        print("MISMATCH", b[0], b[1], b[2], b[3])
    return 0 if stats["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
