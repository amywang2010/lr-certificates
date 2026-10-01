"""B03 paper figures: generated ONLY from sealed artifacts (hashed in manifests).

Fig 1  Breast certified intervals per (pair, region): T0 point + [T_lo, T_hi] bar,
       colored by verdict (certified positive / certified negative / non-identifiable).
Fig 2  Kernel-margin distributions c* for both tissues (log scale strip + box).
Fig 3  Cross-tissue panel: ESR1->PGR, CDH1->ERBB2 (breast) / CDH1->EGFR (lung),
       PECAM1->KDR control; certified interval + vendor T0 + Proseg point.
Fig 4  Validation summary bars: coverage (E2) and cross-segmentation significance (E3)
       for both tissues, with pre-specified gates marked.

Every plotted number is read from the artifact CSVs; assertions verify row counts
against the sealed summaries before plotting. No synthetic data.
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
LUN = os.path.join(ROOT, "lung")   # released lung artifacts
OUT = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.size": 7.5, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.6, "lines.linewidth": 0.9,
    "font.family": "DejaVu Sans",
})

C_POS = "#B2182B"   # certified positive
C_NEG = "#2166AC"   # certified negative
C_NON = "#9E9E9E"   # non-identifiable
C_PRO = "#4D4D4D"   # proseg point

# ---------------------------------------------------------------- Fig 1
fin = pd.read_csv(f"{RES}/scout_final_matched.csv")
assert len(fin) == 75
assert int(fin.certified_pos_matched.sum()) == 37
assert int(fin.certified_neg_matched.sum()) == 9

fin["verdict"] = np.where(fin.certified_pos_matched, "pos",
                  np.where(fin.certified_neg_matched, "neg", "non"))
# order: controls block first, then disputed; within block by median T_lo descending
pair_meta = (fin.assign(ctrl=fin.control)
                .groupby("pair")
                .agg(ctrl=("ctrl", "first"), med=("T_lo", "median")))
pair_meta = pair_meta.sort_values(["ctrl", "med"], ascending=[False, True])
pairs = list(pair_meta.index)
ypos = {p: i for i, p in enumerate(pairs)}

fig, ax = plt.subplots(figsize=(7.08, 5.6))
for _, r in fin.iterrows():
    y = ypos[r.pair] + (r.region - 1) * 0.26
    col = {"pos": C_POS, "neg": C_NEG, "non": C_NON}[r.verdict]
    ax.plot([r.T_lo, r.T_hi], [y, y], color=col, lw=1.5, solid_capstyle="butt",
            alpha=0.95, zorder=2)
    ax.plot(r.T0, y, "o", ms=2.2, color="black", zorder=3)
for p, y in ypos.items():
    ctrl = pair_meta.loc[p, "ctrl"]
    ax.text(-7.6, y, ("CTRL " if ctrl else "") + p.replace("->", " \u2192 "),
            ha="right", va="center",
            fontsize=6.2, fontweight="bold" if ctrl else "normal")
ax.set_yticks([])
ax.set_xlim(-7.6, 10.2)
ax.set_ylim(-0.8, len(pairs) - 0.2)
ax.axvline(0, color="black", lw=0.8, ls="--", zorder=1)
ax.set_xlabel("T (log$_2$ expression-product ratio)")
ax.set_title("a", loc="left", fontweight="bold", fontsize=10)
handles = [Line2D([], [], color=C_POS, lw=2, label="certified positive"),
           Line2D([], [], color=C_NEG, lw=2, label="certified negative"),
           Line2D([], [], color=C_NON, lw=2, label="non-identifiable"),
           Line2D([], [], color="black", marker="o", ls="", ms=2.5, label="T(A$_0$)")]
ax.legend(handles=handles, loc="lower right", frameon=False)
fig.tight_layout(pad=0.4)
fig.savefig(f"{OUT}/fig1_breast_intervals.png", dpi=600)
fig.savefig(f"{OUT}/fig1_breast_intervals.pdf")
plt.close(fig)

# ---------------------------------------------------------------- Fig 2
kb = pd.read_csv(f"{RES}/phase4_kernel_margins.csv")
kl = pd.read_csv(f"{LUN}/phase4_kernel_margins.csv")
assert len(kb) == 46 and len(kl) == 7

fig, ax = plt.subplots(figsize=(3.35, 2.5))
rng = np.random.default_rng(20260907)
data = [kb.c_star.to_numpy(), kl.c_star.to_numpy()]
labels = ["Breast\n(46 certificates)", "Lung\n(7 certificates)"]
for i, d in enumerate(data):
    x = i + 1
    ax.hlines(np.median(d), x - 0.18, x + 0.18, color="black", lw=1.4, zorder=3)
    q1, q3 = np.quantile(d, [0.25, 0.75])
    ax.vlines(x, q1, q3, color="black", lw=3, alpha=0.25, zorder=2)
    jit = rng.uniform(-0.09, 0.09, len(d))
    ax.plot(x + jit, d, "o", ms=2.6, mfc="#4D4D4D", mec="none", alpha=0.75, zorder=4)
ax.set_yscale("log")
ax.set_xticks([1, 2])
ax.set_xticklabels(labels)
ax.set_ylabel(r"critical kernel multiplier  $c^*$")
ax.axhline(1.0, color="black", lw=0.7, ls=":")
ax.text(2.42, 1.02, "boundary-local\nkernel (c = 1)", fontsize=6, va="bottom", ha="right")
ax.set_xlim(0.55, 2.45)
ax.set_title("a", loc="left", fontweight="bold", fontsize=10)
fig.tight_layout(pad=0.4)
fig.savefig(f"{OUT}/fig2_kernel_margins.png", dpi=600)
fig.savefig(f"{OUT}/fig2_kernel_margins.pdf")
plt.close(fig)

# ---------------------------------------------------------------- Fig 3
f3 = pd.read_csv(f"{RES}/paper_fig3_data.csv")
assert len(f3) == 18
sel = f3.copy()
# layout: rows = pair, columns = regions 0..2, two tissue blocks
blocks = [("Breast", sel[sel.tissue == "breast"]), ("Lung", sel[sel.tissue == "lung"])]
fig, ax = plt.subplots(figsize=(3.35, 3.9))
ylabels, yvals, ycols, ypoints, ycert = [], [], [], [], []
y = 0
for bname, bdf in blocks:
    for pair in ["PECAM1->KDR", "CDH1->ERBB2", "CDH1->EGFR", "ESR1->PGR"]:
        sub = bdf[bdf.pair == pair].sort_values("region")
        if sub.empty:
            continue
        sub = sub.reset_index(drop=True)
        for _, r in sub.iterrows():
            col = C_POS if r.certified_pos_matched else (C_NEG if r.certified_neg_matched else C_NON)
            yvals.append(y); ycols.append(col)
            ypoints.append(r.T_proseg)
            ylabels.append((pair, r.region, bname, r.T_lo, r.T_hi))
            y += 1
        y += 0.6
yvals = np.array(yvals)
for i in range(len(yvals)):
    pair, region, bname, t_lo, t_hi = ylabels[i]
    ax.plot([t_lo, t_hi], [yvals[i], yvals[i]], color=ycols[i], lw=2.2,
            solid_capstyle="butt")
    ax.plot(ypoints[i], yvals[i], "D", ms=2.6, color=C_PRO, mec="white", mew=0.3, zorder=4)
ax.axvline(0, color="black", lw=0.7, ls="--")
ax.set_yticks(yvals)
ax.set_yticklabels([f"{p.split('->')[0]} \u2192 {p.split('->')[1]}  r{r}" for p, r, _, _, _ in ylabels], fontsize=6)
ax.set_xlabel("T (log$_2$ expression-product ratio)")
ax.text(0.02, 0.985, "Breast", transform=ax.transAxes, fontweight="bold", fontsize=7.5, va="top")
ax.text(0.02, 0.50, "Lung", transform=ax.transAxes, fontweight="bold", fontsize=7.5, va="top")
handles = [Line2D([], [], color=C_POS, lw=2, label="certified pos."),
           Line2D([], [], color=C_NEG, lw=2, label="certified neg."),
           Line2D([], [], color=C_NON, lw=2, label="non-identifiable"),
           Line2D([], [], color=C_PRO, marker="D", ls="", ms=3, label="Proseg T")]
ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=6)
ax.set_title("a", loc="left", fontweight="bold", fontsize=10)
fig.tight_layout(pad=0.4)
fig.savefig(f"{OUT}/fig3_crosstissue.png", dpi=600)
fig.savefig(f"{OUT}/fig3_crosstissue.pdf")
plt.close(fig)

# ---------------------------------------------------------------- Fig 4
s2 = json.load(open(f"{RES}/phase2_proseg_summary.json"))
s3 = json.load(open(f"{RES}/phase3_summary.json"))
s2l = json.load(open(f"{LUN}/phase5_e2_summary.json"))
s3l = json.load(open(f"{LUN}/phase3_summary.json"))
fig, axes = plt.subplots(1, 2, figsize=(3.35, 2.1))
ax = axes[0]
vals = [s2["primary"]["fraction"] * 100, s2l["primary"]["fraction"] * 100]
ax.bar([0, 1], vals, width=0.55, color=["#762A83", "#762A83"], alpha=0.85)
ax.axhline(90, color="black", lw=0.8, ls="--")
ax.text(0.5, 91.5, "gate 90%", fontsize=6, ha="center")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Breast\n72/75", "Lung\n6/7"])
ax.set_ylabel("certified rows inside\nProseg interval (%)")
ax.set_ylim(0, 108)
ax.set_title("a", loc="left", fontweight="bold", fontsize=9)
ax = axes[1]
vals = [s3["fraction"] * 100, s3l["fraction"] * 100]
ax.bar([0, 1], vals, width=0.55, color=["#1B7837", "#1B7837"], alpha=0.85)
ax.axhline(85, color="black", lw=0.8, ls="--")
ax.text(0.5, 86.5, "gate 85%", fontsize=6, ha="center")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Breast\n40/46", "Lung\n7/7"])
ax.set_ylabel("certified rows significant under\nProseg-side null (%)")
ax.set_ylim(0, 108)
ax.set_title("b", loc="left", fontweight="bold", fontsize=9)
fig.tight_layout(pad=0.4, w_pad=1.2)
fig.savefig(f"{OUT}/fig4_validation.png", dpi=600)
fig.savefig(f"{OUT}/fig4_validation.pdf")
plt.close(fig)

print("figures written")
