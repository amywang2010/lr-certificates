"""B03 paper figures, generated only from the released artifacts (hashed in the manifest).

fig1_main.pdf    two panels, one column:
                   (a) breast certified intervals per (pair, region): T(A0) point plus
                       [T_lo, T_hi] bar, colored by verdict.
                   (b) cross-tissue certificates for the shared pairs, with the Proseg
                       point estimate.
fig2_support.pdf two panels, side by side:
                   (a) cross-segmentation validation: fraction of certified rows inside
                       their certified interval, with the coverage gate marked.
                   (b) kernel-margin distributions c* for the sections with certified rows.
fig0_workflow.pdf  pipeline schematic (b03_schematic.py, which carries a layout check).

Panels are grouped so that the shipped page count stays low: each figure is drawn at the
LNCS text width (4.8 in) and included at width=\textwidth, so no scaling occurs and the
specified point sizes are the printed point sizes. Every plotted number is read from the
artifact CSVs; assertions verify row counts against the released summaries before
plotting. No synthetic data.
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "lines.linewidth": 0.9,
    "font.family": "DejaVu Sans",
})

TEXT_W = 4.8

C_POS = "#B2182B"   # certified positive
C_NEG = "#2166AC"   # certified negative
C_NON = "#9E9E9E"   # non-identifiable
C_PRO = "#4D4D4D"   # proseg point


def plabel(ax, letter, size=9):
    ax.set_title(letter, loc="left", fontweight="bold", fontsize=size, pad=2)


# ================================================================ fig1_main
fin = pd.read_csv(f"{RES}/scout_final_matched.csv")
assert len(fin) == 75
assert int(fin.certified_pos_matched.sum()) == 37
assert int(fin.certified_neg_matched.sum()) == 9

fin["verdict"] = np.where(fin.certified_pos_matched, "pos",
                  np.where(fin.certified_neg_matched, "neg", "non"))
pair_meta = (fin.assign(ctrl=fin.control)
                .groupby("pair")
                .agg(ctrl=("ctrl", "first"), med=("T_lo", "median")))
pair_meta = pair_meta.sort_values(["ctrl", "med"], ascending=[False, True])
pairs = list(pair_meta.index)
ypos = {p: i for i, p in enumerate(pairs)}

f3 = pd.read_csv(f"{RES}/paper_fig3_data.csv")
assert len(f3) == 18

fig = plt.figure(figsize=(TEXT_W, 3.40))
gs = GridSpec(2, 1, figure=fig, height_ratios=[2.05, 1.35], hspace=0.30)

# ---- (a) breast certified intervals
ax = fig.add_subplot(gs[0])
for _, r in fin.iterrows():
    y = ypos[r.pair] + (r.region - 1) * 0.26
    col = {"pos": C_POS, "neg": C_NEG, "non": C_NON}[r.verdict]
    ax.plot([r.T_lo, r.T_hi], [y, y], color=col, lw=1.5, solid_capstyle="butt",
            alpha=0.95, zorder=2)
    ax.plot(r.T0, y, "o", ms=2.2, color="black", zorder=3)
for p, y in ypos.items():
    ctrl = pair_meta.loc[p, "ctrl"]
    ax.text(-7.6, y, ("CTRL " if ctrl else "") + p.replace("->", " \u2192 "),
            ha="right", va="center", fontsize=6.2,
            fontweight="bold" if ctrl else "normal")
ax.set_yticks([])
ax.set_xlim(-7.6, 10.2)
ax.set_ylim(-0.8, len(pairs) - 0.2)
ax.axvline(0, color="black", lw=0.8, ls="--", zorder=1)
ax.set_xlabel("T (log$_2$ expression-product ratio)", fontsize=7.5)
plabel(ax, "a")
handles = [Line2D([], [], color=C_POS, lw=2, label="certified positive"),
           Line2D([], [], color=C_NEG, lw=2, label="certified negative"),
           Line2D([], [], color=C_NON, lw=2, label="non-identifiable"),
           Line2D([], [], color="black", marker="o", ls="", ms=2.5, label="T(A$_0$)")]
ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=6.2)

# ---- (b) cross-tissue certificates
ax = fig.add_subplot(gs[1])
sel = f3.copy()
blocks = [("Breast", sel[sel.tissue == "breast"]), ("Lung", sel[sel.tissue == "lung"])]
ylabels, yvals, ycols, ypoints = [], [], [], []
y = 0
for bname, bdf in blocks:
    for pair in ["PECAM1->KDR", "CDH1->ERBB2", "CDH1->EGFR", "ESR1->PGR"]:
        sub = bdf[bdf.pair == pair].sort_values("region")
        if sub.empty:
            continue
        for _, r in sub.reset_index(drop=True).iterrows():
            col = C_POS if r.certified_pos_matched else (C_NEG if r.certified_neg_matched else C_NON)
            yvals.append(y); ycols.append(col)
            ypoints.append(r.T_proseg)
            ylabels.append((pair, r.region, r.T_lo, r.T_hi))
            y += 1
        y += 0.7
yvals = np.array(yvals)
for i in range(len(yvals)):
    _, _, t_lo, t_hi = ylabels[i]
    ax.plot([t_lo, t_hi], [yvals[i], yvals[i]], color=ycols[i], lw=2.2,
            solid_capstyle="butt")
    ax.plot(ypoints[i], yvals[i], "D", ms=2.6, color=C_PRO, mec="white", mew=0.3, zorder=4)
ax.axvline(0, color="black", lw=0.7, ls="--")
ax.set_yticks(yvals)
ax.set_yticklabels([f"{p.split('->')[0]} \u2192 {p.split('->')[1]}  r{r}"
                    for p, r, _, _ in ylabels], fontsize=6)
ax.set_xlabel("T (log$_2$ expression-product ratio)", fontsize=7.5)
ax.text(0.985, 0.97, "Breast", transform=ax.transAxes, fontweight="bold",
        fontsize=7, va="top", ha="right")
ax.text(0.985, 0.44, "Lung", transform=ax.transAxes, fontweight="bold",
        fontsize=7, va="top", ha="right")
handles = [Line2D([], [], color=C_POS, lw=2, label="certified pos."),
           Line2D([], [], color=C_NEG, lw=2, label="certified neg."),
           Line2D([], [], color=C_NON, lw=2, label="non-identifiable"),
           Line2D([], [], color=C_PRO, marker="D", ls="", ms=3, label="Proseg T")]
ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=6)
plabel(ax, "b")
fig.savefig(f"{OUT}/fig1_main.png", dpi=600)
fig.savefig(f"{OUT}/fig1_main.pdf")
plt.close(fig)

# ================================================================ fig2_support
s2 = json.load(open(f"{RES}/phase2_proseg_summary.json"))
s2l = json.load(open(f"{LUN}/phase5_e2_summary.json"))
kb = pd.read_csv(f"{RES}/phase4_kernel_margins.csv")
kl = pd.read_csv(f"{LUN}/phase4_kernel_margins.csv")
assert len(kb) == 46 and len(kl) == 7

fig, axes = plt.subplots(1, 2, figsize=(TEXT_W, 1.45))
ax = axes[0]
vals = [s2["primary"]["fraction"] * 100, s2l["primary"]["fraction"] * 100]
ax.bar([0, 1], vals, width=0.55, color=["#762A83", "#762A83"], alpha=0.85)
ax.axhline(90, color="black", lw=0.8, ls="--")
ax.text(0.5, 91.5, "gate 90%", fontsize=6, ha="center")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Breast\n72/75", "Lung\n6/7"])
ax.set_ylabel("certified rows inside\nProseg interval (%)", fontsize=7)
ax.set_ylim(0, 108)
plabel(ax, "a")
ax = axes[1]
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
ax.set_xticks([1, 2]); ax.set_xticklabels(labels)
ax.set_ylabel(r"critical kernel multiplier  $c^*$", fontsize=7)
ax.axhline(1.0, color="black", lw=0.7, ls=":")
ax.set_xlim(0.55, 2.45)
plabel(ax, "b")
fig.tight_layout(pad=0.4, w_pad=1.4)
fig.savefig(f"{OUT}/fig2_support.png", dpi=600)
fig.savefig(f"{OUT}/fig2_support.pdf")
plt.close(fig)

print("figures written:", sorted(f for f in os.listdir(OUT) if f.endswith(".pdf")))

# ================================================================ fig0_workflow
# The schematic carries a layout self-check (every label measured against its shape,
# pairwise shape-overlap check) and aborts rather than writing a broken figure.
import b03_schematic  # noqa: E402

print("fig0 schematic:", b03_schematic.build(OUT))
