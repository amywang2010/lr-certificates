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

# Figures are drawn at the LNCS text width (4.8 in) and included at
# width=\textwidth, so no scaling occurs and the specified point sizes are the
# printed point sizes. Drawing wider and scaling down would shrink the type.
TEXT_W = 4.8

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

fig, ax = plt.subplots(figsize=(TEXT_W, 4.0))
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

fig, ax = plt.subplots(figsize=(TEXT_W, 1.95))

# Horizontal layout: one row per tissue, value on a log x-axis. Horizontal
# reading matches the reader's eye and leaves room to annotate each median,
# which a two-column strip on a log y-axis could not.
data = [("Breast", kb.c_star.to_numpy()), ("Lung", kl.c_star.to_numpy())]


def beeswarm(vals, width=0.34, sep=0.052):
    """Deterministic 1-D beeswarm offsets, so points never overplot."""
    order = np.argsort(vals, kind="stable")
    offs = np.zeros(len(vals))
    last, level = None, 0
    for k in order:
        level = 0 if last is None or vals[k] - last > sep else level + 1
        offs[k] = 0.0 if level == 0 else (width if level % 2 else -width) * (
            (level + 1) // 2) / max(1, (len(vals) // 8) + 1)
        last = vals[k]
    return np.clip(offs, -0.42, 0.42)


for row, (name, d) in enumerate(data):
    d = np.asarray(d, float)
    q1, med, q3 = np.quantile(d, [0.25, 0.5, 0.75])
    ax.plot([q1, q3], [row, row], color="#9E9E9E", lw=4.5,
            solid_capstyle="butt", zorder=2)
    ax.plot(d, row + beeswarm(d), "o", ms=3.0, mfc="white", mec="#4D4D4D",
            mew=0.6, zorder=3)
    ax.plot([med, med], [row - 0.30, row + 0.30], color="black", lw=1.6,
            zorder=4)
    ax.annotate(f"median {med:.2f}", (med, row + 0.34),
                ha="center", va="bottom", fontsize=6.2)

ax.axvline(1.0, color="black", lw=0.8, ls=":")
ax.annotate("constructed set, $c = 1$", (1.0, 1.44), ha="center", va="bottom",
            fontsize=6.2)
ax.set_xscale("log")
ax.set_yticks([0, 1])
ax.set_yticklabels([f"{n}\n({len(d)} certificates)" for n, d in data],
                   fontsize=6.6)
ax.set_ylim(-0.55, 1.62)
ax.set_xlabel("critical kernel multiplier  $c^*$  (log scale)", fontsize=6.8)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.grid(axis="x", which="major", color="#E8E8E8", lw=0.5, zorder=0)
fig.tight_layout(pad=0.3)
fig.savefig(f"{OUT}/fig2_kernel_margins.png", dpi=600)
fig.savefig(f"{OUT}/fig2_kernel_margins.pdf")
plt.close(fig)

# ---------------------------------------------------------------- Fig 3
f3 = pd.read_csv(f"{RES}/paper_fig3_data.csv")
assert len(f3) == 18
sel = f3.copy()
# layout: rows = pair, columns = regions 0..2, two tissue blocks
blocks = [("Breast", sel[sel.tissue == "breast"]), ("Lung", sel[sel.tissue == "lung"])]
fig, ax = plt.subplots(figsize=(TEXT_W, 3.4))
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
fig, axes = plt.subplots(1, 2, figsize=(TEXT_W, 1.9))
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

# ------------------------------------------------------- Fig 0: schematic
# Pipeline schematic in the style of a mermaid/flowchart diagram: white shapes
# with thin black outlines, shape chosen by role (parallelogram = data, box =
# process, hexagon = test, stadium = output), straight orthogonal connectors,
# and grouping bars over the columns. The point the figure has to make is that
# the observed assignment and every permutation enter the SAME interval
# machinery, which is what makes the null a matched null.
from matplotlib.patches import Polygon, FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(TEXT_W, 2.45))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

LW = 0.7
BLACK = "#1A1A1A"


def data_shape(x, y, w, h, text, skew=0.13):
    """Parallelogram: an input."""
    s = w * skew
    pts = [(x + s, y), (x + w, y), (x + w - s, y + h), (x, y + h)]
    ax.add_patch(Polygon(pts, closed=True, facecolor="white", edgecolor=BLACK,
                         linewidth=LW))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=5.9,
            linespacing=1.35)


def box_shape(x, y, w, h, text):
    """Rectangle: a process."""
    ax.add_patch(Polygon([(x, y), (x + w, y), (x + w, y + h), (x, y + h)],
                         closed=True, facecolor="white", edgecolor=BLACK,
                         linewidth=LW))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=5.9,
            linespacing=1.35)


def hex_shape(x, y, w, h, text):
    """Hexagon: the test applied to the interval."""
    k = w * 0.18
    pts = [(x + k, y), (x + w - k, y), (x + w, y + h / 2), (x + w - k, y + h),
           (x + k, y + h), (x, y + h / 2)]
    ax.add_patch(Polygon(pts, closed=True, facecolor="white", edgecolor=BLACK,
                         linewidth=LW))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=5.9,
            linespacing=1.35)


def stadium(x, y, w, h, text):
    """Stadium: the reported output."""
    ax.add_patch(FancyBboxPatch((x + h / 2, y), w - h, h,
                                boxstyle=f"round,pad=0,rounding_size={h / 2}",
                                facecolor="white", edgecolor=BLACK,
                                linewidth=LW, mutation_aspect=1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=5.9,
            linespacing=1.35)


def arrow(pts, head_at_end=True):
    """Orthogonal polyline connector with a filled head at the last segment."""
    for i in range(len(pts) - 1):
        (x1, y1), (x2, y2) = pts[i], pts[i + 1]
        last = i == len(pts) - 2
        ax.add_patch(FancyArrowPatch(
            (x1, y1), (x2, y2), arrowstyle="-|>" if last else "-",
            mutation_scale=6 if last else 1, linewidth=LW, color=BLACK,
            shrinkA=0, shrinkB=0, joinstyle="miter"))


# ---- grouping bars over the columns
def group(x0, x1, label, y=90.5):
    ax.plot([x0, x1], [y, y], color=BLACK, lw=0.6)
    ax.plot([x0 + 0.14 * (x1 - x0), x1 - 0.14 * (x1 - x0)], [y, y],
            color="white", lw=2.6)
    ax.text((x0 + x1) / 2, y + 1.8, label, ha="center", va="bottom",
            fontsize=6.6)


group(2, 45, "Uncertainty set")
group(50, 71, "Interval machinery")
group(76, 95, "Inference")

# ---- observed path (upper row)
data_shape(2, 62, 19, 15, "Molecules and\nvendor masks\n(assignment $A_0$)")
box_shape(26, 61, 19, 17,
          "Build $U$\nband $3\\,\\mu$m, donor $6\\,\\mu$m,\ncell caps $+50\\%/\\!-60\\%$")
box_shape(50, 40, 21, 32,
          "Per-gene interval\nmachinery\n\n$[\\min,\\max]$ counts,\nthen $T_{\\mathrm{lo}}$,\n"
          "$T_{\\mathrm{hi}}$ over $U$")
hex_shape(76, 46, 19, 16,
          "Interval sign\nand matched-null\nsign agree?")
stadium(76, 20, 19, 14, "Certificate\ncertified $+$, certified $-$,\nor non-identifiable")

# ---- permuted path (lower row)
data_shape(2, 33, 19, 15,
           "Region $\\times$ decile\nblocks frozen\nfrom $A_0$")
box_shape(26, 32, 19, 15,
          "Matched null\nlabels permuted\nwithin blocks, $B=1000$")

arrow([(21, 69.5), (26, 69.5)])
arrow([(21, 40.5), (23.2, 40.5), (23.2, 39.5), (26, 39.5)])
arrow([(45, 69.5), (47.6, 69.5), (47.6, 66), (50, 66)])
arrow([(45, 39.5), (47.6, 39.5), (47.6, 46), (50, 46)])
arrow([(71, 54), (76, 54)])
arrow([(85.5, 46), (85.5, 34)])

ax.text(1, 4.5,
        "Both inputs enter the same interval machinery, which is what makes the null "
        "matched to the certified object.\nRegions, crop windows and gates are "
        "frozen before any expression value is read.",
        fontsize=5.6, color=BLACK, va="bottom", linespacing=1.45)

fig.savefig(f"{OUT}/fig0_workflow.png", dpi=600, bbox_inches="tight",
            pad_inches=0.02)
fig.savefig(f"{OUT}/fig0_workflow.pdf", bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print("fig0_workflow written")