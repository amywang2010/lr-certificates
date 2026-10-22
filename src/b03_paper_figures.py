"""B03 paper figures, generated only from the released artifacts (hashed in the manifest).

fig1_main.pdf    two panels, one column:
                   (a) breast certified intervals per (pair, region): T(A0) point plus
                       [T_lo, T_hi] bar, colored by verdict, marker shape coding region.
                   (b) cross-tissue certificates for the shared pairs and their tissue
                       analogues, with the Proseg point estimate and the excursion marked.
fig2_support.pdf two panels, side by side:
                   (a) cross-segmentation validation: fraction of certified rows inside
                       their certified interval, with the coverage gate marked.
                   (b) kernel-margin distributions c* for the sections with certified rows.
fig0_workflow.pdf  pipeline schematic (b03_schematic.py, which carries a layout check).

Each figure is drawn at the LNCS text width (4.8 in) and included at \\textwidth, so no
scaling occurs and the stated point sizes are the printed point sizes. Text is embedded
as TrueType (pdf.fonttype 42) so the figures stay crisp and searchable and no Type 3
fonts enter the submission PDF. Every plotted number is read from the artifact CSVs;
assertions verify row counts against the released summaries before plotting.
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
LUN = os.path.join(ROOT, "lung")   # released lung artifacts
OUT = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.linewidth": 0.7, "lines.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "out", "ytick.direction": "out",
    "figure.dpi": 150,
})

TEXT_W = 4.8

C_POS = "#D55E00"    # claim-certified positive
C_NEG = "#0072B2"    # claim-certified negative
C_NON = "#B3B3B3"    # interval sign only, no claim
C_PRO = "#1A1A1A"    # Proseg point estimate
GREY = "#ECECEC"     # alternating pair band
REGION_MARK = {0: "o", 1: "s", 2: "^"}


def plabel(ax, letter, size=9):
    ax.text(0.004, 0.985, letter, transform=ax.transAxes, ha="left", va="top",
            fontweight="bold", fontsize=size)


MIN_LEGEND_GAP_MM = 2.0


def legend_gap_mm(ax, leg):
    """Smallest visual gap, in printed millimetres, between legend and data ink.

    Measured against the full drawn extent of every interval and point estimate,
    inflated by half a line width, so a bar running under the legend in its
    middle counts the same as one touching it at its end. Two clearances are
    demanded, one for each way a legend and a bar can appear to collide: an
    artist that shares the legend's height band must stand clear of it
    horizontally, and an artist that shares its width must stand clear of it
    vertically. Artists that share neither are diagonally separated and cannot
    touch it.

    The figure is included in the manuscript at its natural size, so a
    millimetre here is a millimetre on the page. Non-overlap alone is not
    enough: a legend edge a tenth of a millimetre from the end of a bar reads as
    a collision, which is how the first version of this figure failed.
    """
    fig = ax.figure
    fig.canvas.draw()
    eb = leg.get_window_extent(fig.canvas.get_renderer())
    mm_per_px = 25.4 / fig.dpi
    worst, where = float("inf"), None
    for ln in ax.get_lines():
        xy = ln.get_xydata()
        ls, mk = ln.get_linestyle(), ln.get_marker()
        if mk not in ("None", "none", "", None):
            d, half = ax.transData.transform(xy), ln.get_markersize()
            tag = f"point estimate at ({xy[0][0]:+.2f}, {xy[0][1]:+.2f})"
        elif len(xy) == 2 and ls not in ("None", "none", ""):
            d = ax.transData.transform(xy)
            half = max(1.5, ln.get_linewidth() * fig.dpi / 72.0)
            tag = f"interval [{xy[0][0]:+.2f}, {xy[-1][0]:+.2f}] at y={xy[0][1]:+.2f}"
        else:
            continue
        x0, x1 = d[:, 0].min() - half, d[:, 0].max() + half
        y0, y1 = d[:, 1].min() - half, d[:, 1].max() + half
        x_ov = x1 > eb.x0 and x0 < eb.x1
        y_ov = y1 > eb.y0 and y0 < eb.y1
        if x_ov and y_ov:
            gap = -1.0
        elif y_ov:
            gap = max(eb.x0 - x1, x0 - eb.x1)
        elif x_ov:
            gap = max(eb.y0 - y1, y0 - eb.y1)
        else:
            continue
        if gap < worst:
            worst, where = gap, tag
    return worst * mm_per_px, where


def check_legend(ax, leg, name, min_gap_mm=MIN_LEGEND_GAP_MM):
    """Assert the legend clears every drawn interval and point estimate."""
    gap, where = legend_gap_mm(ax, leg)
    if gap < min_gap_mm:
        raise AssertionError(
            f"{name} legend is {gap:.2f} mm from the nearest drawn artist "
            f"({where}), below the {min_gap_mm:g} mm minimum; no file written")
    return gap

def check_layout(fig, name, pad=1.5):
    """Assert that no drawn text leaves the canvas or collides with other text.

    Overflow and collision are the two ways a crowded figure fails. Both are
    measured from rendered bounding boxes, so the check aborts instead of
    writing a figure that cannot be read.
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    figbox = fig.bbox
    texts = []
    for ax in fig.axes:
        xl, xh = ax.get_xlim()
        yl, yh = ax.get_ylim()
        live = list(ax.texts) + [ax.xaxis.label, ax.yaxis.label]
        # ticks outside the view are not drawn, so ignore them
        live += [t for t in ax.get_xticklabels() if xl <= t.get_position()[0] <= xh]
        live += [t for t in ax.get_yticklabels() if yl <= t.get_position()[1] <= yh]
        for t in live:
            if t.get_text().strip():
                texts.append((t, t.get_window_extent(renderer=rend)))
    for leg in list(fig.legends) + [a.get_legend() for a in fig.axes if a.get_legend()]:
        if leg is not None:
            for t in leg.get_texts():
                if t.get_text().strip():
                    texts.append((t, t.get_window_extent(renderer=rend)))
    problems = []
    for t, b in texts:
        if (b.x0 < figbox.x0 - pad or b.x1 > figbox.x1 + pad
                or b.y0 < figbox.y0 - pad or b.y1 > figbox.y1 + pad):
            problems.append(f"  TEXT OFF-CANVAS: {t.get_text()!r} at "
                            f"({b.x0:.0f},{b.y0:.0f})-({b.x1:.0f},{b.y1:.0f}) "
                            f"vs canvas {figbox.width:.0f}x{figbox.height:.0f}")
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            a, b = texts[i][1], texts[j][1]
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ox > pad and oy > pad:
                problems.append(f"  TEXT OVERLAP {ox:.1f}x{oy:.1f}px: "
                                f"{texts[i][0].get_text()!r} <-> {texts[j][0].get_text()!r}")
    if problems:
        raise AssertionError(f"{name} layout check failed; no file written:\n"
                             + "\n".join(problems))


def verdict_color(pos, neg):
    return C_POS if pos else (C_NEG if neg else C_NON)


# ================================================================ fig1_main
fin = pd.read_csv(f"{RES}/scout_final_matched.csv")
assert len(fin) == 75
assert int(fin.certified_pos_matched.sum()) == 37
assert int(fin.certified_neg_matched.sum()) == 9

fin["verdict"] = np.where(fin.certified_pos_matched, "pos",
                  np.where(fin.certified_neg_matched, "neg", "non"))
pair_meta = (fin.groupby("pair")
                .agg(ctrl=("control", "first"), med=("T_lo", "median")))
pair_meta = pair_meta.sort_values(["ctrl", "med"], ascending=[False, True])
pairs = list(pair_meta.index)
ypos = {p: i for i, p in enumerate(pairs)}

f3 = pd.read_csv(f"{RES}/paper_fig3_data.csv")
assert len(f3) == 18

xlo = float(min(fin.T_lo.min(), f3.T_lo.min()))
xhi = float(max(fin.T_hi.max(), f3.T_hi.max()))
pad = 0.04 * (xhi - xlo)

fig = plt.figure(figsize=(TEXT_W, 4.15))
gs = GridSpec(2, 1, figure=fig, height_ratios=[2.9, 1.15], hspace=0.50,
              left=0.235, right=0.955, top=0.975, bottom=0.115)

# ---- (a) breast certified intervals
ax = fig.add_subplot(gs[0])
for i in range(len(pairs)):
    if i % 2:
        ax.axhspan(i - 0.5, i + 0.5, color=GREY, zorder=0)
for _, r in fin.iterrows():
    y = ypos[r.pair] + (r.region - 1) * 0.30
    ax.plot([r.T_lo, r.T_hi], [y, y], color=verdict_color(r.certified_pos_matched,
            r.certified_neg_matched), lw=1.5, solid_capstyle="butt", zorder=3)
    ax.plot(r.T0, y, marker=REGION_MARK[int(r.region)], ms=2.4, mfc="white",
            mec="black", mew=0.5, ls="none", zorder=4)

tr = blended_transform_factory(ax.transAxes, ax.transData)
for p, y in ypos.items():
    ctrl = bool(pair_meta.loc[p, "ctrl"])
    ax.text(-0.012, y, ("ctrl  " if ctrl else "") + p.replace("->", "\u2192"),
            transform=tr, ha="right", va="center", fontsize=5.9,
            fontweight="bold" if ctrl else "normal", color="#1A1A1A")
ax.set_yticks([])
ax.set_xlim(xlo - pad, xhi + pad)
ax.set_ylim(-0.7, len(pairs) - 0.3)
ax.axvline(0, color="black", lw=0.7, ls=(0, (3, 2)), zorder=2)
ax.set_xlabel("$T$ (log$_2$ expression-product ratio)", fontsize=7.8)
plabel(ax, "a")
handles = [Line2D([], [], color=C_POS, lw=2.2, label="claim-certified $+$"),
           Line2D([], [], color=C_NEG, lw=2.2, label="claim-certified $-"),
           Line2D([], [], color=C_NON, lw=2.2, label="no claim")]
for r in (0, 1, 2):
    handles.append(Line2D([], [], color="black", marker=REGION_MARK[r], ls="none",
                          ms=3.4, mfc="white", mec="black", mew=0.6,
                          label=f"region {r}"))
# Placement is measured with the same metric the build gate uses. In this panel the
# lower-right corner is the worst place for a legend: the widest intervals, which
# are the most negative pairs, sit in the bottom rows and reach furthest right,
# leaving 0.2 mm between the longest bar and the legend edge. Sweeping the anchor
# with the same metric the build gate uses peaks at y = 0.33, where the legend
# clears the band below it (T_hi = +6.13 at y = +4.30) by 6.7 mm and the control
# bands above it by the same margin. check_legend fails the build below 2 mm.
leg_a = ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=6.4,
                  ncol=1, handletextpad=0.5, borderaxespad=0.2,
                  bbox_to_anchor=(1.0, 0.33), bbox_transform=ax.transAxes)

# ---- (b) cross-tissue certificates
ax = fig.add_subplot(gs[1])
rows = [("Breast", ["PECAM1->KDR", "CDH1->ERBB2", "ESR1->PGR"], "breast"),
        ("Lung", ["PECAM1->KDR", "CDH1->EGFR", "ESR1->PGR"], "lung")]
ylabels, ypos_b, group_y = [], [], {}
segs_b = []
y = 0
excursion_xy = None
for tname, tlist, tkey in rows:
    group_y[tname] = []
    for pair in tlist:
        sub = f3[(f3.tissue == tkey) & (f3.pair == pair)].sort_values("region")
        if sub.empty:
            continue
        for _, r in sub.iterrows():
            yy = y + (r.region - 1) * 0.26
            ax.plot([r.T_lo, r.T_hi], [yy, yy],
                    color=verdict_color(r.certified_pos_matched, r.certified_neg_matched),
                    lw=1.9, solid_capstyle="butt", zorder=3)
            ax.plot(r.T_proseg, yy, "D", ms=2.6, color=C_PRO, mec="white", mew=0.35,
                    ls="none", zorder=4)
            if tkey == "lung" and pair == "CDH1->EGFR" and int(r.region) == 0:
                excursion_xy = (float(r.T_proseg), yy)
        ylabels.append(f"{pair.split('->')[0]} \u2192 {pair.split('->')[1]}")
        ypos_b.append(y)
        group_y[tname].append(y)
        y += 1
    y += 1.8
ax.set_yticks(ypos_b)
ax.set_yticklabels(ylabels, fontsize=6.0)
ax.axvline(0, color="black", lw=0.7, ls=(0, (3, 2)))
ax.set_xlim(xlo - pad, xhi + pad)
# group headings sit in the left label column, in the gap between groups
trb = blended_transform_factory(ax.transAxes, ax.transData)
for tname, ys in group_y.items():
    ax.text(-0.012, max(ys) + 1.0, tname, transform=trb, ha="right", va="center",
            fontweight="bold", fontsize=7, color="#1A1A1A")
# mark the single cross-platform excursion
if excursion_xy is not None:
    ax.annotate("excursion beyond\nthe certified interval", xy=excursion_xy,
                xytext=(0.985, 0.97), textcoords="axes fraction", fontsize=6.0,
                ha="right", va="top", color="#1A1A1A",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="#1A1A1A",
                                shrinkA=2, shrinkB=3))
ax.set_xlabel("$T$ (log$_2$ expression-product ratio)", fontsize=7.8)
ax.set_ylim(-1.2, y - 0.5)
leg_b = ax.legend(handles=[Line2D([], [], color=C_PRO, marker="D", ls="none", ms=3.4,
                                 mec="white", mew=0.35, label="Proseg $T$ (point estimate)")],
                  loc="center right", frameon=False, fontsize=6.2, handletextpad=0.5,
                  borderaxespad=0.3)
plabel(ax, "b")
for _ax, _leg, _nm in ((fig.axes[0], leg_a, "fig1_main(a)"),
                       (fig.axes[1], leg_b, "fig1_main(b)")):
    print(f"  {_nm}: legend clears the nearest drawn artist by "
          f"{check_legend(_ax, _leg, _nm):.2f} mm")
check_layout(fig, "fig1_main")
fig.savefig(f"{OUT}/fig1_main.png", dpi=600)
fig.savefig(f"{OUT}/fig1_main.pdf", metadata={"CreationDate": None})
plt.close(fig)

# ================================================================ fig2_support
s2 = json.load(open(f"{RES}/phase2_proseg_summary.json"))
s2l = json.load(open(f"{LUN}/phase5_e2_summary.json"))
kb = pd.read_csv(f"{RES}/phase4_kernel_margins.csv")
kl = pd.read_csv(f"{LUN}/phase4_kernel_margins.csv")
assert len(kb) == 46 and len(kl) == 7

# E2 coverage over the certified rows of each leg, the rows that carry a claim. The breast
# summary stores its primary over all evaluable rows (72/75); the certified-row count is
# 43/46 in the same coverage table, and both are asserted here.
cb = pd.read_csv(f"{RES}/phase2_proseg_coverage.csv")
cl = pd.read_csv(f"{LUN}/phase5_e2_coverage.csv")
cert_b = cb.certified_pos_matched | cb.certified_neg_matched
cert_l = cl.cert_pos_final | cl.cert_neg_final
n_in = [int(cb.inside[cert_b].sum()), int(cl.inside[cert_l].sum())]
n_cert = [int(cert_b.sum()), int(cert_l.sum())]
assert n_cert == [s2["certified_rows"], s2l["certified_rows"]] == [46, 7]
assert n_in == [43, 6]
assert int(cb.inside.sum()) == 72 and len(cb) == 75 and len(cl) == 54

fig, axes = plt.subplots(1, 2, figsize=(TEXT_W, 1.70))
fig.subplots_adjust(left=0.175, right=0.975, top=0.89, bottom=0.255, wspace=0.42)

ax = axes[0]
vals = [100 * a / b for a, b in zip(n_in, n_cert)]
bars = ax.bar([0, 1], vals, width=0.5, color=["#6A3D9A", "#6A3D9A"], alpha=0.9,
              edgecolor="#1A1A1A", linewidth=0.5, zorder=3)
for x, v, n in zip([0, 1], vals, [f"{a}/{b}" for a, b in zip(n_in, n_cert)]):
    ax.text(x, v + 2.2, f"{v:.1f}%", ha="center", va="bottom", fontsize=7)
    ax.text(x, 3, n, ha="center", va="bottom", fontsize=6.4, color="white")
ax.axhline(90, color="#D55E00", lw=0.9, ls=(0, (4, 2)), zorder=4)
ax.text(0.5, 90.8, "gate 90%", fontsize=6.4, ha="center", va="bottom",
        color="#D55E00")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Breast", "Lung"])
ax.set_ylabel("certified rows inside\nProseg interval (%)", fontsize=7.5)
ax.set_ylim(0, 108)
ax.tick_params(axis="x", length=0)
plabel(ax, "a")

ax = axes[1]
rng = np.random.default_rng(20260907)
data = [kb.c_star.to_numpy(), kl.c_star.to_numpy()]
labels = ["Breast\n(46)", "Lung\n(7)"]
for i, d in enumerate(data):
    x = i + 1
    ax.vlines(x, np.quantile(d, 0.25), np.quantile(d, 0.75), color="#1A1A1A",
              lw=4.5, alpha=0.18, zorder=2)
    ax.hlines(np.median(d), x - 0.17, x + 0.17, color="#1A1A1A", lw=1.5, zorder=3)
    jit = rng.uniform(-0.10, 0.10, len(d))
    ax.plot(x + jit, d, "o", ms=2.6, mfc="none", mec="#4D4D4D", mew=0.6, zorder=4)
    ax.text(x, np.max(d) * 1.12, f"med {np.median(d):g}", ha="center", va="bottom",
            fontsize=6.2, color="#1A1A1A")
ax.set_yscale("log")
ax.set_xticks([1, 2]); ax.set_xticklabels(labels)
ax.set_ylabel(r"critical kernel multiplier  $c^*$", fontsize=7.5)
ax.axhline(1.0, color="#1A1A1A", lw=0.7, ls=(0, (1, 2)))
ax.text(0.02, 1.08, r"$c=1$", transform=ax.get_yaxis_transform(), fontsize=6.2,
        va="bottom", ha="left")
ax.set_xlim(0.55, 2.45)
ax.set_ylim(0.85, 60)
ax.tick_params(axis="x", length=0)
plabel(ax, "b")
check_layout(fig, "fig2_support")
fig.savefig(f"{OUT}/fig2_support.png", dpi=600)
fig.savefig(f"{OUT}/fig2_support.pdf", metadata={"CreationDate": None})
plt.close(fig)

print("figures written:", sorted(f for f in os.listdir(OUT) if f.endswith(".pdf")))

# ================================================================ fig0_workflow
# The schematic carries a layout self-check (every label measured against its shape,
# pairwise shape-overlap check) and aborts rather than writing a broken figure.
import b03_schematic  # noqa: E402

print("fig0 schematic:", b03_schematic.build(OUT))
