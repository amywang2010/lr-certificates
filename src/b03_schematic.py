"""Workflow schematic, drawn to a reference style: white shapes, thin black
outlines, orthogonal connectors, horizontal group bars over the column groups and
side labels naming the row groups.

Two failure modes are made impossible by construction rather than by eyeballing:

  1. Every label is short enough to fit its shape. After drawing, each text
     object's RENDERED bounding box is measured with matplotlib's renderer and
     asserted to lie inside its shape's bounding box. Text wider than its shape is
     the failure this check catches, so no label can overflow a stadium.

  2. No two shapes overlap. Shape bounding boxes are checked pairwise.

If either assertion fails the function raises and writes no file, so a broken
figure can never be shipped.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

TEXT_W = 4.8          # inches; LNCS text width, so no scaling at include time
LW = 0.6              # thin black outlines
BLACK = "#1A1A1A"
FS = 5.8              # label font size, points
LEAD = 1.30           # label line spacing

# _shapes holds (patch, text, label) triples so the containment check can run.
_shapes = []


def _register(patch, text, label):
    _shapes.append((patch, text, label))
    return patch


def data_shape(ax, x, y, w, h, label):
    """Parallelogram: an input artifact."""
    s = w * 0.10
    pts = [(x + s, y), (x + w, y), (x + w - s, y + h), (x, y + h)]
    patch = Polygon(pts, closed=True, facecolor="white", edgecolor=BLACK,
                    linewidth=LW)
    ax.add_patch(patch)
    t = ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
                fontsize=FS, linespacing=LEAD)
    return _register(patch, t, "parallelogram")


def box_shape(ax, x, y, w, h, label):
    """Rectangle: a process step."""
    patch = Polygon([(x, y), (x + w, y), (x + w, y + h), (x, y + h)],
                    closed=True, facecolor="white", edgecolor=BLACK,
                    linewidth=LW)
    ax.add_patch(patch)
    t = ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
                fontsize=FS, linespacing=LEAD)
    return _register(patch, t, "rectangle")


def hex_shape(ax, x, y, w, h, label):
    """Hexagon: the test applied to the interval."""
    k = w * 0.16
    pts = [(x + k, y), (x + w - k, y), (x + w, y + h / 2),
           (x + w - k, y + h), (x + k, y + h), (x, y + h / 2)]
    patch = Polygon(pts, closed=True, facecolor="white", edgecolor=BLACK,
                    linewidth=LW)
    ax.add_patch(patch)
    t = ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
                fontsize=FS, linespacing=LEAD)
    return _register(patch, t, "hexagon")


def stadium(ax, x, y, w, h, label):
    """Stadium: the reported output."""
    r = h / 2
    patch = FancyBboxPatch((x + r, y + r), w - 2 * r, h - 2 * r,
                           boxstyle=f"round,pad={r},rounding_size={r}",
                           facecolor="white", edgecolor=BLACK, linewidth=LW,
                           mutation_aspect=1)
    ax.add_patch(patch)
    t = ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
                fontsize=FS, linespacing=LEAD)
    return _register(patch, t, "stadium")


def arrow(ax, pts):
    """Orthogonal polyline with a head on the final segment."""
    for i in range(len(pts) - 1):
        (x1, y1), (x2, y2) = pts[i], pts[i + 1]
        last = i == len(pts) - 2
        ax.add_patch(FancyArrowPatch(
            (x1, y1), (x2, y2), arrowstyle="-|>" if last else "-",
            mutation_scale=5 if last else 1, linewidth=LW, color=BLACK,
            shrinkA=0, shrinkB=0, joinstyle="miter"))


def group_bar(ax, x0, x1, y, label):
    """Horizontal bar over a column group, broken for its label."""
    ax.plot([x0, x1], [y, y], color=BLACK, lw=0.6)
    inset = 0.10 * (x1 - x0)
    ax.plot([x0 + inset, x1 - inset], [y, y], color="white", lw=2.2)
    ax.text((x0 + x1) / 2, y + 1.6, label, ha="center", va="bottom",
            fontsize=6.2)


def side_label(ax, x, y, label):
    """Row-group label down the left margin."""
    ax.text(x, y, label, ha="left", va="center", fontsize=6.2, rotation=0)


def _verify(fig, ax):
    """Measure rendered text and assert it fits; assert shapes do not overlap.

    Overflow is reported in axis grid units (the 0-100 space the layout is
    written in), so the numbers are directly comparable to the box sizes in
    build().
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px_per_unit_x = fig.get_size_inches()[0] * fig.dpi / 100.0
    px_per_unit_y = fig.get_size_inches()[1] * fig.dpi / 100.0
    problems = []

    for patch, text, kind in _shapes:
        tb = text.get_window_extent(renderer=renderer)
        pb = patch.get_window_extent(renderer=renderer)
        # Overflow is how far the TEXT sticks out past the SHAPE on each side.
        over_left = pb.x0 - tb.x0
        over_right = tb.x1 - pb.x1
        over_bot = pb.y0 - tb.y0
        over_top = tb.y1 - pb.y1
        over_w = max(0.0, over_left, over_right) / px_per_unit_x
        over_h = max(0.0, over_bot, over_top) / px_per_unit_y
        if over_w > 0.4 or over_h > 0.4:
            problems.append(
                f"  TEXT OVERFLOWS ITS {kind.upper()} by "
                f"{over_w:.1f} x {over_h:.1f} grid units "
                f"(L={over_left / px_per_unit_x:.1f} R={over_right / px_per_unit_x:.1f} "
                f"B={over_bot / px_per_unit_y:.1f} T={over_top / px_per_unit_y:.1f}) -> "
                f"{text.get_text()!r}")

    boxes = [(p.get_window_extent(renderer=renderer), lab)
             for p, _, lab in _shapes]
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][0], boxes[j][0]
            ox = (min(a.x1, b.x1) - max(a.x0, b.x0)) / px_per_unit_x
            oy = (min(a.y1, b.y1) - max(a.y0, b.y0)) / px_per_unit_y
            if ox > 0.5 and oy > 0.5:
                problems.append(
                    f"  SHAPES OVERLAP by {ox:.1f}x{oy:.1f} grid units: "
                    f"{boxes[i][1]!r} <-> {boxes[j][1]!r}")

    if problems:
        raise AssertionError(
            "schematic layout check failed; no file written:\n"
            + "\n".join(problems))


def build(outdir):
    """Draw, verify, and write the schematic. Returns the check report."""
    global _shapes
    _shapes = []

    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
    })
    fig, ax = plt.subplots(figsize=(TEXT_W, 2.30))
    # Axes must span the whole figure, otherwise grid units are not the layout
    # units, and the containment arithmetic below misbehaves.
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # ---- column groups: horizontal bars across the top
    group_bar(ax, 11, 30, 90, "Inputs")
    group_bar(ax, 33, 52, 90, "Uncertainty set")
    group_bar(ax, 55, 74, 90, "Interval machinery")
    group_bar(ax, 77, 96, 90, "Inference")

    # ---- row groups: side labels down the left margin
    side_label(ax, 0.5, 70, "Observed")
    side_label(ax, 0.5, 46, "Permuted")

    # ---- observed path
    data_shape(ax, 11, 62, 19, 16,
               "Molecules\nvendor masks\nassignment $A_0$")
    box_shape(ax, 33, 60, 19, 20,
              "Build $U$\nband $3\\,\\mu$m\ndonor $6\\,\\mu$m")
    # ---- permuted path
    data_shape(ax, 11, 38, 19, 16,
               "Region $\\times$\ndecile blocks\nfrozen from $A_0$")
    box_shape(ax, 33, 36, 19, 20,
              "Matched null\nlabels permuted\n$B = 1000$")
    # ---- shared machinery (both rows enter here)
    box_shape(ax, 55, 46, 19, 24,
              "Per-gene bounds\n$T_{\\rm lo}$, $T_{\\rm hi}$\nover $U$")
    hex_shape(ax, 77, 48, 19, 20,
              "Interval sign\nand matched\nnull agree?")
    # ---- reported output
    stadium(ax, 55, 14, 41, 15,
            "Certificate" "\n" "certified $+$ or $-$, or" "\n"
            "non-identifiable")

    # ---- orthogonal connectors
    arrow(ax, [(30, 70), (33, 70)])                      # observed: A0 -> U
    arrow(ax, [(30, 46), (33, 46)])                      # permuted: blocks -> null
    arrow(ax, [(52, 70), (53.5, 70), (53.5, 62), (55, 62)])   # U -> machinery
    arrow(ax, [(52, 46), (53.5, 46), (53.5, 54), (55, 54)])   # null -> machinery
    arrow(ax, [(74, 58), (77, 58)])                      # machinery -> test
    arrow(ax, [(86.5, 48), (86.5, 29)])                  # test -> output

    ax.text(11, 4.0,
            "Both inputs enter the same interval machinery, so the null is computed on "
            "the same layer as the observed bound.",
            fontsize=5.2, color=BLACK, va="bottom", ha="left", linespacing=1.4)

    _verify(fig, ax)

    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, "fig0_workflow.png")
    pdf = os.path.join(outdir, "fig0_workflow.pdf")
    fig.savefig(png, dpi=600, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.02,
                metadata={"CreationDate": None})
    plt.close(fig)
    return f"{len(_shapes)} shapes, all labels contained, no shape overlaps"


if __name__ == "__main__":
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print("layout check:", build(os.path.join(ROOT, "figures")))