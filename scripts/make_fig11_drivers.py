#!/usr/bin/env python3
"""
fig11_drivers -- Conceptual matrix of Sea of Marmara coastal-transformation
drivers (physical / socio-economic / institutional) per coastal sub-region.

Ordinal driver intensity (Negligible / Low / Moderate / High) is assigned from
the qualitative evidence in the manuscript (sub-region characteristics table,
reclamation and mucilage-susceptibility sections). Spectral discrete palette
(one distinct colour per class); each cell also carries its level word so the
figure stays readable in grayscale.

Dependencies: matplotlib + numpy only.
Outputs: fig11_drivers.pdf (vector, for LaTeX) and fig11_drivers.png (600 dpi).
"""
import os
import textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, FancyBboxPatch
import matplotlib.transforms as mtransforms

OUTDIR = "."   # change to "figures" to write straight into the LaTeX figures/ dir

# ---- house style (self-contained; no external helper needed) ---------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.linewidth": 0.8,
    "mathtext.default": "regular",
    "savefig.dpi": 600, "figure.dpi": 120,
})

# ----------------------------------------------------------------------------
# Sub-regions (columns) and drivers (rows), grouped by family.
# Intensity codes: 0 Negligible, 1 Low, 2 Moderate, 3 High.
# Order per row = [Istanbul, Izmit, Southern Marmara, Tekirdag/Canakkale].
# ----------------------------------------------------------------------------
regions = ["Istanbul\nmetropolitan", "Gulf of\nIzmit",
           "Southern\nMarmara bays", "Tekirdag /\nCanakkale"]

groups = [
    ("Physical", [
        ("Enclosed, low-flushing geometry", [2, 3, 3, 1]),
        ("Shallow, stratified water (mucilage)", [2, 3, 3, 1]),
        ("Sediment / deltaic supply", [1, 1, 3, 2]),
    ]),
    ("Socio-economic", [
        ("Urban population density", [3, 3, 2, 2]),
        ("Industry, refineries & ports", [2, 3, 1, 2]),
        ("Land reclamation / coastal fills", [3, 3, 1, 0]),
        ("Aquaculture & agriculture", [1, 1, 3, 2]),
    ]),
    ("Institutional", [
        ("Waterfront-park (dolgu) policy", [3, 2, 1, 0]),
        ("Wastewater / effluent loading", [2, 3, 2, 1]),
        ("Strait shipping traffic", [2, 1, 1, 3]),
    ]),
]

row_labels, group_spans, M = [], [], []
r = 0
for gname, drivers in groups:
    start = r
    for dlabel, vals in drivers:
        row_labels.append(dlabel)
        M.append(vals)
        r += 1
    group_spans.append((gname, start, r))            # [start, end)
M = np.array(M, dtype=float)
nrow, ncol = M.shape

# ----------------------------------------------------------------------------
# Spectral palette: one distinct colour per ordinal class.
# ----------------------------------------------------------------------------
levels = ["Negligible", "Low", "Moderate", "High"]
spec = plt.get_cmap("Spectral")
cell_colors = [spec(x) for x in (0.90, 0.63, 0.32, 0.07)]  # blue, green, orange, red
cmap = ListedColormap(cell_colors)


def text_color(rgba):
    """Dark or white label, whichever contrasts with the cell colour."""
    r, g, b = rgba[:3]
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "0.12" if lum > 0.55 else "white"


# ----------------------------------------------------------------------------
# Figure geometry: wide columns, wide left gutter for the group labels.
# ----------------------------------------------------------------------------
FW, FH = 9.4, 5.7
L, Rr, T, B = 0.345, 0.795, 0.895, 0.065             # axes margins (fig fraction)

fig, ax = plt.subplots(figsize=(FW, FH))
fig.subplots_adjust(left=L, right=Rr, top=T, bottom=B)

ax.imshow(M, cmap=cmap, vmin=-0.5, vmax=3.5, aspect="auto",
          interpolation="nearest")

# level word in every cell (redundant with colour -> grayscale-safe)
for i in range(nrow):
    for j in range(ncol):
        lv = int(M[i, j])
        ax.text(j, i, levels[lv], ha="center", va="center",
                fontsize=8.2, color=text_color(cell_colors[lv]))

# thin white gridlines between cells
ax.set_xticks(np.arange(-0.5, ncol, 1), minor=True)
ax.set_yticks(np.arange(-0.5, nrow, 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.4)
ax.tick_params(which="minor", length=0)
ax.tick_params(which="major", length=0)

# column (sub-region) labels on top
ax.xaxis.set_ticks_position("top")
ax.set_xticks(range(ncol))
ax.set_xticklabels(regions, fontsize=9, fontweight="bold")
ax.xaxis.set_label_position("top")

# row (driver) labels on the left
ax.set_yticks(range(nrow))
ax.set_yticklabels(row_labels, fontsize=8.2)

for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_xlim(-0.5, ncol - 0.5)
ax.set_ylim(nrow - 0.5, -0.5)

# ----------------------------------------------------------------------------
# Group brackets + rotated titles in a far-left gutter (fig-x / data-y blend),
# placed well left of the longest row label so nothing overlaps.
# ----------------------------------------------------------------------------
blend = mtransforms.blended_transform_factory(fig.transFigure, ax.transData)
group_colors = {"Physical": "#1b4965", "Socio-economic": "#5f0f40",
                "Institutional": "#0b6e4f"}
BAR_X, LAB_X = 0.052, 0.022                          # figure-fraction x positions
for gname, start, end in group_spans:
    if start != 0:
        ax.axhline(start - 0.5, color="0.25", lw=1.6)  # separator between families
    ymid = (start + end - 1) / 2.0
    ax.plot([BAR_X, BAR_X], [start - 0.42, end - 0.58], transform=blend,
            color=group_colors[gname], lw=3.4, solid_capstyle="butt",
            clip_on=False)
    ax.text(LAB_X, ymid, gname, transform=blend, rotation=90,
            ha="center", va="center", fontsize=9.5, fontweight="bold",
            color=group_colors[gname], clip_on=False)

# ----------------------------------------------------------------------------
# Legend (ordinal intensity), top-right in empty space.
# ----------------------------------------------------------------------------
handles = [Patch(facecolor=cell_colors[k], edgecolor="0.5", lw=0.5,
                 label=levels[k]) for k in range(4)]
leg = ax.legend(handles=handles, title="Driver intensity", loc="upper left",
                bbox_to_anchor=(1.02, 1.0), fontsize=8.5, title_fontsize=9,
                frameon=True, borderpad=0.7, handlelength=1.4, handleheight=1.4,
                labelspacing=0.6)
leg.get_frame().set_linewidth(0.6)
leg.get_frame().set_edgecolor("0.6")

# ----------------------------------------------------------------------------
# Boxed note under the legend: same width as the legend, taller, wrapped,
# same contour line.
# ----------------------------------------------------------------------------
fig.canvas.draw()
inv = fig.transFigure.inverted()
lb = leg.get_window_extent()
(lx0, ly0) = inv.transform((lb.x0, lb.y0))
(lx1, ly1) = inv.transform((lb.x1, lb.y1))
box_w = lx1 - lx0                                    # match legend width exactly
box_left = lx0
box_top = ly0 - 0.028                                # just below the legend

note = ("Dark cells recurring down a column mark sub-regions where physical, "
        "socio-economic and institutional pressures superimpose and reinforce "
        "one another.")
pad_pt = 6.0
usable_pt = box_w * FW * 72 - 2 * pad_pt
ncar = max(12, int(usable_pt / 4.55))                # ~mean glyph advance at 8 pt
wrapped = textwrap.fill(note, width=ncar)
nlines = wrapped.count("\n") + 1
box_h = (nlines * 8.0 * 1.32 + 2 * pad_pt) / (FH * 72)

fbox = FancyBboxPatch((box_left, box_top - box_h), box_w, box_h,
                      boxstyle="square,pad=0", transform=fig.transFigure,
                      facecolor="white", edgecolor="0.6", linewidth=0.6,
                      clip_on=False, zorder=3)
fig.add_artist(fbox)
fig.text(box_left + pad_pt / (FW * 72), box_top - pad_pt / (FH * 72), wrapped,
         ha="left", va="top", fontsize=8.0, color="0.2", zorder=4)

# credit line (8 pt)
fig.text(0.010, 0.010,
         "Software: Python 3 / Matplotlib %s. "
         "Data: manuscript sub-region assessment (authors). Source: authors."
         % matplotlib.__version__,
         fontsize=8, color="0.35", ha="left", va="bottom")

# ----------------------------------------------------------------------------
# Export vector PDF + 600-dpi PNG.
# ----------------------------------------------------------------------------
os.makedirs(OUTDIR, exist_ok=True)
for ext, kw in (("pdf", {}), ("png", {"dpi": 600})):
    fig.savefig(os.path.join(OUTDIR, f"fig11_drivers.{ext}"),
                bbox_inches="tight", pad_inches=0.02, **kw)
print("wrote fig11_drivers.pdf and fig11_drivers.png")
