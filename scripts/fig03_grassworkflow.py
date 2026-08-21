#!/usr/bin/env python3
"""fig03_grassworkflow -- GRASS GIS module-level workflow: the principal modules
at each stage of the Marmara Landsat coastal-transformation chain.
Lane per stage (left spine, numbered, jet-coloured) -> chain of GRASS command
chips (lighter tint of the lane colour). Vector PDF + 300-dpi PNG."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.monospace": ["DejaVu Sans Mono"],
    "mathtext.fontset": "dejavusans",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

INK  = "#1A1A1A"
NOTE = "#404A52"
ARROW = "#6E7A82"

def tint(c, f):   # blend toward white (f in 0..1)
    return tuple(c[i]*(1-f) + f for i in range(3))
def shade(c, f):  # blend toward black
    return tuple(c[i]*(1-f) for i in range(3))

FIGH = 8.7
PPU  = FIGH * 72 / 100.0
def LH(s): return s * 1.2 / PPU

fig, ax = plt.subplots(figsize=(9.0, FIGH))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

# stage lanes: (number, name, [ (module, note), ... ])
stages = [
    (1, "Import &\nregion setup",
        [("g.proj", ""), ("r.import\nr.in.gdal", "Landsat SR"),
         ("v.import", "coast / AOI"), ("g.region\nr.mask", "extent")]),
    (2, "Preprocessing\n& QA",
        [("r.mask", "QA-pixel"), ("r.null", "no-data"),
         ("r.mapcalc", "cloud / shadow mask")]),
    (3, "Composites &\nfeature stack",
        [("r.series", "method=median"), ("i.vi", "NDVI"),
         ("r.mapcalc", "MNDWI \u00b7 NDBI \u00b7 NDWI"), ("i.group", "10-band stack")]),
    (4, "Random Forest\nclassification",
        [("r.learn.train", "fit RF"), ("r.learn.predict", "10 classes"),
         ("r.kappa\nr.report", "accuracy")]),
    (5, "Shoreline &\nreclamation",
        [("r.mapcalc", "MNDWI threshold"), ("r.to.vect", "shoreline"),
         ("r.grow.distance\nv.distance", "EPR (m/yr)"), ("v.overlay", "reclamation")]),
    (6, "Land-cover change\n& urbanization",
        [("r.mapcalc", "2015 \u2194 2025 change"), ("r.stats.zonal", "class areas"),
         ("r.report", "urban gain / loss")]),
    (7, "Mucilage\nsusceptibility",
        [("r.rescale", "shallowness"), ("r.grow.distance", "nutrient proximity"),
         ("r.mapcalc", r"MSI $=\Sigma\, w_i m_i$")]),
    (8, "Composite index\n& spatial analysis",
        [("r.mapcalc", "CTI weighted sum"), ("r.reclass\nr.rescale", "classes"),
         ("r.univar\nr.stats.zonal", "sub-region stats"), ("v.kernel\nr.neighbors", "KDE hotspots")]),
]

# distinct jet colour per stage
JET = plt.get_cmap("jet")
BASE = [JET(v) for v in np.linspace(0.04, 0.96, len(stages))]

# compact vertical layout
top_cy = 89.6
step   = 11.55
lane_h = 10.4
lab_x0, lab_x1 = 2.5, 28.5
chip_x0, chip_x1 = 30.5, 98.6
badge_x = 6.6

def arrow(x1, y1, x2, y2, lw=1.3, color=ARROW):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
        mutation_scale=10, lw=lw, color=color, shrinkA=0, shrinkB=0, zorder=3))

ax.text(2.5, 99.4, "GRASS GIS module-level workflow \u2014 principal modules by stage",
        ha="left", va="top", fontsize=12.5, fontweight="bold", color=shade(BASE[0], 0.35))

cys = [top_cy - i*step for i in range(len(stages))]

for i, ((num, name, chips), cy) in enumerate(zip(stages, cys)):
    b = BASE[i]
    lane_fill = tint(b, 0.62); lane_edge = shade(b, 0.18)
    badge_fill = shade(b, 0.42)
    chip_fill = tint(b, 0.85); chip_edge = tint(b, 0.32)
    mod_col = shade(b, 0.58)

    # lane label
    ax.add_patch(FancyBboxPatch((lab_x0, cy-lane_h/2), lab_x1-lab_x0, lane_h,
        boxstyle="round,pad=0,rounding_size=1.4", facecolor=lane_fill,
        edgecolor=lane_edge, linewidth=1.5, zorder=2))
    ax.add_patch(FancyBboxPatch((badge_x-2.35, cy-2.35), 4.7, 4.7,
        boxstyle="round,pad=0,rounding_size=1.1", facecolor=badge_fill,
        edgecolor="white", linewidth=1.0, zorder=4))
    ax.text(badge_x, cy, str(num), ha="center", va="center", fontsize=11.0,
            fontweight="bold", color="white", zorder=5)
    ax.text((badge_x+2.35+lab_x1)/2 + 0.4, cy, name, ha="center", va="center",
            fontsize=10.0, fontweight="bold", color=INK, zorder=4, linespacing=1.1)

    # chips chain
    n = len(chips); gap = 1.9
    cw = (chip_x1 - chip_x0 - (n-1)*gap) / n
    ch = 8.4
    arrow(lab_x1, cy, chip_x0-0.2, cy, lw=1.3)
    for k, (mod, note) in enumerate(chips):
        x0 = chip_x0 + k*(cw+gap)
        ax.add_patch(FancyBboxPatch((x0, cy-ch/2), cw, ch,
            boxstyle="round,pad=0,rounding_size=1.0", facecolor=chip_fill,
            edgecolor=chip_edge, linewidth=1.2, zorder=2))
        nmod = mod.count("\n")+1
        has_note = bool(note)
        # vertically centre the module(+note) block inside the chip
        block = nmod*LH(8.4) + (LH(7.0)+0.5 if has_note else 0)
        mty = cy + block/2 - 0.2
        ax.text(x0+cw/2, mty, mod, ha="center", va="top", fontsize=8.4,
                family="monospace", fontweight="bold", color=mod_col, zorder=4,
                linespacing=1.08)
        if has_note:
            ny = mty - nmod*LH(8.4) - 0.35
            ax.text(x0+cw/2, ny, note, ha="center", va="top", fontsize=7.0,
                    color=NOTE, zorder=4, linespacing=1.1)
        if k < n-1:
            arrow(x0+cw, cy, x0+cw+gap+0.1, cy, lw=1.1)

    # vertical spine to next lane (coloured by the current stage)
    if i < len(stages)-1:
        arrow(badge_x, cy-lane_h/2-0.15, badge_x, cys[i+1]+lane_h/2+0.15,
              lw=1.8, color=lane_edge)

fig.subplots_adjust(left=0.004, right=0.996, top=0.997, bottom=0.005)
fig.savefig("/mnt/user-data/outputs/fig03_grassworkflow.pdf")
fig.savefig("/mnt/user-data/outputs/fig03_grassworkflow.png", dpi=300)
print("done")
