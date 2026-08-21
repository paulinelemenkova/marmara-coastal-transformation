#!/usr/bin/env python3
"""fig02_pipeline -- Landsat 2015-2025 processing & analysis pipeline (open-source
GRASS GIS) for the Sea of Marmara coastal-transformation study.
Stages: acquisition -> preprocessing -> composites -> classification ->
shoreline/reclamation -> mucilage susceptibility -> composite index.
Vector PDF + 300-dpi PNG. Pure matplotlib (no external data)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from matplotlib.patheffects import withStroke

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "mathtext.fontset": "dejavusans",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

# ---- palette (muted, colour-blind aware; distinguished by hue + border) -------
INK   = "#22303B"     # text
ARROW = "#5A6673"
C = {
    "input": ("#E4EAF1", "#7C8CA0"),   # open-access inputs
    "prep":  ("#D2E7E0", "#3E8E7E"),   # preprocessing / composites (GRASS)
    "class": ("#CBDCEF", "#3B6FA0"),   # classification
    "shore": ("#DAEAC8", "#5C8C3C"),   # shoreline / reclamation
    "urban": ("#F1E1CA", "#BE8B3E"),   # land-cover change / urban gain
    "muci":  ("#F4D9D6", "#C0655E"),   # mucilage susceptibility
    "index": ("#E4D5EE", "#7E56A3"),   # composite index (emphasis)
    "out":   ("#ECEEF0", "#9AA4AE"),   # outputs
}

fig, ax = plt.subplots(figsize=(7.4, 7.4))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_aspect("equal")
ax.axis("off")

PPU = 7.4 * 72 / 100.0          # points per y data-unit (fig height 7.4in, yspan 100)
def LH(s):                       # rendered line height in data units
    return s * 1.2 / PPU

def box(cx, cy, w, h, key, title, sub=None, eq=None, badge=None,
        tsize=9.0, ssize=6.6, lw=1.1, emph=False):
    fc, ec = C[key]
    ax.add_patch(FancyBboxPatch((cx-w/2, cy-h/2), w, h,
        boxstyle="round,pad=0.0,rounding_size=1.6",
        linewidth=1.8 if emph else lw, facecolor=fc, edgecolor=ec, zorder=2))
    top = cy + h/2
    ty = top - 2.2
    n_t = title.count("\n") + 1
    ax.text(cx, ty, title, ha="center", va="top", fontsize=tsize,
            fontweight="bold", color=INK, zorder=4, linespacing=1.2)
    y = ty - n_t*LH(tsize) - 0.9
    if sub:
        n_s = sub.count("\n") + 1
        ax.text(cx, y, sub, ha="center", va="top", fontsize=ssize,
                color=INK, zorder=4, linespacing=1.3)
        y = y - n_s*LH(ssize)*1.08 - 0.7
    if eq:
        ax.text(cx, y, eq, ha="center", va="top", fontsize=ssize+0.6,
                color="#3A2B47", zorder=4)
    if badge is not None:
        bx, by = cx-w/2+2.9, top-2.9
        ax.add_patch(Circle((bx, by), 1.9, facecolor=ec, edgecolor="white",
                            linewidth=0.8, zorder=5))
        ax.text(bx, by, str(badge), ha="center", va="center", fontsize=7.2,
                fontweight="bold", color="white", zorder=6)

def arrow(x1, y1, x2, y2, lw=1.3, ls="-", color=ARROW, rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
        arrowstyle="-|>", mutation_scale=11, lw=lw, linestyle=ls,
        color=color, shrinkA=0, shrinkB=0, zorder=1,
        connectionstyle=f"arc3,rad={rad}"))

# ---- inputs panel -------------------------------------------------------------
ax.add_patch(FancyBboxPatch((3, 85.2), 94, 12.6,
    boxstyle="round,pad=0,rounding_size=1.4", facecolor="#F5F7F9",
    edgecolor="#C4CCD4", linewidth=1.0, zorder=1))
ax.text(50, 96.4, "Open-access inputs", ha="center", va="top",
        fontsize=8.2, fontweight="bold", color="#54606B", zorder=4)
inputs = [
    ("Landsat 8 OLI\n2015 epoch"), ("Landsat 8/9 OLI\n2025 epoch"),
    ("DEM\nGLO-30 / SRTM"), ("Bathymetry\nEMODnet / GEBCO"),
    ("Reference LC\nWorldCover / CORINE"),
]
xs = [14, 32, 50, 68, 86]
for x, lab in zip(xs, inputs):
    fc, ec = C["input"]
    ax.add_patch(FancyBboxPatch((x-8.0, 89.0-3.3), 16.0, 6.6,
        boxstyle="round,pad=0,rounding_size=1.2", facecolor=fc,
        edgecolor=ec, linewidth=0.9, zorder=2))
    ax.text(x, 89.0, lab, ha="center", va="center", fontsize=6.2,
            color=INK, zorder=4, linespacing=1.3)

# ---- main spine ---------------------------------------------------------------
box(50, 78.0, 62, 8.4, "prep", "Acquisition & footprints", badge=1, tsize=9.0,
    sub="Landsat C2 L2 surface reflectance, 30 m  \u00b7  WRS-2 \u2248 3 scenes\n20 km coastal buffer  \u00b7  2015 & 2025 epochs",
    ssize=6.4)
box(50, 66.8, 62, 8.4, "prep", "Preprocessing (GRASS GIS)", badge=2,
    sub="r.import / r.in.gdal  \u00b7  reproject to Turkish CRS  \u00b7  QA cloud & water masking")
box(50, 55.6, 62, 8.4, "prep", "Median composites & feature stack", badge=3,
    sub="r.series median per epoch  \u00b7  6 bands + NDVI, MNDWI, NDBI, NDWI  (\u2248 10-band stack)")
box(50, 44.4, 62, 8.4, "class", "Random Forest classification", badge=4,
    sub="r.learn.train / r.learn.predict  \u00b7  10-class coastal scheme  \u00b7  accuracy validation")

# ---- three analysis branches --------------------------------------------------
box(18.5, 28.5, 30, 15.0, "shore", "Shoreline &\nreclamation", badge=5, tsize=8.4,
    sub="MNDWI water mask \u00b7 r.to.vect \u00b7\nend-point rate; seaward gain\nflagged as reclamation",
    eq=r"$\mathrm{EPR}=\dfrac{d_{t_2}-d_{t_1}}{t_2-t_1}\;\rightarrow\;R'$", ssize=6.2)
box(50, 28.5, 30, 15.0, "urban", "Land-cover change\n& urbanization", badge=6, tsize=8.4,
    sub="2015 vs 2025 class change \u00b7\nurban / built-up gain \u00b7\nnatural-land loss",
    eq=r"$\rightarrow\;U'$", ssize=6.2)
box(81.5, 28.5, 30, 15.0, "muci", "Mucilage\nsusceptibility", badge=7, tsize=8.4,
    sub="shallowness + weak flushing +\nnutrient / urban proximity \u00b7\noptional FAI detection",
    eq=r"$\mathrm{MSI}=\sum_i w_i m_i\;\rightarrow\;M'$", ssize=6.2)

# ---- composite index ----------------------------------------------------------
box(50, 12.8, 66, 9.6, "index", "Composite coastal-transformation & vulnerability index",
    badge=8, tsize=9.0, emph=True,
    sub="normalize + weight components  \u00b7  reclass into vulnerability classes",
    eq=r"$\mathrm{CTI}=w_r R' + w_u U' + w_m M'$", ssize=6.5)

# ---- outputs ------------------------------------------------------------------
ax.add_patch(FancyBboxPatch((3, 0.6), 94, 6.6,
    boxstyle="round,pad=0,rounding_size=1.2", facecolor=C["out"][0],
    edgecolor=C["out"][1], linewidth=1.0, zorder=1))
ax.text(50, 6.4, "Outputs & spatial analysis", ha="center", va="top",
        fontsize=7.4, fontweight="bold", color="#54606B", zorder=4)
ax.text(50, 4.3,
        "land-cover change maps \u00b7 shoreline-change & reclamation maps \u00b7 mucilage-susceptibility map \u00b7 composite CTI map + sub-region comparison\n"
        "statistical & spatial analysis:  change rate $R$ \u00b7 Mann\u2013Kendall NDVI trend \u00b7 Moran's $I$ \u00b7 KDE hotspots",
        ha="center", va="top", fontsize=6.1, color=INK, linespacing=1.3, zorder=4)

# ---- connectors ---------------------------------------------------------------
arrow(50, 85.2, 50, 82.3)                 # inputs -> stage1
arrow(50, 73.8, 50, 71.05)                # 1 -> 2
arrow(50, 62.6, 50, 59.85)                # 2 -> 3
arrow(50, 51.4, 50, 48.65)                # 3 -> 4
# fan-out from classification to the three branches
for bx, br in [(18.5, 0.28), (50.0, 0.0), (81.5, -0.28)]:
    arrow(50, 40.2, bx, 36.2, rad=br, lw=1.2)
# fan-in from branches to composite index
for bx, br in [(18.5, -0.28), (50.0, 0.0), (81.5, 0.28)]:
    arrow(bx, 21.0, 50, 17.7, rad=br, lw=1.2)
arrow(50, 8.0, 50, 7.4)                    # index -> outputs

fig.subplots_adjust(left=0.01, right=0.99, top=0.995, bottom=0.005)
fig.savefig("/mnt/user-data/outputs/fig02_pipeline.pdf")
fig.savefig("/mnt/user-data/outputs/fig02_pipeline.png", dpi=300)
print("done")
