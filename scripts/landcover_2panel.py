import sys, numpy as np
from PIL import Image
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, "."); import figstyle
figstyle.apply_style()

U = "/mnt/user-data/uploads/"
CROP = (6, 160, 2272, 1668)   # tight: lat labels..map frame right, drop title/legend/credit
img15 = np.array(Image.open(U + "marmara_landcover_utm35n_2015.png").convert("RGB").crop(CROP))
img25 = np.array(Image.open(U + "marmara_landcover_utm35n.png").convert("RGB").crop(CROP))

CLASSES = [  # colour, short name, CLC code, 2015%, 2025%
    ((0,0,128),   "Water",                "CLC 5.1–5.2",    25.3, 25.3),
    ((0,0,255),   "Forest",               "CLC 3.1",         7.3,  5.5),
    ((0,97,255),  "Dense vegetation",     "CLC 3.2.4",       7.3,  7.0),
    ((0,213,255), "Shrub / transitional", "CLC 3.2",         9.4,  8.1),
    ((77,255,170),"Irrigated crops",      "CLC 2.1.2–2.4.2", 5.1,  6.3),
    ((170,255,77),"Cropland",             "CLC 2.1.1",      10.2,  9.8),
    ((255,230,0), "Urban",                "CLC 1.1–1.2",     7.5,  6.6),
    ((255,122,0), "Agric. mosaic",        "CLC 2.4",        11.8, 11.1),
    ((255,19,0),  "Dry grass / fallow",   "CLC 3.2–3.3",    14.0, 12.0),
    ((128,0,0),   "Bare / sparse",        "CLC 3.3.3",       2.2,  8.2),
]

fig = plt.figure(figsize=(13.6, 11.4))
gs = fig.add_gridspec(2, 2, width_ratios=[3.35, 0.82], height_ratios=[1, 1],
                      left=0.006, right=0.994, top=0.928, bottom=0.072,
                      hspace=0.015, wspace=0.004)
for gsi, img, tag, yr in [(gs[0,0], img15, "a", "2015"), (gs[1,0], img25, "b", "2025")]:
    ax = fig.add_subplot(gsi); ax.imshow(img); ax.axis("off")
    ax.set_title(f"({tag}) {yr}", fontsize=14, fontweight="bold", pad=2, loc="left", x=0.01)

# ---- shared legend, packed tight -------------------------------------------
axl = fig.add_subplot(gs[:, 1]); axl.axis("off"); axl.set_xlim(0, 1); axl.set_ylim(0, 1)
PN, PL = 0.66, 0.89         # x of 2015 and 2025 percent columns (right-aligned)
axl.text(0.0, 0.995, "Land cover classes (CORINE Level-3)",
         fontsize=12, fontweight="bold", va="top")
axl.text(PN, 0.958, "2015", fontsize=8.5, style="italic", color="0.35", ha="right", va="top")
axl.text(PL, 0.958, "2025", fontsize=8.5, style="italic", color="0.35", ha="right", va="top")
y = 0.905; dy = 0.092
for (rgb, name, clc, p15, p25) in CLASSES:
    axl.add_patch(Rectangle((0.0, y-0.030), 0.075, 0.060,
                            facecolor=np.array(rgb)/255, edgecolor="0.3", lw=0.5,
                            transform=axl.transAxes, clip_on=False))
    axl.text(0.10, y+0.012, name, fontsize=10.5, va="center")
    axl.text(0.10, y-0.026, clc, fontsize=8, va="center", color="0.45")
    axl.text(PN, y, f"{p15:.1f}", fontsize=9.5, va="center", ha="right")
    axl.text(PL, y, f"{p25:.1f}", fontsize=9.5, va="center", ha="right", color="0.4")
    y -= dy

# ---- title + 2-line, full-width credit -------------------------------------
fig.text(0.5, 0.99, "Sea of Marmara — Land cover, 2015 vs 2025 (WRS-2 180/032 + 181/032)",
         ha="center", va="top", fontsize=17, fontweight="bold")
fig.text(0.5, 0.952, "Unsupervised 10-class classification (ISOCLUST) — UTM Zone 35N / WGS 84",
         ha="center", va="top", fontsize=11)
fig.text(0.006, 0.045,
    "Scenes (Collection-2 L1):  2015 — West 181/032 L8 2015-07-27, East 180/032 L8 2015-07-20;   "
    "2025 — West 181/032 L8 2025-07-22, East 180/032 L9 2025-07-07.\n"
    "Projection UTM 35N / WGS 84, 250 m; clouds masked.   Method: unsupervised ISOCLUST 10-class "
    "(GRASS i.cluster / i.maxlik equivalent).   Source: authors.",
    ha="left", va="bottom", fontsize=12.5, color="0.2", linespacing=1.4)

for e in ("png", "pdf"):
    fig.savefig(f"/mnt/user-data/outputs/fig_landcover_2panel_2015_2025.{e}",
                bbox_inches="tight", pad_inches=0.03, dpi=200 if e == "png" else None)
print("wrote figure")
