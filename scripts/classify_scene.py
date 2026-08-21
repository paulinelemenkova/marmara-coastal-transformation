#!/usr/bin/env python3
"""
classify_scene.py
-----------------
Unsupervised classification of a Landsat true-colour scene into 1 water + 9 land
classes (k-means; GRASS GIS i.cluster/i.maxlik equivalent), rendered as a map PNG
with a true lat/lon graticule (from scene corner metadata), a single legend block
(colour + name + CORINE-referenced explanation) and place labels. Clouds are
masked (excluded from clustering, drawn as no-data).

Configured scenes: 2015 (Landsat 8) and 2025 (Landsat 9), southern Marmara
(WRS-2 180/032). Select with SCENE below. Corner coordinates and Scene IDs come
from the USGS Collection-2 metadata (table.docx).
"""
import glob
import numpy as np
from PIL import Image
from scipy import ndimage
from sklearn.cluster import KMeans
from pyproj import Transformer
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib import font_manager as fm

# --- Helvetica (free metric-compatible TeX Gyre Heros; real Helvetica if present) --
for _f in glob.glob("/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreheros-*.otf"):
    fm.fontManager.addfont(_f)
mpl.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Helvetica", "TeX Gyre Heros", "Nimbus Sans",
                                         "Arial", "DejaVu Sans"]})

# ============================ scene configuration ============================ #
SCENES = {
    "2025": dict(
        src="/mnt/user-data/uploads/LC09_L1TP_180032_20250707_20250709_02_T1.jpg",
        out="LC09_180032_20250707_lc10.png",
        sat="Landsat 9", date_t="2025-07-07", date_a="2025/07/07",
        scene_id="LC91800322025188LGN00", cloud=False, swaps=[(8, 9)],
        corners=[(27.29545, 41.39095), (30.04212, 41.35110),
                 (27.28641, 39.28822), (29.94927, 39.25120)]),
    "2015": dict(
        src="/mnt/user-data/uploads/LC08_L1TP_180032_20150720_20200908_02_T1.jpg",
        out="LC08_180032_20150720_lc10.png",
        sat="Landsat 8", date_t="2015-07-20", date_a="2015/07/20",
        scene_id="LC81800322015201LGN01", cloud=True, swaps=[], match_to="2025",
        corners=[(27.30262, 41.39093), (30.04570, 41.35101),
                 (27.29337, 39.28820), (29.95275, 39.25112)]),
}
SCENE = "2015"                       # <- scene to render
cfg = SCENES[SCENE]
K, WATER_N, BORDER = 10, 2, 24
CLOUD_BR, CLOUD_SAT = 0.50, 0.09     # cloud = bright & low-saturation (neutral)

# short class names + CORINE Level-3 references
NAMES = {1: "Water", 2: "Forest", 3: "Dense vegetation", 4: "Shrub / transitional",
         5: "Irrigated crops", 6: "Cropland", 7: "Urban areas",
         8: "Agric. mosaic", 9: "Dry grass / fallow", 10: "Bare / sparse"}
EXPL = {1: "Sea and ocean, water bodies, water courses (CLC 5.1.1-5.2.3)",
        2: "Broad-leaved / mixed forest (CLC 3.1.1-3.1.3)",
        3: "Transitional woodland-shrub, dense (CLC 3.2.4)",
        4: "Natural grasslands & transitional shrub (CLC 3.2.1-3.2.4)",
        5: "Permanently irrigated land / complex cultivation (CLC 2.1.2-2.4.2)",
        6: "Non-irrigated arable land (CLC 2.1.1)",
        7: "Urban fabric, industrial & commercial units (CLC 1.1.1-1.2.1)",
        8: "Complex cultivation / land principally agriculture (CLC 2.4.2-2.4.3)",
        9: "Natural grasslands, sparsely vegetated (CLC 3.2.1-3.3.3)",
        10: "Sparsely vegetated / bare soil, harvested fields (CLC 3.3.3)"}

# --- classification helpers ------------------------------------------------- #
def cluster_water_land(V, swaps):
    """Water = 2 darkest clusters; 9 land clusters brightness-ordered; +swaps."""
    km0 = KMeans(K, random_state=0, n_init=10).fit(V)
    dark = np.argsort(km0.cluster_centers_.sum(1))[:WATER_N]
    isw = np.isin(km0.labels_, dark)
    kl = KMeans(K - 1, random_state=0, n_init=10).fit(V[~isw])
    lo = np.argsort(kl.cluster_centers_.sum(1))
    rm = {o: n for n, o in enumerate(lo)}
    fin = np.empty(V.shape[0], int); fin[isw] = 0
    fin[~isw] = np.vectorize(rm.get)(kl.labels_) + 1
    for a, bnum in swaps:
        va, vb = a - 1, bnum - 1
        ma, mb = fin == va, fin == vb
        fin[ma], fin[mb] = vb, va
    return fin

def read_valid(src, cloud):
    rgb = np.asarray(Image.open(src).convert("RGB"), np.float32) / 255.0
    Hh, Ww, _ = rgb.shape; fl = rgb.reshape(-1, 3)
    val = fl.sum(1) > (BORDER / 255.0); scn = val.reshape(Hh, Ww).copy()
    cl2 = np.zeros((Hh, Ww), bool)
    if cloud:
        br = fl.mean(1); sat = fl.max(1) - fl.min(1)
        c = val & (br > CLOUD_BR) & (sat < CLOUD_SAT)
        cl2 = ndimage.binary_dilation(c.reshape(Hh, Ww), iterations=2)
        val = val & ~cl2.ravel()
    return rgb, fl, val, scn, cl2, Hh, Ww

def reference_model(ref):
    """Classify the reference scene; return class prototypes + its (mu, sd)."""
    _, fl, val, _, _, _, _ = read_valid(ref["src"], ref["cloud"])
    fin = cluster_water_land(fl[val], ref.get("swaps", []))
    proto = np.array([fl[val][fin == c].mean(0) for c in range(K)])
    return proto, fl[val].mean(0), fl[val].std(0)

# --- read scene ------------------------------------------------------------- #
rgb = np.asarray(Image.open(cfg["src"]).convert("RGB"), np.float32) / 255.0
H, W, _ = rgb.shape
flat = rgb.reshape(-1, 3)
valid = flat.sum(1) > (BORDER / 255.0)       # non-border scene mask
scene2d = valid.reshape(H, W).copy()

# --- cloud mask (bright + neutral); dilate; exclude from clustering --------- #
cloud2d = np.zeros((H, W), bool)
cloud_frac = 0.0
if cfg["cloud"]:
    br = flat.mean(1); sat = flat.max(1) - flat.min(1)
    cl = valid & (br > CLOUD_BR) & (sat < CLOUD_SAT)
    cloud2d = ndimage.binary_dilation(cl.reshape(H, W), iterations=2)
    cloud_frac = 100.0 * cloud2d.sum() / valid.sum()
    valid = valid & ~cloud2d.ravel()

# --- classify: independent (reference) or matched to a reference scene ------- #
if cfg.get("match_to"):
    # reclassify this scene INTO the reference scene's classes, after linear
    # radiometric normalization to the reference (removes global brightness/
    # phenology offset -> consistent classes -> plausible change).
    proto, mu_ref, sd_ref = reference_model(SCENES[cfg["match_to"]])
    mu, sd = flat[valid].mean(0), flat[valid].std(0)
    norm = (flat - mu) / sd * sd_ref + mu_ref
    d = ((norm[:, None, :] - proto[None, :, :]) ** 2).sum(2)
    cls = d.argmin(1)
    lab = np.full(flat.shape[0], -1, int); lab[valid] = cls[valid]
else:
    fin = cluster_water_land(flat[valid], cfg.get("swaps", []))
    lab = np.full(flat.shape[0], -1, int); lab[valid] = fin
labels = lab.reshape(H, W)

# fill cloud-masked speckles with the nearest classified class (no white spots)
if cloud2d.any():
    known = labels >= 0
    ind = ndimage.distance_transform_edt(~known, return_distances=False, return_indices=True)
    nn = labels[tuple(ind)]
    labels[cloud2d] = nn[cloud2d]

tot = scene2d.sum()
pct = [100.0 * np.sum((labels == k) & scene2d) / tot for k in range(K)]
colors = plt.get_cmap("jet")(np.linspace(0, 1, K))[:, :3]

out = np.ones((H, W, 3), np.float32)
for k in range(K):
    out[labels == k] = colors[k]

# --- exact georeferencing from scene corner metadata ------------------------- #
CLL = cfg["corners"]                          # UL, UR, LL, LR (lon, lat)
CPX = [(0, 0), (W - 1, 0), (0, H - 1), (W - 1, H - 1)]
Tf = Transformer.from_crs("EPSG:4326", "EPSG:32635", always_xy=True)
EN = [Tf.transform(lon, lat) for lon, lat in CLL]
M = np.array([[c, r, 1] for c, r in CPX])
cE = np.linalg.lstsq(M, np.array([e for e, n in EN]), rcond=None)[0]
cN = np.linalg.lstsq(M, np.array([n for e, n in EN]), rcond=None)[0]
A = np.array([[cE[0], cE[1]], [cN[0], cN[1]]]); b = np.array([cE[2], cN[2]])
Ai = np.linalg.inv(A)

def ll2px(lon, lat):
    E, N = Tf.transform(lon, lat)
    return tuple(Ai @ (np.array([E, N]) - b))

# --- figure ----------------------------------------------------------------- #
fig = plt.figure(figsize=(12.2, 8.8))
ax = fig.add_axes([0.055, 0.085, 0.545, 0.83])
ax.imshow(out, interpolation="nearest")

lon_t = [27.5, 28.0, 28.5, 29.0, 29.5, 30.0]; lat_t = [39.5, 40.0, 40.5, 41.0]
la_s = np.linspace(39.25, 41.39, 60); lo_s = np.linspace(27.29, 30.04, 60)
for L in lon_t:
    xy = np.array([ll2px(L, a) for a in la_s]); ax.plot(xy[:, 0], xy[:, 1], color="0.45", lw=0.4, alpha=0.7)
for P in lat_t:
    xy = np.array([ll2px(o, P) for o in lo_s]); ax.plot(xy[:, 0], xy[:, 1], color="0.45", lw=0.4, alpha=0.7)
ax.set_xticks([ll2px(L, 39.28)[0] for L in lon_t]); ax.set_xticklabels([f"{L:.1f}\u00b0E" for L in lon_t], fontsize=8.5)
ax.set_yticks([ll2px(27.30, P)[1] for P in lat_t]); ax.set_yticklabels([f"{P:.1f}\u00b0N" for P in lat_t], fontsize=8.5)
ax.set_xlim(0, W); ax.set_ylim(H, 0)
ax.tick_params(length=3, pad=2)

halo_w = [pe.withStroke(linewidth=2.2, foreground="black")]
halo_b = [pe.withStroke(linewidth=2.4, foreground="white")]
sx, sy = ll2px(28.28, 40.71)
ax.text(sx, sy, "Sea of Marmara", color="white", fontsize=12, style="italic",
        fontweight="bold", ha="center", va="center", rotation=0, path_effects=halo_w)
ix, iy = ll2px(29.40, 40.43)
ax.text(ix, iy, "Iznik Lake", color="white", fontsize=9, fontweight="bold",
        ha="left", va="center", path_effects=halo_w)
cx, cy = ll2px(28.98, 41.01)
ax.plot(cx, cy, marker="s", color="red", markersize=11, markeredgecolor="white",
        markeredgewidth=1.4, path_effects=[pe.withStroke(linewidth=2.8, foreground="black")])
ax.text(cx + 16, cy, "Istanbul", color="black", fontsize=11, fontweight="bold",
        ha="left", va="center", path_effects=halo_b)

fig.text(0.328, 0.938, f"{cfg['sat']} OLI \u2014 {cfg['date_t']} \u2014 Path 180 / Row 032 (southern Marmara)",
         ha="center", fontsize=13, fontweight="bold")

# --- single merged legend block --------------------------------------------- #
import textwrap
axL = fig.add_axes([0, 0, 1, 1]); axL.axis("off"); axL.set_xlim(0, 1); axL.set_ylim(0, 1)
x_sq, x_tx = 0.615, 0.652
figH = fig.get_size_inches()[1]
lh = 8.4 * 1.42 / 72.0 / figH
y = 0.882
axL.text(x_sq, 0.912, "Land cover classes (CORINE Level-3)", fontsize=10.5,
         fontweight="bold", va="top", ha="left")
for k in range(K):
    txt = f"{k+1}. {NAMES[k+1]} \u2014 {EXPL[k+1]}  ({pct[k]:.1f}%)"
    wrapped = textwrap.fill(txt, width=34)
    nlines = wrapped.count("\n") + 1
    axL.add_patch(Rectangle((x_sq, y - lh * 0.95), 0.022, lh * 0.92,
                            facecolor=tuple(np.clip(colors[k], 0, 1)),
                            edgecolor="black", lw=0.5))
    axL.text(x_tx, y, wrapped, fontsize=8.4, va="top", ha="left", linespacing=1.35)
    y -= nlines * lh + 0.011

# bottom annotation (2 rows)
fig.text(0.055, 0.05,
         f"Data: USGS {cfg['sat']} OLI/TIRS C2, Scene ID: {cfg['scene_id']} ({cfg['date_a']}).",
         fontsize=9, color="0.25", ha="left")
meth = "Method: k-means clustering (GRASS GIS i.cluster, i.maxlik)"
meth += ", clouds masked. Source: authors." if cfg["cloud"] else ". Source: authors."
fig.text(0.055, 0.028, meth, fontsize=9, color="0.25", ha="left")

fig.savefig(cfg["out"], dpi=300, bbox_inches="tight", pad_inches=0.08)
print("saved", cfg["out"], "| cloud masked %.2f%%" % cloud_frac, "| legend ends y=%.3f" % y)
for k in range(K):
    print(f"  {k+1:2d} {NAMES[k+1]:20s} {pct[k]:5.1f}%")
