#!/usr/bin/env python3
"""
SEAMLESS MNDWI mosaic of the Sea of Marmara, 2025 (WRS-2 181/032 + 180/032).

Removes the visible scene border by:
  1. Radiometric matching  - fit a linear model  west ~= a*east + b  on the
     overlap pixels and apply it to the whole east scene, so the two scenes
     share one radiometric scale (kills the tonal step).
  2. Distance feathering    - blend across the overlap with a weight that ramps
     smoothly from all-west to all-east (kills the hard edge).

Memory-lean (3 GB box): feather weights are computed on a decimated grid and
upsampled; only three full-size float32 buffers are held at once.
"""
import numpy as np, rasterio, matplotlib as mpl, os
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from scipy.ndimage import distance_transform_edt, zoom

S, O, RES = 0.0000275, -0.2, 30.0
U = "/mnt/user-data/uploads/"
NAME = "mndwi"
SCENES = {   # name: (Green B3, SWIR1 B6)
    "west": (f"{U}LC08_L2SP_181032_20250722_20250730_02_T1_SR_B3.TIF",
             f"{U}LC08_L2SP_181032_20250722_20250730_02_T1_SR_B6.TIF"),
    "east": (f"{U}LC09_L2SP_180032_20250707_20250709_02_T1_SR_B3.TIF",
             f"{U}LC09_L2SP_180032_20250707_20250709_02_T1_SR_B6.TIF"),
}


def index_from_bands(p1, p2):
    with rasterio.open(p1) as d1:
        a_dn = d1.read(1); bounds = d1.bounds; crs = d1.crs
    with rasterio.open(p2) as d2:
        b_dn = d2.read(1)
    valid = (a_dn > 0) & (b_dn > 0)
    a = np.clip(a_dn.astype("float32") * S + O, 0, 1)
    b = np.clip(b_dn.astype("float32") * S + O, 0, 1)
    del a_dn, b_dn
    out = np.full(a.shape, np.nan, "float32")
    den = a + b
    ok = valid & (den > 0)
    out[ok] = (a[ok] - b[ok]) / den[ok]     # (B3-B6)/(B3+B6)
    return out, bounds, crs


# --- union grid --------------------------------------------------------------
meta = {}
for k, (p1, p2) in SCENES.items():
    with rasterio.open(p1) as ds:
        meta[k] = ds.bounds; crs = ds.crs
uL = min(b.left for b in meta.values());  uR = max(b.right for b in meta.values())
uT = max(b.top for b in meta.values());   uB = min(b.bottom for b in meta.values())
UW = int(round((uR - uL) / RES)); UH = int(round((uT - uB) / RES))
print(f"union {UW}x{UH}")

# --- place each scene onto union grid (west_u, east_u) -----------------------
def place(name):
    a, b, _ = index_from_bands(*SCENES[name])
    u = np.full((UH, UW), np.nan, "float32")
    coff = int(round((b.left - uL) / RES)); roff = int(round((uT - b.top) / RES))
    h, w = a.shape
    u[roff:roff + h, coff:coff + w] = a
    return u

west_u = place("west")
east_u = place("east")
has_w = np.isfinite(west_u); has_e = np.isfinite(east_u)
ov = has_w & has_e
print(f"overlap pixels: {int(ov.sum()):,}")

# --- 1. radiometric match  west ~= a*east + b  (fit on overlap subsample) -----
wi, ei = west_u[ov], east_u[ov]
rng = np.random.default_rng(0)
sub = rng.choice(wi.size, min(400_000, wi.size), replace=False)
a_coef, b_coef = np.polyfit(ei[sub], wi[sub], 1)
before = float(np.mean(wi - ei))
east_u = (a_coef * east_u + b_coef).astype("float32")     # apply to whole east
after = float(np.mean(west_u[ov] - east_u[ov]))
print(f"match: east' = {a_coef:.4f}*east + {b_coef:+.4f}  "
      f"| overlap mean diff {before:+.4f} -> {after:+.4f}")

# --- 2. distance-feather weight (decimated, then upsampled) -------------------
D = 8
wo = (has_w & ~has_e)[::D, ::D]        # west-only (small)
eo = (has_e & ~has_w)[::D, ::D]        # east-only (small)
dW = distance_transform_edt(~wo)       # dist to nearest west-only
dE = distance_transform_edt(~eo)       # dist to nearest east-only
w_small = (dE / (dW + dE + 1e-6)).astype("float32")   # ~1 near west, ~0 near east
del wo, eo, dW, dE
w_full = zoom(w_small, (UH / w_small.shape[0], UW / w_small.shape[1]),
              order=1).astype("float32")
w_full = w_full[:UH, :UW]
if w_full.shape != (UH, UW):           # pad safety
    tmp = np.zeros((UH, UW), "float32"); tmp[:w_full.shape[0], :w_full.shape[1]] = w_full
    w_full = tmp
del w_small

# --- blend (reuse west_u as output buffer) ----------------------------------
mosaic = west_u                        # alias
mosaic[ov] = (w_full[ov] * west_u[ov] + (1.0 - w_full[ov]) * east_u[ov]).astype("float32")
eonly = has_e & ~has_w
mosaic[eonly] = east_u[eonly]
del east_u, w_full
vm = mosaic[np.isfinite(mosaic)]
print(f"MOSAIC valid {vm.size:,}  min/mean/median/max "
      f"{vm.min():.3f}/{vm.mean():.3f}/{np.median(vm):.3f}/{vm.max():.3f}")

# --- write full-res + 100 m + PNG -------------------------------------------
T = from_origin(uL, uT, RES, RES)
prof = dict(driver="GTiff", height=UH, width=UW, count=1, dtype="float32",
            crs=crs, transform=T, nodata=np.nan, compress="deflate",
            predictor=3, tiled=True, blockxsize=256, blockysize=256)
full = f"/mnt/user-data/outputs/{NAME}_mosaic_2025_seamless.tif"
with rasterio.open(full, "w", **prof) as d:
    d.write(mosaic, 1); d.set_band_description(1, "MNDWI seamless mosaic 2025")
print("wrote", full, round(os.path.getsize(full) / 1e6, 1), "MB")

res2 = 100.0
nw = int((uR - uL) // res2); nh = int((uT - uB) // res2)
T2 = from_origin(uL, uT, res2, res2)
dst = np.full((nh, nw), np.nan, "float32")
reproject(mosaic, dst, src_transform=T, src_crs=crs, dst_transform=T2,
          dst_crs=crs, src_nodata=np.nan, dst_nodata=np.nan,
          resampling=Resampling.average)
scaled = np.where(np.isfinite(dst), np.round(dst * 10000), -32768).astype("int16")
p2 = dict(driver="GTiff", height=nh, width=nw, count=1, dtype="int16", crs=crs,
          transform=T2, nodata=-32768, compress="deflate", predictor=2,
          tiled=True, blockxsize=256, blockysize=256)
small = f"/mnt/user-data/outputs/{NAME}_mosaic_2025_seamless_100m.tif"
with rasterio.open(small, "w", **p2) as d:
    d.write(scaled, 1); d.set_band_description(1, "MNDWI*10000 int16 @100m seamless")
print("wrote", small, round(os.path.getsize(small) / 1e6, 2), "MB")

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["DejaVu Sans", "Arial"],
                     "axes.labelpad": 2, "axes.titlepad": 3})
disp = mosaic[::8, ::8]
ext = [uL / 1e3, uR / 1e3, uB / 1e3, uT / 1e3]
fig, ax = plt.subplots(figsize=(9.2, 6.6), constrained_layout=True)
cmap = plt.get_cmap("RdBu").copy(); cmap.set_bad("white")
im = ax.imshow(np.ma.masked_invalid(disp), extent=ext, origin="upper",
               cmap=cmap, vmin=-0.7, vmax=0.7, interpolation="nearest")
ax.set_xlabel("UTM 35N easting (km)"); ax.set_ylabel("UTM 35N northing (km)")
ax.set_title("MNDWI seamless mosaic - Sea of Marmara 2025 "
             "(181/032 + 180/032, matched + feathered)",
             fontsize=11, fontweight="bold")
ax.xaxis.set_minor_locator(AutoMinorLocator()); ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.grid(which="major", lw=0.4, alpha=0.35); ax.grid(which="minor", lw=0.25, alpha=0.2)
ax.tick_params(which="both", top=True, right=True, direction="out")
ax.set_aspect("equal")
cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02, extend="both")
cb.set_label("MNDWI  (water > 0, land < 0)", fontsize=9.5); cb.ax.tick_params(labelsize=8)
ax.text(0.995, -0.12,
        "Data: USGS Landsat 8/9 OLI C2 L2 SR (B3,B6), 2025-07-22 & 2025-07-07. "
        f"Radiometric match (east'={a_coef:.3f}*east{b_coef:+.3f}) + distance feather. "
        f"Software: rasterio {rasterio.__version__}, scipy, matplotlib. Source: authors.",
        transform=ax.transAxes, ha="right", va="top", fontsize=6.5, color="0.35")
for e in ("png", "pdf"):
    fig.savefig(f"/mnt/user-data/outputs/fig_{NAME}_mosaic_2025_seamless.{e}",
                bbox_inches="tight", pad_inches=0.03, dpi=600 if e == "png" else None)
print("wrote figure")
