#!/usr/bin/env python3
"""
NDVI mosaic of the Sea of Marmara, 2025, from Landsat Collection-2 Level-2
surface reflectance (bands B4=Red, B5=NIR):

  WEST  WRS-2 181/032  Landsat 8  2025-07-22
  EAST  WRS-2 180/032  Landsat 9  2025-07-07

C2 L2 scaling: reflectance = DN*0.0000275 - 0.2 ; fill = 0.
NDVI = (NIR - Red)/(NIR + Red), computed on reflectance.
Both scenes share EPSG:32635, 30 m, and the same pixel grid (origin = 15 mod 30),
so the mosaic is placed by integer offsets and averaged in the overlap.
"""
import numpy as np, rasterio, matplotlib as mpl, os
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

S, O, RES = 0.0000275, -0.2, 30.0
U = "/mnt/user-data/uploads/"
SCENES = {
    "west_181032": (f"{U}LC08_L2SP_181032_20250722_20250730_02_T1_SR_B4.TIF",
                    f"{U}LC08_L2SP_181032_20250722_20250730_02_T1_SR_B5.TIF"),
    "east_180032": (f"{U}LC09_L2SP_180032_20250707_20250709_02_T1_SR_B4.TIF",
                    f"{U}LC09_L2SP_180032_20250707_20250709_02_T1_SR_B5.TIF"),
}


def ndvi_from_bands(red_path, nir_path):
    with rasterio.open(red_path) as dr:
        red_dn = dr.read(1); bounds = dr.bounds; crs = dr.crs
    with rasterio.open(nir_path) as dn_:
        nir_dn = dn_.read(1)
    valid = (red_dn > 0) & (nir_dn > 0)                     # exclude fill only
    red = np.clip(red_dn.astype("float32") * S + O, 0, 1)
    nir = np.clip(nir_dn.astype("float32") * S + O, 0, 1)
    del red_dn, nir_dn
    ndvi = np.full(red.shape, np.nan, "float32")
    den = nir + red
    ok = valid & (den > 0)
    ndvi[ok] = (nir[ok] - red[ok]) / den[ok]
    del red, nir, den
    return ndvi, bounds, crs


# --- compute per-scene NDVI --------------------------------------------------
scene_arrays = {}
for name, (r, n) in SCENES.items():
    a, b, crs = ndvi_from_bands(r, n)
    scene_arrays[name] = (a, b)
    v = a[np.isfinite(a)]
    print(f"{name}: {np.isfinite(a).sum():,} valid px, NDVI mean {v.mean():.3f}")

# --- union grid (integer-aligned) -------------------------------------------
lefts = [b.left for _, b in scene_arrays.values()]
rights = [b.right for _, b in scene_arrays.values()]
tops = [b.top for _, b in scene_arrays.values()]
bots = [b.bottom for _, b in scene_arrays.values()]
uL, uR, uT, uB = min(lefts), max(rights), max(tops), min(bots)
UW = int(round((uR - uL) / RES)); UH = int(round((uT - uB) / RES))
print(f"union grid {UW} x {UH}  ({uL:.0f},{uB:.0f})-({uR:.0f},{uT:.0f})")

ssum = np.zeros((UH, UW), "float32")
cnt = np.zeros((UH, UW), "uint8")
for name, (a, b) in scene_arrays.items():
    coff = int(round((b.left - uL) / RES))
    roff = int(round((uT - b.top) / RES))
    h, w = a.shape
    sl = (slice(roff, roff + h), slice(coff, coff + w))
    m = np.isfinite(a)
    ssum[sl][m] += a[m]
    cnt[sl][m] += 1
overlap_px = int((cnt >= 2).sum())
mosaic = np.where(cnt > 0, ssum / np.maximum(cnt, 1), np.nan).astype("float32")
del ssum
print(f"overlap (averaged) pixels: {overlap_px:,}")
vm = mosaic[np.isfinite(mosaic)]
print(f"MOSAIC valid {vm.size:,}  NDVI min/mean/median/max "
      f"{vm.min():.3f}/{vm.mean():.3f}/{np.median(vm):.3f}/{vm.max():.3f}")

from rasterio.transform import from_origin
T = from_origin(uL, uT, RES, RES)
prof = dict(driver="GTiff", height=UH, width=UW, count=1, dtype="float32",
            crs=crs, transform=T, nodata=np.nan, compress="deflate",
            predictor=3, tiled=True, blockxsize=256, blockysize=256)
full = "/mnt/user-data/outputs/ndvi_mosaic_2025.tif"
with rasterio.open(full, "w", **prof) as d:
    d.write(mosaic, 1); d.set_band_description(1, "NDVI mosaic 2025 (181/032+180/032)")
print("wrote", full, round(os.path.getsize(full) / 1e6, 1), "MB")

# --- downloadable 100 m Int16 version ---------------------------------------
from rasterio.warp import reproject, Resampling
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
small = "/mnt/user-data/outputs/ndvi_mosaic_2025_100m.tif"
with rasterio.open(small, "w", **p2) as d:
    d.write(scaled, 1); d.set_band_description(1, "NDVI*10000 int16 @100m")
print("wrote", small, round(os.path.getsize(small) / 1e6, 2), "MB")

# --- quick-look PNG ----------------------------------------------------------
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["DejaVu Sans", "Arial"],
                     "axes.labelpad": 2, "axes.titlepad": 3})
dec = 8
disp = mosaic[::dec, ::dec]
ext = [uL / 1e3, uR / 1e3, uB / 1e3, uT / 1e3]
fig, ax = plt.subplots(figsize=(9.2, 6.6), constrained_layout=True)
cmap = plt.get_cmap("viridis").copy(); cmap.set_bad("white")
im = ax.imshow(np.ma.masked_invalid(disp), extent=ext, origin="upper",
               cmap=cmap, vmin=-0.1, vmax=0.8, interpolation="nearest")
ax.set_xlabel("UTM 35N easting (km)"); ax.set_ylabel("UTM 35N northing (km)")
ax.set_title("NDVI mosaic - Sea of Marmara 2025 "
             "(WRS-2 181/032 + 180/032, Landsat 8/9)",
             fontsize=11, fontweight="bold")
ax.xaxis.set_minor_locator(AutoMinorLocator()); ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.grid(which="major", lw=0.4, alpha=0.35); ax.grid(which="minor", lw=0.25, alpha=0.2)
ax.tick_params(which="both", top=True, right=True, direction="out")
ax.set_aspect("equal")
cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02, extend="both")
cb.set_label("NDVI", fontsize=10); cb.ax.tick_params(labelsize=8)
ax.text(0.995, -0.12,
        "Data: USGS Landsat 8/9 OLI C2 L2 SR (B4,B5), 2025-07-22 & 2025-07-07. "
        f"Software: rasterio {rasterio.__version__}/GDAL {rasterio.__gdal_version__}, "
        f"matplotlib {mpl.__version__}. Overlap averaged. Source: authors.",
        transform=ax.transAxes, ha="right", va="top", fontsize=7, color="0.35")
for e in ("png", "pdf"):
    fig.savefig(f"/mnt/user-data/outputs/fig_ndvi_mosaic_2025.{e}",
                bbox_inches="tight", pad_inches=0.03, dpi=600 if e == "png" else None)
print("wrote quick-look figure")
