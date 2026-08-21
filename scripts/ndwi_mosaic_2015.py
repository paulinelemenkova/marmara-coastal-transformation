#!/usr/bin/env python3
"""
NDWI mosaic of the Sea of Marmara, 2015 (WRS-2 181/032 + 180/032, Landsat 8),
from Collection-2 Level-2 surface reflectance B3 (Green) and B5 (NIR):

  WEST  WRS-2 181/032  Landsat 8  2015-07-27
  EAST  WRS-2 180/032  Landsat 8  2015-07-20

C2 L2 scaling: reflectance = DN*0.0000275 - 0.2 ; fill = 0 ; clip to [0,1].
NDWI = (Green - NIR)/(Green + NIR)  [McFeeters 1996]:  water > 0, land < 0.
Integer-offset mosaic + linear east->west match on the overlap + distance
feather. Masking DN>0. Colour ramp: jet. Style: shared figstyle.
"""
import os, sys, gc
import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from scipy.ndimage import distance_transform_edt
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle

U, OUT = "/mnt/user-data/uploads/", "/mnt/user-data/outputs/"
YEAR, NAME = "2015", "ndwi"
S, O, RES = 0.0000275, -0.2, 30.0
VMIN, VMAX = -0.6, 0.6
SCENES = {  # name: (Green B3, NIR B5)  ->  NDWI = (B3-B5)/(B3+B5)
    "west": (f"{U}LC08_L2SP_181032_20150727_20200908_02_T1_SR_B3.TIF",
             f"{U}LC08_L2SP_181032_20150727_20200908_02_T1_SR_B5.TIF"),
    "east": (f"{U}LC08_L2SP_180032_20150720_20200908_02_T1_SR_B3.TIF",
             f"{U}LC08_L2SP_180032_20150720_20200908_02_T1_SR_B5.TIF"),
}


def index_from_bands(pa, pb):
    with rasterio.open(pa) as d:
        a_dn = d.read(1); bnds = d.bounds; crs = d.crs
    with rasterio.open(pb) as d:
        b_dn = d.read(1)
    valid = (a_dn > 0) & (b_dn > 0)
    a = np.clip(a_dn.astype("float32") * S + O, 0, 1)
    b = np.clip(b_dn.astype("float32") * S + O, 0, 1)
    del a_dn, b_dn
    out = np.full(a.shape, np.nan, "float32")
    den = a + b
    ok = valid & (den > 0)
    out[ok] = (a[ok] - b[ok]) / den[ok]
    return out, bnds, crs


idx, bnd = {}, {}
crs = None
for k, (pa, pb) in SCENES.items():
    idx[k], bnd[k], crs = index_from_bands(pa, pb)
    print(f"{k}: {idx[k].shape}  valid {np.isfinite(idx[k]).mean()*100:.1f}%")

uL = min(b.left for b in bnd.values()); uR = max(b.right for b in bnd.values())
uT = max(b.top for b in bnd.values());  uB = min(b.bottom for b in bnd.values())
Wpx = int(round((uR - uL) / RES)); Hpx = int(round((uT - uB) / RES))
T = from_origin(uL, uT, RES, RES)
print(f"union grid: {Wpx} x {Hpx}")


def off(b):
    return int(round((uT - b.top) / RES)), int(round((b.left - uL) / RES))


mos = np.full((Hpx, Wpx), np.nan, "float32")
rW, cW = off(bnd["west"]); hW, wW = idx["west"].shape
m = np.isfinite(idx["west"]); mos[rW:rW + hW, cW:cW + wW][m] = idx["west"][m]; del m
rE, cE = off(bnd["east"]); hE, wE = idx["east"].shape

r0, r1 = max(rW, rE), min(rW + hW, rE + hE)
c0, c1 = max(cW, cE), min(cW + wW, cE + wE)
Wob = mos[r0:r1, c0:c1].copy()
Eov = idx["east"][r0 - rE:r1 - rE, c0 - cE:c1 - cE]
both = np.isfinite(Wob) & np.isfinite(Eov)
xe = Eov[both].astype("float64"); yw = Wob[both].astype("float64")
a, b = np.polyfit(xe, yw, 1)
print(f"overlap pixels: {both.sum():,}  match west=a*east+b: a={a:.4f} b={b:+.4f}")
print(f"overlap mean(west-east): {np.mean(yw-xe):+.4f} -> after match {np.mean(yw-(a*xe+b)):+.4f}")
del Eov, both, xe, yw; gc.collect()

east_m = a * idx["east"] + b
east_m[~np.isfinite(idx["east"])] = np.nan
del idx; gc.collect()

eb = mos[rE:rE + hE, cE:cE + wE]
fill = np.isfinite(east_m) & ~np.isfinite(eb); eb[fill] = east_m[fill]; del fill

Eob = east_m[r0 - rE:r1 - rE, c0 - cE:c1 - cE]
del east_m; gc.collect()
wcov = np.isfinite(Wob); ecov = np.isfinite(Eob); both_ov = wcov & ecov
dW = distance_transform_edt(wcov); dE = distance_transform_edt(ecov)
with np.errstate(invalid="ignore", divide="ignore"):
    wgt = (dW / (dW + dE)).astype("float32")
del dW, dE
blend = np.where(both_ov, wgt * np.nan_to_num(Wob) + (1 - wgt) * np.nan_to_num(Eob),
                 np.where(wcov, Wob, Eob)).astype("float32")
mos[r0:r1, c0:c1] = blend
del Wob, Eob, wcov, ecov, both_ov, wgt, blend; gc.collect()

v = mos[np.isfinite(mos)]
print(f"mosaic valid {v.size:,}  NDWI mean/median {v.mean():.3f}/{np.median(v):.3f}  water(>0) {100*np.mean(v>0):.1f}%")

with rasterio.open(f"{OUT}{NAME}_mosaic_{YEAR}_30m.tif", "w", driver="GTiff",
                   height=Hpx, width=Wpx, count=1, dtype="float32", crs=crs, transform=T,
                   nodata=np.nan, compress="deflate", predictor=3, tiled=True,
                   blockxsize=256, blockysize=256) as dst:
    dst.write(mos, 1); dst.set_band_description(1, "NDWI mosaic 2015 (181/032+180/032)")
print("wrote 30m")

FACT = 100.0 / RES
H2 = int(np.ceil(Hpx / FACT)); W2 = int(np.ceil(Wpx / FACT))
T100 = from_origin(uL, uT, 100.0, 100.0)
m100 = np.full((H2, W2), np.nan, "float32")
reproject(mos, m100, src_transform=T, src_crs=crs, dst_transform=T100, dst_crs=crs,
          resampling=Resampling.average)
i16 = np.where(np.isfinite(m100), np.clip(np.round(m100 * 1e4), -32767, 32767), -32768).astype("int16")
with rasterio.open(f"{OUT}{NAME}_mosaic_{YEAR}_100m.tif", "w", driver="GTiff",
                   height=H2, width=W2, count=1, dtype="int16", crs=crs, transform=T100,
                   nodata=-32768, compress="deflate", predictor=2, tiled=True) as dst:
    dst.write(i16, 1); dst.set_band_description(1, "NDWI x1e4 (Int16), 100 m")
print("wrote 100m", round(os.path.getsize(f'{OUT}{NAME}_mosaic_{YEAR}_100m.tif')/1e6,1), "MB")

# ---- quick-look figure ------------------------------------------------------
figstyle.apply_style()
DEC = 5
disp = mos[::DEC, ::DEC]; dx = RES * DEC
ext = [uL/1e3, (uL + disp.shape[1]*dx)/1e3, (uT - disp.shape[0]*dx)/1e3, uT/1e3]
fig, ax = plt.subplots(figsize=(8.4, 6.6))
im = ax.imshow(np.ma.masked_invalid(disp), extent=ext, origin="upper",
               cmap="jet", vmin=VMIN, vmax=VMAX, interpolation="nearest")
ax.set_xlabel("Easting (km, UTM 35N)"); ax.set_ylabel("Northing (km, UTM 35N)")
ax.set_aspect("equal")
ax.xaxis.set_minor_locator(AutoMinorLocator()); ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which="major", length=4, labelsize=8); ax.tick_params(which="minor", length=2)
ax.grid(which="major", color="white", lw=0.3, alpha=0.5)
cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02, extend="both")
cb.set_label("NDWI  (water > 0, land < 0)", fontsize=9); cb.ax.tick_params(labelsize=8)
cb.locator = plt.MultipleLocator(0.2); cb.update_ticks()
figstyle.add_fitted_title(fig, ax, cb.ax,
    "NDWI mosaic - Sea of Marmara 2015 (WRS-2 181/032 + 180/032, Landsat 8)")
figstyle.add_wrapped_credit(fig, ax,
    "Landsat 8 OLI Collection-2 Level-2 surface reflectance: WRS-2 181/032 "
    "(2015-07-27) and 180/032 (2015-07-20). UTM zone 35N / EPSG:32635, 30 m, "
    "identical pixel grid (integer-offset mosaic). NDWI = (B3-B5)/(B3+B5) "
    "[McFeeters 1996] on C2 L2 reflectance (DN x 0.0000275 - 0.2, clipped to [0,1]); "
    "fill (DN=0) masked; east scene radiometrically matched to west on the overlap "
    "and distance-blended. Colour: jet. Data courtesy USGS/NASA. Map: authors.")
for e in ("png", "pdf"):
    fig.savefig(f"{OUT}fig_{NAME}_mosaic_{YEAR}.{e}", bbox_inches="tight",
                pad_inches=0.04, dpi=300 if e == "png" else None)
print("wrote figure")
