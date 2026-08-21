#!/usr/bin/env python3
"""
Two-panel MNDWI comparison, Sea of Marmara:
    (a) 2015  (WRS-2 181/032 + 180/032, Landsat 8)
    (b) 2025  (WRS-2 181/032 + 180/032, Landsat 8/9)
Both panels share one RdBu colour scale [-0.7, 0.7] (water > 0, land < 0) and a
single common colorbar. Both are drawn on the same geographic window (the union
of the two mosaics) so they are spatially registered and directly comparable.
Concise title fitted so its text spans the width of the two map panels.

Inputs: the two MNDWI mosaic GeoTIFFs (100 m is plenty for a 2-panel figure).
"""
import os, sys, textwrap
import numpy as np
import rasterio
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle

OUT = "/mnt/user-data/outputs/"
TIF = {
    "2015": f"{OUT}mndwi_mosaic_2015_100m.tif",
    "2025": "/mnt/user-data/uploads/mndwi_mosaic_2025_100m.tif",
}
VMIN, VMAX, CMAP = -0.7, 0.7, "RdBu"


def load(path):
    with rasterio.open(path) as ds:
        a = ds.read(1).astype("float32")
        nd = ds.nodata
        if nd is not None:
            a[a == nd] = np.nan
        if np.nanmax(np.abs(a)) > 5:      # Int16 stored as MNDWI*1e4
            a = a / 1e4
        T = ds.transform
        L, T0 = T.c, T.f
        R = L + ds.width * T.a
        B = T0 + ds.height * T.e
    return a, [L/1e3, R/1e3, B/1e3, T0/1e3]   # array, extent in km


A15, e15 = load(TIF["2015"])
A25, e25 = load(TIF["2025"])

# common geographic window (union), km
xL = min(e15[0], e25[0]); xR = max(e15[1], e25[1])
yB = min(e15[2], e25[2]); yT = max(e15[3], e25[3])

figstyle.apply_style()
fig, (axa, axb) = plt.subplots(1, 2, figsize=(12.6, 5.6))
fig.subplots_adjust(left=0.06, right=0.88, bottom=0.13, top=0.86, wspace=0.10)

for ax, arr, ext, tag, yr in [(axa, A15, e15, "a", "2015"),
                              (axb, A25, e25, "b", "2025")]:
    im = ax.imshow(np.ma.masked_invalid(arr), extent=ext, origin="upper",
                   cmap=CMAP, vmin=VMIN, vmax=VMAX, interpolation="nearest")
    ax.set_xlim(xL, xR); ax.set_ylim(yB, yT)
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (km, UTM 35N)")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="major", length=4, labelsize=8)
    ax.tick_params(which="minor", length=2)
    ax.grid(which="major", color="0.5", lw=0.3, alpha=0.5)
    ax.set_title(f"({tag}) {yr}", fontsize=10, fontweight="bold", pad=4)

axa.set_ylabel("Northing (km, UTM 35N)")
axb.tick_params(labelleft=False)

# one common colorbar for both panels
cax = fig.add_axes([0.90, 0.13, 0.016, 0.73])
cb = fig.colorbar(im, cax=cax, extend="both")
cb.set_label("MNDWI  (water > 0, land < 0)", fontsize=9)
cb.ax.tick_params(labelsize=8)
cb.locator = plt.MultipleLocator(0.2); cb.update_ticks()


def fitted_title_over_panels(fig, ax_left, ax_right, text, max_pt=15.0, min_pt=8.0):
    """Bold title centred over, and fitted to, the two-panel span."""
    fig.canvas.draw()
    l = ax_left.get_position().x0
    r = ax_right.get_position().x1
    span_px = (r - l) * fig.get_figwidth() * fig.dpi
    ytop = max(ax_left.get_position().y1, ax_right.get_position().y1)
    t = fig.text(0.5*(l+r), ytop + 0.045, text, ha="center", va="bottom",
                 fontweight="bold", fontsize=max_pt)
    rend = fig.canvas.get_renderer(); pt = max_pt
    while pt > min_pt and t.get_window_extent(renderer=rend).width > span_px:
        pt -= 0.25; t.set_fontsize(pt)
    return t


def wrapped_credit_over_panels(fig, ax_left, ax_right, text, pt=8.0):
    fig.canvas.draw()
    l = ax_left.get_position().x0
    r = ax_right.get_position().x1
    width_px = (r - l) * fig.get_figwidth() * fig.dpi
    rend = fig.canvas.get_renderer()
    probe = fig.text(0, 0, "x", fontsize=pt); cw = probe.get_window_extent(renderer=rend).width
    probe.remove()
    ncols = max(40, int(width_px / max(cw, 1e-6)))
    lines = textwrap.wrap(text, width=ncols)
    for _ in range(40):
        tmp = fig.text(0, 0, "\n".join(lines), fontsize=pt)
        w = tmp.get_window_extent(renderer=rend).width; tmp.remove()
        if w <= width_px or ncols <= 40:
            break
        ncols -= 2; lines = textwrap.wrap(text, width=ncols)
    return fig.text(l, 0.055, "\n".join(lines), ha="left", va="top",
                    fontsize=pt, linespacing=1.25)


fitted_title_over_panels(fig, axa, axb,
    "MNDWI mosaic - Sea of Marmara 2015-2025")
wrapped_credit_over_panels(fig, axa, axb,
    "Landsat OLI Collection-2 Level-2 surface reflectance, WRS-2 181/032 + 180/032, "
    "UTM 35N / EPSG:32635. (a) 2015 Landsat 8 (2015-07-27 / 2015-07-20); "
    "(b) 2025 Landsat 8/9 (2025-07-22 / 2025-07-07). MNDWI = (B3-B6)/(B3+B6); "
    "east scene radiometrically matched to west on the overlap and distance-blended; "
    "shared RdBu scale [-0.7, 0.7]. Data courtesy USGS/NASA. Map: authors.")

for e in ("png", "pdf"):
    fig.savefig(f"{OUT}fig_mndwi_2panel_2015_2025.{e}",
                bbox_inches="tight", pad_inches=0.04, dpi=300 if e == "png" else None)
print("wrote fig_mndwi_2panel_2015_2025.png/.pdf")
