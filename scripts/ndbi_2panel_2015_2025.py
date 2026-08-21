#!/usr/bin/env python3
"""
Two-panel NDBI comparison, Sea of Marmara, stacked vertically:
    (a) 2015  (WRS-2 181/032 + 180/032, Landsat 8)      [top]
    (b) 2025  (WRS-2 181/032 + 180/032, Landsat 8/9)    [bottom]
Shared GRASS 'ryg' scale [-0.5, 0.5] (0 at yellow): built-up/bare (NDBI>0) green,
vegetation (NDBI<0) red. Single common colorbar spanning both panels.
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
    "2015": f"{OUT}ndbi_mosaic_2015_100m.tif",
    "2025": "/mnt/user-data/uploads/ndbi_mosaic_2025_100m.tif",
}
VMIN, VMAX = -0.5, 0.5
CMAP = figstyle.ryg_cmap()


def load(path):
    with rasterio.open(path) as ds:
        a = ds.read(1).astype("float32")
        if ds.nodata is not None:
            a[a == ds.nodata] = np.nan
        if np.nanmax(np.abs(a)) > 5:
            a = a / 1e4
        T = ds.transform
        ext = [T.c/1e3, (T.c + ds.width*T.a)/1e3,
               (T.f + ds.height*T.e)/1e3, T.f/1e3]
    return a, ext


A15, e15 = load(TIF["2015"])
A25, e25 = load(TIF["2025"])
xL = min(e15[0], e25[0]); xR = max(e15[1], e25[1])
yB = min(e15[2], e25[2]); yT = max(e15[3], e25[3])

figstyle.apply_style()
fig, (axa, axb) = plt.subplots(2, 1, figsize=(9.4, 10.2))
fig.subplots_adjust(left=0.09, right=0.86, bottom=0.10, top=0.90, hspace=0.14)

for ax, arr, ext, tag, yr in [(axa, A15, e15, "a", "2015"),
                              (axb, A25, e25, "b", "2025")]:
    im = ax.imshow(np.ma.masked_invalid(arr), extent=ext, origin="upper",
                   cmap=CMAP, vmin=VMIN, vmax=VMAX, interpolation="nearest")
    ax.set_xlim(xL, xR); ax.set_ylim(yB, yT)
    ax.set_aspect("equal")
    ax.set_ylabel("Northing (km, UTM 35N)")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="major", length=4, labelsize=8)
    ax.tick_params(which="minor", length=2)
    ax.grid(which="major", color="0.5", lw=0.3, alpha=0.5)
    ax.set_title(f"({tag}) {yr}", fontsize=10, fontweight="bold", pad=4)

axb.set_xlabel("Easting (km, UTM 35N)")
axa.tick_params(labelbottom=False)

p_top = axa.get_position(); p_bot = axb.get_position()
cax = fig.add_axes([0.88, p_bot.y0, 0.016, p_top.y1 - p_bot.y0])
cb = fig.colorbar(im, cax=cax, extend="both")
cb.set_label("NDBI  (built-up/bare > 0, veg < 0)", fontsize=9)
cb.ax.tick_params(labelsize=8)
cb.locator = plt.MultipleLocator(0.25); cb.update_ticks()


def fitted_title_over_panels(fig, ax_left, ax_right, text, max_pt=15.0, min_pt=8.0):
    fig.canvas.draw()
    l = ax_left.get_position().x0; r = ax_right.get_position().x1
    span_px = (r - l) * fig.get_figwidth() * fig.dpi
    ytop = max(ax_left.get_position().y1, ax_right.get_position().y1)
    t = fig.text(0.5*(l+r), ytop + 0.02, text, ha="center", va="bottom",
                 fontweight="bold", fontsize=max_pt)
    rend = fig.canvas.get_renderer(); pt = max_pt
    while pt > min_pt and t.get_window_extent(renderer=rend).width > span_px:
        pt -= 0.25; t.set_fontsize(pt)
    return t


def wrapped_credit_over_panels(fig, ax_left, right_edge, text, pt=8.0):
    fig.canvas.draw()
    l = ax_left.get_position().x0
    r = right_edge if isinstance(right_edge, (int, float)) else right_edge.get_position().x1
    width_px = (r - l) * fig.get_figwidth() * fig.dpi
    rend = fig.canvas.get_renderer()
    probe = fig.text(0, 0, "x", fontsize=pt); cw = probe.get_window_extent(renderer=rend).width
    probe.remove()
    ncols = max(40, int(width_px / max(cw, 1e-6)))
    lines = textwrap.wrap(text, width=ncols)
    for _ in range(60):
        tmp = fig.text(0, 0, "\n".join(lines), fontsize=pt)
        w = tmp.get_window_extent(renderer=rend).width; tmp.remove()
        if w <= width_px or ncols <= 40:
            break
        ncols -= 2; lines = textwrap.wrap(text, width=ncols)
    return fig.text(l, 0.055, "\n".join(lines), ha="left", va="top",
                    fontsize=pt, linespacing=1.3)


fitted_title_over_panels(fig, axa, axb, "NDBI mosaic - Sea of Marmara 2015-2025")
wrapped_credit_over_panels(fig, axa, 0.965,
    "Landsat OLI C2 L2 surface reflectance, WRS-2 181/032 + 180/032, UTM 35N. "
    "(a) 2015 Landsat 8; (b) 2025 Landsat 8/9. NDBI = (B6-B5)/(B6+B5), east matched "
    "and feathered. Colour: GRASS 'ryg' [-0.5, 0.5] (built-up green, veg red). "
    "USGS/NASA; map: authors.")

for e in ("png", "pdf"):
    fig.savefig(f"{OUT}fig_ndbi_2panel_2015_2025.{e}",
                bbox_inches="tight", pad_inches=0.04, dpi=300 if e == "png" else None)
print("wrote fig_ndbi_2panel_2015_2025.png/.pdf")
