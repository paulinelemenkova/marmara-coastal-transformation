#!/usr/bin/env python3
"""
mosaic_map.py — georeferenced UTM 35N mosaic map of the Marmara Landsat scenes
(WRS-2 180/032 + 181/032, 2025) with a white lat/lon graticule and a cartographic
metadata block on the black free space. Font: Arial (Liberation Sans equivalent).
Input: /tmp/mos_utm.npz (out RGB, valid mask, UTM extent) built by the warp step.
Output: marmara_mosaic_utm35n.png
"""
import glob
import numpy as np
from pyproj import Transformer
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrow

# --- Arial (Liberation Sans = free metric-compatible equivalent) ------------ #
for f in glob.glob("/usr/share/fonts/truetype/liberation/LiberationSans-*.ttf"):
    mpl.font_manager.fontManager.addfont(f)
mpl.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"]})

d = np.load("/tmp/mos_utm.npz")
out = d["out"]; Emin, Emax, Nmin, Nmax = d["ext"]
Tf = Transformer.from_crs("EPSG:4326", "EPSG:32635", always_xy=True)
def ll2utm(lon, lat): return Tf.transform(lon, lat)

WHITE = "white"
fig = plt.figure(figsize=(14.5, 9.6), facecolor="black")
ax = fig.add_axes([0.055, 0.11, 0.90, 0.80]); ax.set_facecolor("black")
ax.imshow(out, extent=[Emin, Emax, Nmin, Nmax], origin="upper", interpolation="nearest")
ax.set_xlim(Emin, Emax); ax.set_ylim(Nmin, Nmax); ax.set_aspect("equal")

# --- white lat/lon graticule ------------------------------------------------ #
lon_t = [26.0, 26.5, 27.0, 27.5, 28.0, 28.5, 29.0, 29.5, 30.0]
lat_t = [39.5, 40.0, 40.5, 41.0]
la_s = np.linspace(39.25, 41.40, 80); lo_s = np.linspace(25.7, 30.05, 80)
for L in lon_t:
    xy = np.array([ll2utm(L, a) for a in la_s]); ax.plot(xy[:, 0], xy[:, 1], color=WHITE, lw=0.5, alpha=0.9)
for P in lat_t:
    xy = np.array([ll2utm(o, P) for o in lo_s]); ax.plot(xy[:, 0], xy[:, 1], color=WHITE, lw=0.5, alpha=0.9)
# labels on black margins
for L in lon_t:
    x, y = ll2utm(L, 39.28)
    if Emin < x < Emax:
        ax.annotate(f"{L:.1f}\u00b0E", (x, Nmin), xytext=(0, -6), textcoords="offset points",
                    ha="center", va="top", color=WHITE, fontsize=8.5, clip_on=False)
for P in lat_t:
    x, y = ll2utm(25.72, P)
    if Nmin < y < Nmax:
        ax.annotate(f"{P:.1f}\u00b0N", (Emin, y), xytext=(-6, 0), textcoords="offset points",
                    ha="right", va="center", color=WHITE, fontsize=8.5, clip_on=False)
for s in ax.spines.values():
    s.set_color(WHITE); s.set_linewidth(0.6)
ax.tick_params(colors="none", length=0)
ax.set_xticks([]); ax.set_yticks([])

# --- title (top black strip) ------------------------------------------------ #
fig.text(0.50, 0.965, "Sea of Marmara \u2014 Landsat mosaic (WRS-2 180/032 + 181/032)",
         ha="center", color=WHITE, fontsize=16, fontweight="bold")
fig.text(0.50, 0.928, "UTM Zone 35N / WGS 84 (EPSG:32635)", ha="center", color=WHITE, fontsize=10)

# --- metadata block (bottom-left black space) ------------------------------- #
meta = ("Scenes (Collection-2 Level-1, 2025):\n"
        "  West  \u2014 Path 181 / Row 032 \u2014 Landsat 8 OLI/TIRS \u2014 2025-07-22\n"
        "           Scene ID  LC81810322025203LGN00\n"
        "  East  \u2014 Path 180 / Row 032 \u2014 Landsat 9 OLI/TIRS \u2014 2025-07-07\n"
        "           Scene ID  LC91800322025188LGN00\n"
        "Projection: UTM 35N / WGS 84    Pixel: 250 m (resampled)\n"
        "Data: USGS Landsat.  Mosaic & cartography: authors.")
fig.text(0.062, 0.075, meta, ha="left", va="top", color=WHITE, fontsize=8.6,
         family="sans-serif", linespacing=1.5)

# --- scale bar (bottom-right black strip, figure coords) -------------------- #
km = 50.0
frac = km * 1000.0 / (Emax - Emin) * 0.90        # 0.90 = axes width fraction
x0, y0 = 0.74, 0.055
fig.add_artist(mpl.lines.Line2D([x0, x0 + frac], [y0, y0], color=WHITE, lw=2.2))
for xx, lab in [(x0, "0"), (x0 + frac / 2, "25"), (x0 + frac, "50 km")]:
    fig.add_artist(mpl.lines.Line2D([xx, xx], [y0, y0 + 0.008], color=WHITE, lw=1.2))
    fig.text(xx, y0 - 0.018, lab, ha="center", va="top", color=WHITE, fontsize=8)

# --- north arrow (top-right black wedge, data coords) ----------------------- #
nx, ny = 0.945, 0.78
ax.annotate("N", xy=(nx, ny + 0.07), xytext=(nx, ny), xycoords="axes fraction",
            ha="center", va="center", color=WHITE, fontsize=13, fontweight="bold",
            arrowprops=dict(arrowstyle="-|>", color=WHITE, lw=1.8))

# --- place names (added) ---------------------------------------------------- #
halo_dark = [pe.withStroke(linewidth=2.0, foreground="black")]   # for white text
halo_light = [pe.withStroke(linewidth=1.8, foreground="white")]  # for black text

# water bodies (white font)
def _wlabel(lon, lat, s, fs, ha="center"):
    x, y = ll2utm(lon, lat)
    ax.text(x, y, s, color="white", fontsize=fs, fontweight="bold", style="italic",
            ha=ha, va="center", path_effects=halo_dark, zorder=6)
_wlabel(28.05, 40.75, "Sea of Marmara", 15)
_wlabel(29.52, 40.43, "Iznik Lake", 10)

# Istanbul: red square symbol + black name
xi, yi = ll2utm(28.978, 41.008)
ax.plot(xi, yi, marker="s", color="red", markersize=9, markeredgecolor="white",
        markeredgewidth=1.0, zorder=7,
        path_effects=[pe.withStroke(linewidth=2.2, foreground="black")])
ax.text(xi + 3200, yi, "Istanbul", color="black", fontsize=10, fontweight="bold",
        ha="left", va="center", path_effects=halo_light, zorder=7)

# cities under the Landsat footprints: yellow circles + black names
_cities = {
    "Tekirda\u011f": (27.511, 40.978), "\u00c7orlu": (27.803, 41.159),
    "Gebze": (29.431, 40.803), "\u0130zmit": (29.918, 40.765),
    "Yalova": (29.277, 40.655), "Band\u0131rma": (27.976, 40.352),
    "Bursa": (29.061, 40.183), "\u0130neg\u00f6l": (29.512, 40.078),
    "Gelibolu": (26.670, 40.410), "Biga": (27.243, 40.228),
    "\u00c7anakkale": (26.414, 40.155), "\u00c7an": (27.053, 40.032),
    "G\u00f6nen": (27.653, 40.103), "Mustafakemalpa\u015fa": (28.410, 40.038),
    "Bal\u0131kesir": (27.887, 39.649), "Edremit": (27.024, 39.596),
    "\u0130znik": (29.720, 40.429),
}
for _name, (_lon, _lat) in _cities.items():
    _x, _y = ll2utm(_lon, _lat)
    ax.plot(_x, _y, marker="o", color="yellow", markersize=4.5, markeredgecolor="black",
            markeredgewidth=0.5, zorder=6)
    _left = _lon >= 29.2                    # east cities: label to the left
    ax.text(_x + (-2600 if _left else 2600), _y, _name, color="black", fontsize=7.5,
            ha="right" if _left else "left", va="center",
            path_effects=halo_light, zorder=6)

fig.savefig("marmara_mosaic_utm35n.png", dpi=200, facecolor="black",
            bbox_inches="tight", pad_inches=0.05)
print("saved marmara_mosaic_utm35n.png")
