#!/usr/bin/env python3
"""
Shared figure style for the Marmara index figures so every figure matches:
  - font family : Nimbus Sans (URW Helvetica = LaTeX 'helvet'; the TeX Helvetica)
  - title       : bold, font size auto-fitted so it spans the figure width
                  (map left edge -> colorbar right edge), centred over that span,
                  placed above the map with guaranteed clearance
  - credit      : 8 pt, word-wrapped to exactly the map (axes) width

Usage:
    from figstyle import apply_style, add_fitted_title, add_wrapped_credit
    apply_style()
    ... build fig, ax, colorbar cb ...
    add_fitted_title(fig, ax, cb.ax, "My title")
    add_wrapped_credit(fig, ax, "Data: ... Map: authors.")
"""
import os, textwrap
import matplotlib as mpl
import matplotlib.font_manager as fm

_HERE = os.path.dirname(os.path.abspath(__file__))
CREDIT_PT = 8.0


def apply_style():
    """Register Nimbus Sans (fallback DejaVu Sans) and set global rcParams."""
    fam = "DejaVu Sans"
    for fn in ("NimbusSans-Regular.otf", "NimbusSans-Bold.otf",
               "NimbusSans-Italic.otf", "NimbusSans-BoldItalic.otf"):
        p = os.path.join(_HERE, fn)
        if os.path.exists(p):
            fm.fontManager.addfont(p)
    if any("Nimbus Sans" in f.name for f in fm.fontManager.ttflist):
        fam = "Nimbus Sans"
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [fam, "DejaVu Sans", "Arial"],
        "axes.labelpad": 2, "axes.titlepad": 3,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
    })
    return fam


def _axspan_fig_fraction(fig, ax, right_ax):
    """Left edge of `ax` to right edge of `right_ax`, in figure-fraction x."""
    fig.canvas.draw()
    l = ax.get_position().x0
    r = right_ax.get_position().x1
    return l, r


def add_fitted_title(fig, ax, right_ax, text, max_pt=13.0, min_pt=7.0, pad=0.012):
    """Bold title whose font size is shrunk until it spans <= the map->colorbar
    width, centred over that span, sitting just above the map top edge."""
    fig.canvas.draw()
    l, r = _axspan_fig_fraction(fig, ax, right_ax)
    span_px = (r - l) * fig.get_figwidth() * fig.dpi
    xc = 0.5 * (l + r)
    ytop = ax.get_position().y1
    t = fig.text(xc, ytop + pad, text, ha="center", va="bottom",
                 fontweight="bold", fontsize=max_pt)
    rend = fig.canvas.get_renderer()
    pt = max_pt
    while pt > min_pt:
        w = t.get_window_extent(renderer=rend).width
        if w <= span_px:
            break
        pt -= 0.25
        t.set_fontsize(pt)
    return t


def add_wrapped_credit(fig, ax, text, pt=CREDIT_PT, gap=0.055):
    """8 pt credit, word-wrapped to the map (axes) width, left-aligned under it."""
    fig.canvas.draw()
    pos = ax.get_position()
    width_px = pos.width * fig.get_figwidth() * fig.dpi
    rend = fig.canvas.get_renderer()
    # binary-ish search on wrap width (chars) so the longest line ~= map width
    probe = fig.text(0, 0, "x", fontsize=pt)
    char_px = probe.get_window_extent(renderer=rend).width
    probe.remove()
    ncols = max(20, int(width_px / max(char_px, 1e-6)))
    lines = textwrap.wrap(text, width=ncols)
    # shrink columns until the widest wrapped line fits the map width
    for _ in range(40):
        tmp = fig.text(0, 0, "\n".join(lines), fontsize=pt)
        w = tmp.get_window_extent(renderer=rend).width
        tmp.remove()
        if w <= width_px or ncols <= 20:
            break
        ncols -= 2
        lines = textwrap.wrap(text, width=ncols)
    c = fig.text(pos.x0, pos.y0 - gap, "\n".join(lines),
                 ha="left", va="top", fontsize=pt, linespacing=1.25)
    return c


# --- GRASS GIS "ndvi" colour table (authoritative value->RGB, domain [-1,1]) --
# Source: OSGeo/grass lib/gis/colors/ndvi. Absolute NDVI values, so render with
# vmin=-1, vmax=1 to keep colours pinned to NDVI values (good for comparisons).
_GRASS_NDVI = [
    (-1.000, (5, 24, 82)),   (-0.300, (5, 24, 82)),   (-0.180, (255, 255, 255)),
    (0.000, (255, 255, 255)),(0.025, (206, 197, 180)),(0.075, (191, 163, 124)),
    (0.125, (179, 174, 96)), (0.150, (163, 181, 80)), (0.175, (144, 170, 60)),
    (0.233, (166, 195, 29)), (0.266, (135, 183, 3)),  (0.333, (121, 175, 1)),
    (0.366, (101, 163, 0)),  (0.433, (78, 151, 0)),   (0.466, (43, 132, 4)),
    (0.550, (0, 114, 0)),    (0.650, (0, 90, 1)),     (0.750, (0, 73, 0)),
    (0.850, (0, 56, 0)),     (0.950, (0, 31, 0)),     (1.000, (0, 0, 0)),
]


def ndvi_cmap():
    """matplotlib colormap reproducing the GRASS 'ndvi' cpt over NDVI [-1, 1]."""
    from matplotlib.colors import LinearSegmentedColormap
    stops = [((v + 1) / 2.0, (r/255, g/255, b/255)) for v, (r, g, b) in _GRASS_NDVI]
    return LinearSegmentedColormap.from_list("ndvi", stops)


def byg_cmap():
    """GRASS 'byg' colour table: blue (0%) -> yellow (50%) -> green (100%)."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(
        "byg", [(0.0, (0, 0, 1)), (0.5, (1, 1, 0)), (1.0, (0, 1, 0))])
