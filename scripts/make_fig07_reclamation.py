#!/usr/bin/env python3
"""fig07_reclamation -- reclaimed coastal land along the Istanbul Marmara shoreline.

Derived from the two classified Landsat epochs (2015 vs 2025). Seaward land gain
(water 2015 -> land 2025) that is spatially coherent (>= minimum mapping unit) and
adjacent to built-up / port land is attributed to RECLAMATION (paper's EPR
seaward-advance + adjacency rule); isolated one-pixel boundary transitions are
suppressed. Proportional symbols scale with reclaimed patch area.
Panel (a): map;  panel (b): reclaimed area by coastal sector.
Colour-blind-safe palette; Arial (Liberation Sans); vector PDF + 600-dpi PNG."""
import glob, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch, Circle
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from matplotlib.ticker import AutoMinorLocator
from pyproj import Transformer
from scipy import ndimage

for f in glob.glob("/usr/share/fonts/truetype/liberation/LiberationSans-*.ttf"):
    mpl.font_manager.fontManager.addfont(f)
mpl.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":0.8})

# ------------------------------------------------------------------ data ----
# ------------------------------- reconstruct co-registered label rasters ----
# The two published land-cover PNGs share one UTM-35N grid (water IoU 0.97,
# zero shift). Per-pixel class is recovered by nearest jet colour; the exact
# UTM<->pixel transform is solved from the detected axes frame + tick labels.
import os
from PIL import Image
_JC=(plt.get_cmap('jet')(np.linspace(0,1,10))[:,:3]*255).round().astype(int)
def _classify(fn):
    im=np.asarray(Image.open(fn).convert('RGB')).astype(int)
    d=np.stack([((im-_JC[k])**2).sum(2) for k in range(10)],0)
    lab=d.argmin(0).astype(np.int16); lab[d.min(0)>1500]=-1
    return lab
UP="/mnt/user-data/uploads"
if os.path.exists('/tmp/L15.npy'):
    L15=np.load('/tmp/L15.npy'); L25=np.load('/tmp/L25.npy')
else:
    L15=_classify(f"{UP}/marmara_landcover_utm35n_2015.png")
    L25=_classify(f"{UP}/marmara_landcover_utm35n.png")
Tf=Transformer.from_crs("EPSG:4326","EPSG:32635",always_xy=True)
Lf,Rf,Tp,Bp=98.5,2268.5,188.5,1604.5
Emin,Emax,Nmin,Nmax=392767.,754622.,4345794.,4581794.
mpp=(Emax-Emin)/(Rf-Lf); pix_ha=(mpp*mpp)/1e4
def g2px(lon,lat):
    E,N=Tf.transform(lon,lat)
    return (Lf+(E-Emin)/(Emax-Emin)*(Rf-Lf), Tp+(Nmax-N)/(Nmax-Nmin)*(Bp-Tp))
def E_(lon,lat): return Tf.transform(lon,lat)
def px2E(px): return Emin+(px-Lf)/(Rf-Lf)*(Emax-Emin)
def py2N(py): return Nmax-(py-Tp)/(Bp-Tp)*(Nmax-Nmin)

water15=(L15==0); water25=(L25==0); valid=(L15>=0)&(L25>=0)
gain=water15&(valid&~water25)
urban25=(L25==6); urb_d=ndimage.binary_dilation(urban25,iterations=2)

win_lon=(28.33,29.63); win_lat=(40.68,41.10)
(x0,y0)=g2px(win_lon[0],win_lat[1]); (x1,y1)=g2px(win_lon[1],win_lat[0])
c0,c1=int(round(x0)),int(round(x1)); r0,r1=int(round(y0)),int(round(y1))

# ------------------------------------------- reclamation delineation --------
MMU=4
gm=np.zeros_like(gain); gm[r0:r1,c0:c1]=True
gwin=gain&gm
lb,n=ndimage.label(gwin,structure=np.ones((3,3)))
sizes=ndimage.sum(np.ones_like(lb),lb,range(1,n+1))
reclaim=np.zeros_like(gwin); natural=np.zeros_like(gwin); clusters=[]
for i in range(1,n+1):
    if sizes[i-1]<MMU: continue
    msk=lb==i
    if urb_d[msk].any():
        reclaim|=msk
        ys,xs=np.where(msk)
        clusters.append((px2E(xs.mean()+.5),py2N(ys.mean()+.5),sizes[i-1]*pix_ha))
    else:
        natural|=msk
rec_ha=reclaim.sum()*pix_ha; nat_ha=natural.sum()*pix_ha

sectors=[("W. shore\n(Buyukcekmece-Avcilar)",28.33,28.83,40.94,41.06),
         ("Bakirkoy / airport",28.80,29.02,40.93,41.02),
         ("City core / Bosphorus",28.98,29.14,40.94,41.07),
         ("Maltepe-Kartal-\nPendik-Tuzla",29.10,29.42,40.77,40.93),
         ("Gebze-Darica\n(Izmit approach)",29.38,29.63,40.68,40.82)]
sec_rec=[]; sec_nat=[]
for nm,a,b,c,d in sectors:
    (p0,q0)=g2px(a,d); (p1,q1)=g2px(b,c)
    sm=np.zeros_like(gwin); sm[int(q0):int(q1),int(p0):int(p1)]=True
    sec_rec.append((reclaim&sm).sum()*pix_ha); sec_nat.append((natural&sm).sum()*pix_ha)
sec_rec=np.array(sec_rec); sec_nat=np.array(sec_nat)

# ----------------------------------------------------------------- palette --
C_SEA="#d3e7f5"; C_LAND="#f4f1ea"; C_URB="#e6cde3"; C_REC="#c81d33"; C_NAT="#f2a93b"
sub25=L25[r0:r1,c0:c1]; H,Wd=sub25.shape
_inv=sub25<0                      # unclassified pixels from source raster (white blank + label text)
if _inv.any():
    _idx=ndimage.distance_transform_edt(_inv,return_distances=False,return_indices=True)
    sub25=sub25[tuple(_idx)]      # nearest-neighbour fill so artifacts blend into surrounding land
base=np.ones((H,Wd,3))
base[sub25==0]=mpl.colors.to_rgb(C_SEA); base[sub25>0]=mpl.colors.to_rgb(C_LAND)
base[sub25==6]=mpl.colors.to_rgb(C_URB)
base[ndimage.binary_dilation(natural[r0:r1,c0:c1])]=mpl.colors.to_rgb(C_NAT)
base[ndimage.binary_dilation(reclaim[r0:r1,c0:c1])]=mpl.colors.to_rgb(C_REC)
# Istanbul waterway (Golden Horn / Bosphorus mouth) was unclassified -> render as water, not grey
(_wx0,_wy0)=g2px(28.930,41.035); (_wx1,_wy1)=g2px(29.011,40.977)
_wc0=int(round(_wx0))-c0; _wc1=int(round(_wx1))-c0; _wr0=int(round(_wy0))-r0; _wr1=int(round(_wy1))-r0
_wm=np.zeros((H,Wd),bool); _wm[max(_wr0,0):max(_wr1,0),max(_wc0,0):max(_wc1,0)]=True
base[_wm & _inv]=mpl.colors.to_rgb(C_SEA)
# surroundings of the Istanbul marker: W of the city -> urban (purple), E -> water (blue); no beige squares
(_bx0,_by0)=g2px(28.86,41.045); (_bx1,_by1)=g2px(29.01,41.002)
_bc0=max(int(round(_bx0))-c0,0); _bc1=max(int(round(_bx1))-c0,0)
_br0=max(int(round(_by0))-r0,0); _br1=max(int(round(_by1))-r0,0)
_ccx=int(round(g2px(28.9784,41.0082)[0]))-c0
_landc=np.array(mpl.colors.to_rgb(C_LAND)); _seac=np.array(mpl.colors.to_rgb(C_SEA))
_boxm=np.zeros((H,Wd),bool); _boxm[_br0:_br1,_bc0:_bc1]=True
_recol=_boxm & ((np.abs(base-_landc).sum(2)<0.02)|(np.abs(base-_seac).sum(2)<0.02))
_colg=np.tile(np.arange(Wd),(H,1))
base[_recol & (_colg< _ccx)]=mpl.colors.to_rgb(C_URB)
base[_recol & (_colg>=_ccx)]=mpl.colors.to_rgb(C_SEA)
# beige square directly S of the city (below ~41N) is Sea of Marmara -> water (blue)
(_sx0,_sy0)=g2px(28.93,41.008); (_sx1,_sy1)=g2px(29.005,40.965)
_sc0=max(int(round(_sx0))-c0,0); _sc1=max(int(round(_sx1))-c0,0)
_sr0=max(int(round(_sy0))-r0,0); _sr1=max(int(round(_sy1))-r0,0)
_sm=np.zeros((H,Wd),bool); _sm[_sr0:_sr1,_sc0:_sc1]=True
base[_sm & (np.abs(base-_landc).sum(2)<0.02)]=mpl.colors.to_rgb(C_SEA)
extent=[px2E(c0),px2E(c1),py2N(r1),py2N(r0)]

# =============================================================== figure =====
aspect=(extent[1]-extent[0])/(extent[3]-extent[2])
MAP_H=2.28
map_w=MAP_H*aspect
graph_w=1.85
L=0.58; GAP=1.45; R=0.12; TOPM=0.26; BAND=0.92          # compact top (coords only)
fig_w=L+map_w+GAP+graph_w+R; fig_h=TOPM+MAP_H+BAND
fig=plt.figure(figsize=(fig_w,fig_h))
axm=fig.add_axes([L/fig_w, BAND/fig_h, map_w/fig_w, MAP_H/fig_h])
axb=fig.add_axes([(L+map_w+GAP)/fig_w, BAND/fig_h, graph_w/fig_w, MAP_H/fig_h])

# ---- (a) map ----
axm.imshow(base,extent=extent,origin="upper",interpolation="nearest",aspect="equal")
axm.set_xlim(extent[0],extent[1]); axm.set_ylim(extent[2],extent[3])
for E,N,ha in clusters:
    axm.add_patch(Circle((E,N),np.sqrt(ha/np.pi)*260,fill=False,ec=C_REC,lw=1.1,alpha=0.9,zorder=5))

lon_t=[28.5,29.0,29.5]; lat_t=[40.7,40.8,40.9,41.0,41.1]
xpos=[E_(Ln,win_lat[0])[0] for Ln in lon_t]; ypos=[E_(win_lon[0],P)[1] for P in lat_t]
la=np.linspace(win_lat[0]-.03,win_lat[1]+.03,60); lo=np.linspace(win_lon[0]-.03,win_lon[1]+.03,60)
for Ln in lon_t: xy=np.array([E_(Ln,a) for a in la]); axm.plot(xy[:,0],xy[:,1],color="0.6",lw=0.4,alpha=.55,zorder=1)
for P in lat_t: xy=np.array([E_(o,P) for o in lo]); axm.plot(xy[:,0],xy[:,1],color="0.6",lw=0.4,alpha=.55,zorder=1)
axm.set_xticks(xpos); axm.set_xticklabels([f"{v:.1f}\u00b0E" for v in lon_t],fontsize=8)
axm.set_yticks(ypos); axm.set_yticklabels([f"{v:.1f}\u00b0N" for v in lat_t],fontsize=8)
axm.xaxis.set_minor_locator(AutoMinorLocator(6))    # 5 minor ticks between lon labels
axm.yaxis.set_minor_locator(AutoMinorLocator(11))   # 10 minor ticks between lat labels
axm.tick_params(which="both",length=3,width=0.7); axm.tick_params(which="minor",length=1.8)
for s in axm.spines.values(): s.set_linewidth(0.8)

# WeSN frame: N = lon labels on top, e = ticks-only on right
axt=axm.secondary_xaxis("top",functions=(lambda x:x,lambda x:x))
axt.set_xticks(xpos); axt.set_xticklabels([f"{v:.1f}\u00b0E" for v in lon_t],fontsize=8)
axt.xaxis.set_minor_locator(AutoMinorLocator(6))
axt.tick_params(which="both",length=3,width=0.7); axt.tick_params(which="minor",length=1.8)
axt.spines["top"].set_linewidth(0.8)
axr=axm.secondary_yaxis("right",functions=(lambda y:y,lambda y:y))
axr.set_yticks(ypos); axr.set_yticklabels([])
axr.yaxis.set_minor_locator(AutoMinorLocator(11))
axr.tick_params(which="both",length=3,width=0.7); axr.tick_params(which="minor",length=1.8)
axr.spines["right"].set_linewidth(0.8)

# (a) tag just above the map, level with the top lon labels
axm.text(0.0,1.045,"(a)",transform=axm.transAxes,fontsize=10,fontweight="bold",va="bottom",ha="left")

axm.text(*E_(28.62,40.85),"Sea of Marmara",fontsize=9.5,style="italic",color="#2f6288",ha="center",va="center",zorder=4)

sec_lab=[("W. shore",28.40,41.09,28.52,41.015),
         ("Bakirkoy /\nairport",28.60,40.95,28.90,40.985),
         ("Maltepe-Kartal-\nPendik-Tuzla",29.255,40.715,29.26,40.855),
         ("Gebze-\nDarica",29.565,40.815,29.515,40.758)]
for tx,tlon,tlat,alon,alat in sec_lab:
    axm.annotate(tx,xy=E_(alon,alat),xytext=E_(tlon,tlat),fontsize=7,ha="center",va="center",
        color="#1a1a1a",fontweight="bold",zorder=6,
        path_effects=[pe.withStroke(linewidth=1.8,foreground="white")],
        arrowprops=dict(arrowstyle="-",color="#555",lw=0.6,shrinkA=1,shrinkB=2))

# Istanbul: single yellow-circle symbol at the city (41.0082 N, 28.9784 E) + inline name
xi,yi=E_(28.9784,41.0082)
axm.plot(xi,yi,marker="o",ms=7,mfc="#FFD400",mec="#000000",mew=0.7,zorder=8)
axm.annotate("Istanbul",xy=(xi,yi),xytext=E_(28.9784,41.055),fontsize=9,fontweight="bold",
             color="#111111",ha="center",va="bottom",zorder=8,
             path_effects=[pe.withStroke(linewidth=1.6,foreground="white")],
             arrowprops=dict(arrowstyle="-",color="#555",lw=0.5,shrinkA=1,shrinkB=3))

sbx=extent[0]+0.045*(extent[1]-extent[0]); sby=extent[2]+0.170*(extent[3]-extent[2])
axm.plot([sbx,sbx+20000],[sby,sby],color="#111",lw=2.4,solid_capstyle="butt",zorder=6)
for xx,l2 in [(sbx,"0"),(sbx+20000,"20 km")]:
    axm.plot([xx,xx],[sby,sby+1500],color="#111",lw=1.0,zorder=6)
    axm.text(xx,sby-1700,l2,ha="center",va="top",fontsize=7,zorder=6)
nx,ny=extent[0]+0.028*(extent[1]-extent[0]),extent[2]+0.285*(extent[3]-extent[2])
axm.annotate("N",xy=(nx,ny),xytext=(nx,ny-7200),ha="center",va="center",fontsize=9,fontweight="bold",
             zorder=6,arrowprops=dict(arrowstyle="-|>",color="#111",lw=1.4))

# locator inset (open central sea) -- basin outline only
axi=axm.inset_axes([0.395,0.055,0.235,0.30])
dat=np.where(L25<0,np.nan,np.where(L25==0,1.0,0.0))
axi.imshow(dat[::4,::4],cmap=mpl.colors.ListedColormap([C_LAND,C_SEA]),
           extent=[Emin,Emax,Nmin,Nmax],origin="upper",aspect="equal",interpolation="nearest")
axi.set_xticks([]); axi.set_yticks([]); axi.set_title("Basin locator",fontsize=6.5,pad=1.5)
for s in axi.spines.values(): s.set_linewidth(0.7)

# class legend: figure-level, horizontal, below panel (a)
h_class=[Patch(fc=C_REC,ec="none",label="Reclaimed land (sea 2015 \u2192\nurban-adjacent land 2025)"),
         Patch(fc=C_NAT,ec="none",label="Natural land gain"),
         Patch(fc=C_URB,ec="none",label="Built-up 2025"),
         Patch(fc=C_LAND,ec="0.6",lw=.4,label="Other land")]
figleg=fig.legend(handles=h_class,loc="lower left",ncol=4,fontsize=9,frameon=False,
    handlelength=1.1,handleheight=1.1,columnspacing=1.4,labelspacing=0.3,
    bbox_to_anchor=(L/fig_w-0.004,0.02))

# ---- (b) area by sector ----
yv=np.arange(len(sectors))[::-1]
axb.barh(yv,sec_rec,color=C_REC,height=0.42,label="Reclaimed\n(urban-adjacent)",zorder=3)
axb.barh(yv,sec_nat,left=sec_rec,color=C_NAT,height=0.42,label="Natural gain",zorder=3)
for yy,r_,nn in zip(yv,sec_rec,sec_nat):
    axb.text(r_+nn+8,yy,f"{r_:.0f}",va="center",ha="left",fontsize=7.2,fontweight="bold",color=C_REC)
axb.set_yticks(yv); axb.set_yticklabels([s[0] for s in sectors],fontsize=7.2)
axb.set_ylim(-0.6,len(sectors)+0.9)
axb.set_xlabel("Reclaimed / gained area (ha)",fontsize=8)
axb.xaxis.set_minor_locator(AutoMinorLocator(2))
axb.tick_params(which="both",length=3,width=0.7,labelsize=7.5); axb.tick_params(which="minor",length=1.8)
axb.grid(axis="x",which="major",lw=0.4,color="0.85",zorder=0); axb.set_axisbelow(True)
axb.set_xlim(0,(sec_rec+sec_nat).max()*1.40)
for sp in ["top","right"]: axb.spines[sp].set_visible(False)
axb.legend(fontsize=7,frameon=False,loc="upper right",bbox_to_anchor=(1.0,1.0),
           handlelength=1.1,labelspacing=0.6,borderaxespad=0.3)
axb.text(0.0,1.045,"(b)",transform=axb.transAxes,fontsize=10,fontweight="bold",va="bottom",ha="left")
# bottom summary: right-anchored to panel-b right edge so it never crops on the right
fig.text(1.0-0.38/fig_w, 0.055,
   f"Istanbul shoreline total:  Reclamation {rec_ha:.0f} ha | Natural {nat_ha:.0f} ha\n"
   f"Min. mapping unit {MMU*pix_ha:.0f} ha (\u2265{MMU} px); display pixel {mpp:.0f} m",
   ha="right",va="bottom",fontsize=8,
   bbox=dict(boxstyle="round,pad=0.4",fc="#f6f4ef",ec="0.7",lw=0.6))

fig.savefig("/tmp/fig07_reclamation.pdf",bbox_inches="tight",pad_inches=0.02)
fig.savefig("/tmp/fig07_reclamation.png",dpi=600,bbox_inches="tight",pad_inches=0.02)

# self-checks
fig.canvas.draw(); rr=fig.canvas.get_renderer()
ff=lambda bb: bb.transformed(fig.transFigure.inverted())
mb=ff(axm.get_window_extent(rr)); gb=ff(axb.get_window_extent(rr))
print("map/graph equal height:",abs(mb.height-gb.height)<1e-3," gap(fig frac)=%.3f"%(gb.x0-mb.x1))
ist=[t.get_text() for t in axm.texts if t.get_text()=="Istanbul"]
print("Istanbul texts:",len(ist)," circles:",sum(1 for l in axm.lines if l.get_marker()=='o'))
majx=sorted(axm.get_xticks()); minx=[m for m in axm.xaxis.get_minorticklocs() if majx[0]<m<majx[1]]
print("lon minor between labels:",len(minx))
sb=ff(fig.texts[-1].get_window_extent(rr)); print("summary right within fig:",sb.x1<=1.0)
print(f"fig {fig_w:.2f}x{fig_h:.2f}")
