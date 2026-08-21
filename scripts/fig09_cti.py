#!/usr/bin/env python3
"""fig09_cti -- composite coastal-transformation index (CTI) for the Marmara littoral.

CTI = 0.40*R' + 0.35*U' + 0.25*M' (Eq. CTI), each component min-max normalised to [0,1]
across the four coastal sub-regions:
  R' reclamation intensity  (water 2015 -> urban-adjacent land 2025, per coastal-zone area);
  U' urban/built-up density (2025 built-up per coastal-zone area);
  M' mucilage susceptibility (mean MSI of the sub-region's characteristic coastal water,
     from fig08).  Screening product: 167 m grid, no bathymetry (see fig07/fig08 caveats)."""
import glob, numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.patheffects as pe
from matplotlib.ticker import AutoMinorLocator
from pyproj import Transformer
from scipy import ndimage

for f in glob.glob("/usr/share/fonts/truetype/liberation/LiberationSans-*.ttf"):
    mpl.font_manager.fontManager.addfont(f)
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":0.8})
FS_PANEL, FS_MID, FS_SM = 11, 9, 8

# ---- MSI + sea from fig08, land-cover epochs ----
ns={}; exec(open('/tmp/fig08_mucilage.py').read().split('# ==')[0], ns)
MSI=ns['MSI']; sea=ns['sea']
L25=np.load('/tmp/L25.npy'); L15=np.load('/tmp/L15.npy')
mpp=166.8; pkm2=(mpp*mpp)/1e6
Tf=Transformer.from_crs("EPSG:4326","EPSG:32635",always_xy=True)
Ti=Transformer.from_crs("EPSG:32635","EPSG:4326",always_xy=True)
Lf,Rf,Tp,Bp=98.5,2268.5,188.5,1604.5
Emin,Emax,Nmin,Nmax=392767.,754622.,4345794.,4581794.
def g2pxf(lo,la):
    E,N=Tf.transform(lo,la); return (Lf+(E-Emin)/(Emax-Emin)*(Rf-Lf),Tp+(Nmax-N)/(Nmax-Nmin)*(Bp-Tp))
def g2px(lo,la):
    x,y=g2pxf(lo,la); return int(round(x)),int(round(y))
def E_(lo,la): return Tf.transform(lo,la)
def px2E(px): return Emin+(px-Lf)/(Rf-Lf)*(Emax-Emin)
def py2N(py): return Nmax-(py-Tp)/(Bp-Tp)*(Nmax-Nmin)

# ---- CTI components per sub-region ----
H,W=L25.shape; ys,xs=np.mgrid[0:H,0:W]
subs=[("Istanbul coast",28.95,41.00,(28.30,29.30,40.90,41.10)),
      ("Gulf of Izmit",29.68,40.75,(29.35,30.00,40.66,40.83)),
      ("Southern Marmara",28.25,40.40,(28.35,29.20,40.34,40.55)),
      ("Tekirdag/Canakkale",27.05,40.75,(26.55,27.95,40.55,41.02))]
names=[s[0] for s in subs]
cpx=[g2pxf(s[1],s[2]) for s in subs]
region=np.stack([(xs-c[0])**2+(ys-c[1])**2 for c in cpx]).argmin(0)
dist_sea=ndimage.distance_transform_edt(~sea)*mpp/1000.0
coastland=(L25>0)&(dist_sea<6)
valid=(L15>=0)&(L25>=0); gain=(L15==0)&valid&(L25>0)
lb,n=ndimage.label(gain,np.ones((3,3))); sz=ndimage.sum(np.ones_like(lb),lb,range(1,n+1))
reclaim=np.isin(lb,[i+1 for i in range(n) if sz[i]>=4])&ndimage.binary_dilation(L25==6,iterations=2)
urb25=(L25==6)
R=[];U=[];M=[]
for i,(nm,clo,cla,(a,b,c,d)) in enumerate(subs):
    cz=(region==i)&coastland; nz=max(cz.sum(),1)
    R.append(reclaim[cz].sum()/nz); U.append(urb25[cz].sum()/nz)
    (p0,q0)=g2px(a,d);(p1,q1)=g2px(b,c)
    box=np.zeros_like(sea); box[q0:q1,p0:p1]=True; ms=box&sea
    M.append(float(np.nanmean(MSI[ms])) if ms.sum() else 0.0)
R,U,M=map(np.array,(R,U,M))
mm=lambda x:(x-x.min())/(x.max()-x.min()+1e-9)
Rp,Up,Mp=mm(R),mm(U),mm(M)
wr,wu,wm=0.40*Rp,0.35*Up,0.25*Mp
CTI=wr+wu+wm
print("CTI:",{names[i]:round(float(CTI[i]),2) for i in range(4)})

# ---- basemap: sea coloured by sub-region CTI ----
win=[26.45,30.02,40.20,41.22]
(c0,r0)=g2px(win[0],win[3]); (c1,r1)=g2px(win[1],win[2])
sub=lambda A:A[r0:r1,c0:c1]
Ls,seaS,regS=sub(L25),sub(sea),sub(region)
cmap=plt.get_cmap("plasma"); norm=Normalize(0,1)
Hs,Ws=Ls.shape; rgb=np.empty((Hs,Ws,3))
rgb[:]=mpl.colors.to_rgb("#dddcd6")                       # monochrome land
ctiimg=np.choose(regS,CTI)                                # CTI per pixel by region
rgb[seaS]=cmap(norm(ctiimg[seaS]))[:,:3]
extent=[px2E(c0),px2E(c1),py2N(r1),py2N(r0)]

# =============================================================== figure =====
aspect=(extent[1]-extent[0])/(extent[3]-extent[2])
MAP_H=2.62; map_w=MAP_H*aspect; graph_w=2.15
L=0.62; GAP=0.95; Rm=0.14; TOPM=0.42; BAND=1.40
fig_w=L+map_w+GAP+graph_w+Rm; fig_h=TOPM+MAP_H+BAND
fig=plt.figure(figsize=(fig_w,fig_h))
axm=fig.add_axes([L/fig_w,BAND/fig_h,map_w/fig_w,MAP_H/fig_h])
axb=fig.add_axes([(L+map_w+GAP)/fig_w,BAND/fig_h,graph_w/fig_w,MAP_H/fig_h])

# (a) map
axm.imshow(rgb,extent=extent,origin="upper",interpolation="nearest",aspect="equal")
axm.set_xlim(extent[0],extent[1]); axm.set_ylim(extent[2],extent[3])
lon_t=[27,28,29]; lat_t=[40.4,40.6,40.8,41.0,41.2]
la=np.linspace(win[2]-.05,win[3]+.05,80); lo=np.linspace(win[0]-.05,win[1]+.05,80)
for Ln in lon_t: xy=np.array([E_(Ln,a) for a in la]); axm.plot(xy[:,0],xy[:,1],color="0.7",lw=0.3,alpha=.5,zorder=1)
for P in lat_t: xy=np.array([E_(o,P) for o in lo]); axm.plot(xy[:,0],xy[:,1],color="0.7",lw=0.3,alpha=.5,zorder=1)
xpos=[E_(Ln,win[2])[0] for Ln in lon_t]; ypos=[E_(win[0],P)[1] for P in lat_t]
axm.set_xticks(xpos); axm.set_xticklabels([f"{v:.0f}\u00b0E" for v in lon_t],fontsize=FS_SM)
axm.set_yticks(ypos); axm.set_yticklabels([f"{v:.1f}\u00b0N" for v in lat_t],fontsize=FS_SM)
axm.xaxis.set_minor_locator(AutoMinorLocator(5)); axm.yaxis.set_minor_locator(AutoMinorLocator(4))
axm.tick_params(which="both",length=3,width=0.7); axm.tick_params(which="minor",length=1.8)
axt=axm.secondary_xaxis("top",functions=(lambda x:x,lambda x:x))
axt.set_xticks(xpos); axt.set_xticklabels([f"{v:.0f}\u00b0E" for v in lon_t],fontsize=FS_SM)
axt.xaxis.set_minor_locator(AutoMinorLocator(5)); axt.tick_params(which="both",length=3,width=0.7); axt.tick_params(which="minor",length=1.8); axt.spines["top"].set_linewidth(0.8)
axr=axm.secondary_yaxis("right",functions=(lambda y:y,lambda y:y))
axr.set_yticks(ypos); axr.set_yticklabels([]); axr.yaxis.set_minor_locator(AutoMinorLocator(4)); axr.tick_params(which="both",length=3,width=0.7); axr.tick_params(which="minor",length=1.8); axr.spines["right"].set_linewidth(0.8)
axm.text(0.0,1.075,"(a)",transform=axm.transAxes,fontsize=FS_PANEL,fontweight="bold",va="bottom",ha="left")

# sub-region labels + CTI values (white text over coloured sea, halo for safety)
labpos={"Istanbul coast":(28.75,40.93),"Gulf of Izmit":(29.62,40.79),
        "Southern Marmara":(28.55,40.44),"Tekirdag/Canakkale":(27.05,40.66)}
for i,nm in enumerate(names):
    lo0,la0=labpos[nm]
    axm.text(*E_(lo0,la0),f"{nm}\nCTI {CTI[i]:.2f}",fontsize=FS_SM,fontweight="bold",ha="center",va="center",
             color="white",zorder=6,path_effects=[pe.withStroke(linewidth=1.8,foreground="#333")])

sbx=extent[0]+0.05*(extent[1]-extent[0]); sby=extent[2]+0.10*(extent[3]-extent[2])
axm.plot([sbx,sbx+50000],[sby,sby],color="#111",lw=2.2,solid_capstyle="butt",zorder=6)
for xx,l2 in [(sbx,"0"),(sbx+50000,"50 km")]:
    axm.plot([xx,xx],[sby,sby+2600],color="#111",lw=1.0,zorder=6); axm.text(xx,sby-3800,l2,ha="center",va="top",fontsize=FS_MID,zorder=6)
nx,ny=extent[0]+0.045*(extent[1]-extent[0]),extent[2]+0.955*(extent[3]-extent[2])
axm.annotate("N",xy=(nx,ny),xytext=(nx,ny-22000),ha="center",va="center",fontsize=FS_MID,fontweight="bold",
             zorder=6,arrowprops=dict(arrowstyle="-|>",color="#111",lw=1.6))
for sp in axm.spines.values(): sp.set_linewidth(0.8)

# CTI colour scale below the map
mcx=(L+0.5*map_w)/fig_w; bcx=(L+map_w+GAP+0.5*graph_w)/fig_w; mapbot=BAND/fig_h
cbw=0.55*map_w/fig_w
cax=fig.add_axes([mcx-cbw/2, mapbot-0.135, cbw, 0.026])
cb=fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation="horizontal")
cb.set_label("Composite transformation index, CTI (dimensionless, 0\u20131)",fontsize=FS_MID,labelpad=2)
cb.set_ticks([0,0.25,0.5,0.75,1.0]); cb.ax.tick_params(labelsize=FS_SM,length=2)

# (b) stacked component bars, ranked by CTI
order=np.argsort(-CTI)
COL={"R":"#c81d33","U":"#3b6ea5","M":"#e2a24a"}
yv=np.arange(len(names))[::-1]
axb.barh(yv,wr[order],color=COL["R"],height=0.6,label="Reclamation ($w$=0.40)",zorder=3)
axb.barh(yv,wu[order],left=wr[order],color=COL["U"],height=0.6,label="Urban/natural loss ($w$=0.35)",zorder=3)
axb.barh(yv,wm[order],left=wr[order]+wu[order],color=COL["M"],height=0.6,label="Mucilage ($w$=0.25)",zorder=3)
for j,yy in enumerate(yv):
    axb.text(CTI[order][j]+0.015,yy,f"{CTI[order][j]:.2f}",va="center",ha="left",fontsize=FS_MID,fontweight="bold",color="#333")
axb.set_yticks(yv); axb.set_yticklabels([names[k] for k in order],fontsize=FS_SM)
axb.set_ylim(-0.6,len(names)-0.4)
axb.set_xlabel("Composite index, CTI\n(weighted component contributions)",fontsize=FS_MID,labelpad=6)
axb.set_xlim(0,1.02)
axb.xaxis.set_minor_locator(AutoMinorLocator(2)); axb.tick_params(which="both",length=3,width=0.7,labelsize=FS_SM); axb.tick_params(which="minor",length=1.8)
axb.grid(axis="x",which="major",lw=0.4,color="0.85",zorder=0); axb.set_axisbelow(True)
for sp in ["top","right"]: axb.spines[sp].set_visible(False)
axb.text(0.0,1.075,"(b)",transform=axb.transAxes,fontsize=FS_PANEL,fontweight="bold",va="bottom",ha="left")

h,l=axb.get_legend_handles_labels()
fig.legend(h,l,loc="upper center",bbox_to_anchor=(bcx,mapbot-0.05),ncol=1,fontsize=FS_MID,frameon=False,handlelength=1.2,labelspacing=0.45,handletextpad=0.5,alignment="left")
fig.text(mcx,0.03,
   "Screening CTI (167 m grid): components min-max normalised across the four sub-regions.\n"
   "R' = reclamation intensity; U' = built-up density; M' = mean coastal MSI (fig08).",
   ha="center",va="bottom",fontsize=FS_SM,bbox=dict(boxstyle="round,pad=0.4",fc="#f6f4ef",ec="0.7",lw=0.6))

fig.savefig("/tmp/fig09_cti.pdf",bbox_inches="tight",pad_inches=0.02)
fig.savefig("/tmp/fig09_cti.png",dpi=600,bbox_inches="tight",pad_inches=0.02)
print("ranking:", " > ".join(names[k] for k in order))
print("saved fig %.2f x %.2f in"%(fig_w,fig_h))
