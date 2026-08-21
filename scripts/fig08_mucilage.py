#!/usr/bin/env python3
"""fig08_mucilage -- screening mucilage-susceptibility (MSI) map of Marmara coastal waters.

MSI = 0.5*shallowness + 0.3*(poor flushing) + 0.2*(nutrient proximity), each rescaled to
[0,1] (Eq. MSI). Open-data proxies derived from the classified Landsat raster:
 shallowness  <- distance from shore (no bathymetry available; shore-distance proxy);
 flushing     <- geodesic distance within the sea from the open basin and the straits
                 (enclosed gulfs flush poorly); nutrient proximity <- distance to urban.
Screening product: coarse 167 m display grid; Bandirma-Erdek bay is not resolved as
water in the source classification and is therefore not assessed."""
import glob, numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pyproj import Transformer
from scipy import ndimage
from skimage.morphology import closing, disk, remove_small_objects
from skimage.graph import MCP_Geometric

for f in glob.glob("/usr/share/fonts/truetype/liberation/LiberationSans-*.ttf"):
    mpl.font_manager.fontManager.addfont(f)
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":0.8})

# ------------------------------------------------------------- georeference --
L25=np.load('/tmp/L25.npy')
Tf=Transformer.from_crs("EPSG:4326","EPSG:32635",always_xy=True)
Ti=Transformer.from_crs("EPSG:32635","EPSG:4326",always_xy=True)
Lf,Rf,Tp,Bp=98.5,2268.5,188.5,1604.5
Emin,Emax,Nmin,Nmax=392767.,754622.,4345794.,4581794.
mpp=(Emax-Emin)/(Rf-Lf); pkm2=(mpp*mpp)/1e6
def g2px(lo,la):
    E,N=Tf.transform(lo,la); return (int(round(Lf+(E-Emin)/(Emax-Emin)*(Rf-Lf))),int(round(Tp+(Nmax-N)/(Nmax-Nmin)*(Bp-Tp))))
def E_(lo,la): return Tf.transform(lo,la)
def px2E(px): return Emin+(px-Lf)/(Rf-Lf)*(Emax-Emin)
def py2N(py): return Nmax-(py-Tp)/(Bp-Tp)*(Nmax-Nmin)

# --------------------------------------------------------------- MSI model ---
water=(L25==0); land=(L25>0)
w=remove_small_objects(water,300); c=closing(w,disk(6)); f=ndimage.binary_fill_holes(c)
lb,n=ndimage.label(f); sz=ndimage.sum(np.ones_like(lb),lb,range(1,n+1))
sea=np.isin(lb,[i+1 for i in range(n) if sz[i]*pkm2>60])
dist_land=ndimage.distance_transform_edt(~land)*mpp/1000.0
seed=(sea&(dist_land>6)).copy()
for lo,la in [(29.05,41.08),(29.02,41.05),(26.62,40.32),(26.70,40.38),(26.50,40.22)]:
    x,y=g2px(lo,la); seed[max(y-3,0):y+4,max(x-3,0):x+4]|=sea[max(y-3,0):y+4,max(x-3,0):x+4]
cum,_=MCP_Geometric(np.where(sea,1.0,np.inf)).find_costs(list(zip(*np.where(seed))))
geo=np.where(np.isfinite(cum),cum,0)*mpp/1000.0
dist_urb=ndimage.distance_transform_edt(~(L25==6))*mpp/1000.0
def nz(a):
    v=a[sea]; lo,hi=np.percentile(v,2),np.percentile(v,98); return np.clip((a-lo)/(hi-lo),0,1)
shallow=nz(-dist_land); flush=nz(geo); nutrient=nz(-dist_urb)
Craw=0.5*shallow+0.3*flush+0.2*nutrient
mv=Craw[sea]; MSI=np.where(sea,(Craw-mv.min())/(mv.max()-mv.min()),np.nan)
THR=0.6
print("sea %.0f km2  high-susc(MSI>%.2f) %.0f km2 (%.0f%%)"%(sea.sum()*pkm2,THR,(sea&(MSI>THR)).sum()*pkm2,100*(sea&(MSI>THR)).sum()/sea.sum()))

# -------------------------------------------------------------- basemap RGB --
win=[26.45,30.02,40.20,41.22]
(cx0,cy0)=g2px(win[0],win[3]); (cx1,cy1)=g2px(win[1],win[2])
c0,c1,r0,r1=cx0,cx1,cy0,cy1
C_LAND="#e7e3da"; C_URB="#e6cde3"
cmap=LinearSegmentedColormap.from_list("msi",["#cfe6ef","#7 fcdbb".replace(" ",""),"#fee08b","#f28f2a","#c81d33"]) if False else plt.get_cmap("cividis")
sub=lambda A: A[r0:r1,c0:c1]
Ls,seaS,MS=sub(L25),sub(sea),sub(MSI)
H,Wd=Ls.shape
rgb=np.empty((H,Wd,3))
rgb[:]= mpl.colors.to_rgb("#dddcd6")                 # ALL land / non-sea = one monochrome tone
norm=mpl.colors.Normalize(0,1)
rgb[seaS]=cmap(norm(np.nan_to_num(MS[seaS])))[:,:3]
extent=[px2E(c0),px2E(c1),py2N(r1),py2N(r0)]

# =============================================================== figure =====
# Font policy (scientific-plotting SKILL): exactly 3 sizes, all within 8-12 pt.
FS_PANEL, FS_MID, FS_SM = 11, 9, 8      # panel tags / annotations+axis+legend / ticks+notes
import matplotlib.cm as cm
PAL=[mpl.colors.to_rgb(c) for c in ("#fb9a99","#fdbf6f","#b2df8a")]  # Paired light: red, orange, green

aspect=(extent[1]-extent[0])/(extent[3]-extent[2])
MAP_H=2.62; map_w=MAP_H*aspect; graph_w=2.05
L=0.62; GAP=0.92; R=0.14; TOPM=0.42; BAND=1.40
fig_w=L+map_w+GAP+graph_w+R; fig_h=TOPM+MAP_H+BAND
fig=plt.figure(figsize=(fig_w,fig_h))
axm=fig.add_axes([L/fig_w,BAND/fig_h,map_w/fig_w,MAP_H/fig_h])
axb=fig.add_axes([(L+map_w+GAP)/fig_w,BAND/fig_h,graph_w/fig_w,MAP_H/fig_h])

# ---- (a) map ----
axm.imshow(rgb,extent=extent,origin="upper",interpolation="nearest",aspect="equal")
axm.set_xlim(extent[0],extent[1]); axm.set_ylim(extent[2],extent[3])
lon_t=[27,28,29]; lat_t=[40.4,40.6,40.8,41.0,41.2]
la=np.linspace(win[2]-.05,win[3]+.05,80); lo=np.linspace(win[0]-.05,win[1]+.05,80)
for Ln in lon_t: xy=np.array([E_(Ln,a) for a in la]); axm.plot(xy[:,0],xy[:,1],color="0.7",lw=0.3,alpha=.5,zorder=1)
for P in lat_t: xy=np.array([E_(o,P) for o in lo]); axm.plot(xy[:,0],xy[:,1],color="0.7",lw=0.3,alpha=.5,zorder=1)
axm.set_xticks([E_(Ln,win[2])[0] for Ln in lon_t]); axm.set_xticklabels([f"{v:.0f}\u00b0E" for v in lon_t],fontsize=FS_SM)
axm.set_yticks([E_(win[0],P)[1] for P in lat_t]); axm.set_yticklabels([f"{v:.1f}\u00b0N" for v in lat_t],fontsize=FS_SM)
axm.xaxis.set_minor_locator(AutoMinorLocator(5)); axm.yaxis.set_minor_locator(AutoMinorLocator(4))
axm.tick_params(which="both",length=3,width=0.7); axm.tick_params(which="minor",length=1.8)
# WeSN frame: N = lon labels on top, e = ticks-only on right
xpos=[E_(Ln,win[2])[0] for Ln in lon_t]; ypos=[E_(win[0],P)[1] for P in lat_t]
axt=axm.secondary_xaxis("top",functions=(lambda x:x,lambda x:x))
axt.set_xticks(xpos); axt.set_xticklabels([f"{v:.0f}\u00b0E" for v in lon_t],fontsize=FS_SM)
axt.xaxis.set_minor_locator(AutoMinorLocator(5)); axt.tick_params(which="both",length=3,width=0.7); axt.tick_params(which="minor",length=1.8); axt.spines["top"].set_linewidth(0.8)
axr=axm.secondary_yaxis("right",functions=(lambda y:y,lambda y:y))
axr.set_yticks(ypos); axr.set_yticklabels([]); axr.yaxis.set_minor_locator(AutoMinorLocator(4)); axr.tick_params(which="both",length=3,width=0.7); axr.tick_params(which="minor",length=1.8); axr.spines["right"].set_linewidth(0.8)
axm.text(0.0,1.075,"(a)",transform=axm.transAxes,fontsize=FS_PANEL,fontweight="bold",va="bottom",ha="left")

def lab(txt,tlo,tla,alo,ala,c="#111"):
    axm.annotate(txt,xy=E_(alo,ala),xytext=E_(tlo,tla),fontsize=FS_MID,fontweight="bold",ha="center",va="center",
        color=c,zorder=7,path_effects=[pe.withStroke(linewidth=2.0,foreground="white")],
        arrowprops=dict(arrowstyle="-",color="#444",lw=0.6,shrinkA=1,shrinkB=2))
lab("Gulf of Izmit",29.70,41.02,29.68,40.76)
lab("Gulf of Gemlik",29.42,40.44,29.02,40.42)
lab("Istanbul coast",28.55,41.17,28.75,41.02)
lab("Bandirma-Erdek\n(not assessed)",27.55,40.24,27.85,40.36,c="#666")
lab("Bosphorus",29.22,41.13,29.00,41.00)
lab("Dardanelles",26.62,40.55,26.62,40.33)
# Central basin: white Helvetica-analogue, no halo, 9 pt (sits on dark-blue low-MSI water)
axm.text(*E_(28.02,40.72),"Central basin\n(low)",fontsize=FS_MID,family="sans-serif",color="white",
         ha="center",va="center",zorder=6)
# Sea of Marmara: white, no halo, at 28.00 E / 40.85 N
axm.text(*E_(28.00,40.85),"Sea of Marmara",fontsize=FS_MID,style="italic",family="sans-serif",color="white",
         ha="center",va="center",zorder=6)

sbx=extent[0]+0.05*(extent[1]-extent[0]); sby=extent[2]+0.10*(extent[3]-extent[2])
axm.plot([sbx,sbx+50000],[sby,sby],color="#111",lw=2.2,solid_capstyle="butt",zorder=6)
for xx,l2 in [(sbx,"0"),(sbx+50000,"50 km")]:
    axm.plot([xx,xx],[sby,sby+2600],color="#111",lw=1.0,zorder=6); axm.text(xx,sby-3800,l2,ha="center",va="top",fontsize=FS_MID,zorder=6)
nx,ny=extent[0]+0.045*(extent[1]-extent[0]),extent[2]+0.955*(extent[3]-extent[2])
axm.annotate("N",xy=(nx,ny),xytext=(nx,ny-22000),ha="center",va="center",fontsize=FS_MID,fontweight="bold",
             zorder=6,arrowprops=dict(arrowstyle="-|>",color="#111",lw=1.6))
for sp in axm.spines.values(): sp.set_linewidth(0.8)

# ---- (b) component contributions by sub-region ----
regions=[("Gulf of Izmit",29.35,30.0,40.66,40.83),
         ("Gulf of Gemlik",28.72,29.20,40.34,40.50),
         ("Istanbul coast",28.35,29.30,40.90,41.10),
         ("Central basin",27.10,28.40,40.55,40.85)]
sh_c=[];fl_c=[];nu_c=[]
for nm,a,b,cc,d in regions:
    (p0,q0)=g2px(a,d);(p1,q1)=g2px(b,cc)
    m=np.zeros_like(sea); m[q0:q1,p0:p1]=True; m&=sea
    sh_c.append(0.5*shallow[m].mean()); fl_c.append(0.3*flush[m].mean()); nu_c.append(0.2*nutrient[m].mean())
sh_c,fl_c,nu_c=map(np.array,(sh_c,fl_c,nu_c))
yv=np.arange(len(regions))[::-1]
axb.barh(yv,sh_c,color=PAL[0],height=0.55,label="Shallowness (0.5)",zorder=3)
axb.barh(yv,fl_c,left=sh_c,color=PAL[1],height=0.55,label="Poor flushing (0.3)",zorder=3)
axb.barh(yv,nu_c,left=sh_c+fl_c,color=PAL[2],height=0.55,label="Nutrient proximity (0.2)",zorder=3)
for yy,tot in zip(yv,sh_c+fl_c+nu_c): axb.text(tot+0.015,yy,f"{tot:.2f}",va="center",ha="left",fontsize=FS_MID,fontweight="bold",color="#333")
axb.set_yticks(yv); axb.set_yticklabels([r[0] for r in regions],fontsize=FS_SM)
axb.set_ylim(-0.6,len(regions)-0.35)
axb.set_xlabel("Mean weighted contribution\nto MSI (dimensionless)",fontsize=FS_MID,labelpad=6); axb.set_xlim(0,1.06)
axb.xaxis.set_minor_locator(AutoMinorLocator(2)); axb.tick_params(which="both",length=3,width=0.7,labelsize=FS_SM); axb.tick_params(which="minor",length=1.8)
axb.grid(axis="x",which="major",lw=0.4,color="0.85",zorder=0); axb.set_axisbelow(True)
for sp in ["top","right"]: axb.spines[sp].set_visible(False)
axb.text(0.0,1.02,"(b)",transform=axb.transAxes,fontsize=FS_PANEL,fontweight="bold",va="bottom",ha="left")

# MSI colour scale placed BELOW the map (a), horizontal
mcx=(L+0.5*map_w)/fig_w; bcx=(L+map_w+GAP+0.5*graph_w)/fig_w; mapbot=BAND/fig_h
cbw=0.55*map_w/fig_w
cax=fig.add_axes([mcx-cbw/2, mapbot-0.135, cbw, 0.026])
cb=fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation="horizontal")
cb.set_label("Mucilage-susceptibility index, MSI (dimensionless, 0\u20131)",fontsize=FS_MID,labelpad=2)
cb.set_ticks([0,0.25,0.5,0.75,1.0]); cb.ax.tick_params(labelsize=FS_SM,length=2)
# component legend BELOW panel (b), stacked
h,l=axb.get_legend_handles_labels()
fig.legend(h,l,loc="upper center",bbox_to_anchor=(bcx,mapbot-0.115),
           ncol=1,fontsize=FS_MID,frameon=False,handlelength=1.2,labelspacing=0.45,handletextpad=0.5,alignment="left")
# note below the colour scale (under the map)
fig.text(mcx,0.03,
   f"Screening MSI (167 m grid); high susceptibility (MSI > {THR:.1f}) = {(sea&(MSI>THR)).sum()*pkm2:.0f} km$^2$.\n"
   "Shallowness = shore-distance proxy (no bathymetry); Bandirma-Erdek not resolved as water.",
   ha="center",va="bottom",fontsize=FS_SM,bbox=dict(boxstyle="round,pad=0.4",fc="#f6f4ef",ec="0.7",lw=0.6))

fig.savefig("/tmp/fig08_mucilage.pdf",bbox_inches="tight",pad_inches=0.02)
fig.savefig("/tmp/fig08_mucilage.png",dpi=600,bbox_inches="tight",pad_inches=0.02)
print("region totals (MSI):",[round(float(s+f+nu),2) for s,f,nu in zip(sh_c,fl_c,nu_c)])
print("saved fig %.2f x %.2f in"%(fig_w,fig_h))
