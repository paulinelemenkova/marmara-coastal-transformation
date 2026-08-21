#!/usr/bin/env python3
"""Basin-wide 10-class land-cover map of the georeferenced Marmara mosaic
(UTM 35N, WRS-2 180/032 + 181/032, 2025). Classified by nearest 2025-east class
prototype (consistent convention); jet cpt; CORINE Level-3 legend; white lat/lon
graticule; metadata on the black free space; Arial (Liberation Sans)."""
import glob
import numpy as np
from PIL import Image
from sklearn.cluster import KMeans
from pyproj import Transformer
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
import textwrap

for f in glob.glob("/usr/share/fonts/truetype/liberation/LiberationSans-*.ttf"):
    mpl.font_manager.fontManager.addfont(f)
mpl.rcParams.update({"font.family":"sans-serif",
                     "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"]})

NAMES={1:"Water",2:"Forest",3:"Dense vegetation",4:"Shrub / transitional",5:"Irrigated crops",
       6:"Cropland",7:"Urban areas",8:"Agric. mosaic",9:"Dry grass / fallow",10:"Bare / sparse"}
EXPL={1:"Sea, water bodies, water courses (CLC 5.1-5.2)",2:"Broad-leaved / mixed forest (CLC 3.1)",
      3:"Transitional woodland-shrub, dense (CLC 3.2.4)",4:"Natural grasslands & shrub (CLC 3.2)",
      5:"Permanently irrigated / complex cultivation (CLC 2.1.2-2.4.2)",6:"Non-irrigated arable land (CLC 2.1.1)",
      7:"Urban fabric, industrial units (CLC 1.1-1.2)",8:"Complex cultivation / agriculture (CLC 2.4)",
      9:"Natural grasslands, sparsely vegetated (CLC 3.2-3.3)",10:"Sparsely vegetated / bare soil (CLC 3.3.3)"}
K=10; colors=plt.get_cmap("jet")(np.linspace(0,1,K))[:,:3]

# --- 2025-east reference prototypes ----------------------------------------- #
e=np.asarray(Image.open("/mnt/user-data/uploads/LC09_L1TP_180032_20250707_20250709_02_T1.jpg").convert("RGB"),np.float32)/255
fe=e.reshape(-1,3); ve=fe.sum(1)>24/255; V=fe[ve]
km0=KMeans(10,random_state=0,n_init=10).fit(V);dk=np.argsort(km0.cluster_centers_.sum(1))[:2];isw=np.isin(km0.labels_,dk)
kl=KMeans(9,random_state=0,n_init=10).fit(V[~isw]);lo=np.argsort(kl.cluster_centers_.sum(1));rm={o:n for n,o in enumerate(lo)}
fin=np.empty(V.shape[0],int);fin[isw]=0;fin[~isw]=np.vectorize(rm.get)(kl.labels_)+1
for a,b in [(8,9)]:
    va,vb=a-1,b-1;ma,mb=fin==va,fin==vb;fin[ma],fin[mb]=vb,va
proto=np.array([V[fin==c].mean(0) for c in range(10)])

# --- classify mosaic -------------------------------------------------------- #
d=np.load("/tmp/mos_utm.npz");out=d["out"].astype(np.float32)/255;valid=d["valid"]
Emin,Emax,Nmin,Nmax=d["ext"];OH,OW,_=out.shape
flat=out.reshape(-1,3)
cls=(((flat[:,None,:]-proto[None,:,:])**2).sum(2)).argmin(1)
L=np.full(flat.shape[0],-1,int);L[valid.ravel()]=cls[valid.ravel()];L=L.reshape(OH,OW)
tot=valid.sum();pct=[100*(L==c).sum()/tot for c in range(K)]
img=np.ones((OH,OW,3),np.float32)
for c in range(K): img[L==c]=colors[c]

Tf=Transformer.from_crs("EPSG:4326","EPSG:32635",always_xy=True)
def ll2utm(lo,la): return Tf.transform(lo,la)
W="black"; GRID="0.3"

fig=plt.figure(figsize=(17.5,9.4),facecolor="white")
ax=fig.add_axes([0.035,0.10,0.62,0.82]); ax.set_facecolor("white")
ax.imshow(img,extent=[Emin,Emax,Nmin,Nmax],origin="upper",interpolation="nearest")
ax.set_xlim(Emin,Emax);ax.set_ylim(Nmin,Nmax);ax.set_aspect("equal")

lon_t=[26,26.5,27,27.5,28,28.5,29,29.5,30];lat_t=[39.5,40,40.5,41]
la=np.linspace(39.25,41.40,80);loo=np.linspace(25.7,30.05,80)
for Ln in lon_t:
    xy=np.array([ll2utm(Ln,a) for a in la]);ax.plot(xy[:,0],xy[:,1],color=GRID,lw=0.5,alpha=0.8)
for P in lat_t:
    xy=np.array([ll2utm(o,P) for o in loo]);ax.plot(xy[:,0],xy[:,1],color=GRID,lw=0.5,alpha=0.8)
for Ln in lon_t:
    x,_=ll2utm(Ln,39.28)
    if Emin<x<Emax: ax.annotate(f"{Ln:.1f}\u00b0E",(x,Nmin),xytext=(0,-6),textcoords="offset points",ha="center",va="top",color=W,fontsize=8.5,clip_on=False)
for P in lat_t:
    _,y=ll2utm(25.72,P)
    if Nmin<y<Nmax: ax.annotate(f"{P:.1f}\u00b0N",(Emin,y),xytext=(-6,0),textcoords="offset points",ha="right",va="center",color=W,fontsize=8.5,clip_on=False)
for s in ax.spines.values(): s.set_color(W);s.set_linewidth(0.6)
ax.set_xticks([]);ax.set_yticks([])

# place labels
hw=[pe.withStroke(linewidth=2.2,foreground="black")]; hb=[pe.withStroke(linewidth=2.4,foreground="white")]
def txt(lon,lat,s,fs=11,col=W,it=True,he=hw,ha="center"):
    x,y=ll2utm(lon,lat); ax.text(x,y,s,color=col,fontsize=fs,style="italic" if it else "normal",
        fontweight="bold",ha=ha,va="center",path_effects=he)
txt(28.28,40.71,"Sea of Marmara",12,col="white")
txt(29.40,40.43,"Iznik Lake",9,col="white",ha="left")
txt(26.40,40.30,"Dardanelles",9,col="white")
x,y=ll2utm(28.98,41.01)
ax.plot(x,y,marker="s",color="red",markersize=10,markeredgecolor="white",markeredgewidth=1.3,
        path_effects=[pe.withStroke(linewidth=2.6,foreground="black")])
ax.text(x+2500,y,"Istanbul",color="black",fontsize=10,fontweight="bold",ha="left",va="center",path_effects=hb)

# title
fig.text(0.345,0.965,"Sea of Marmara \u2014 Land cover (WRS-2 180/032 + 181/032, 2025)",ha="center",color=W,fontsize=16,fontweight="bold")
fig.text(0.345,0.930,"Unsupervised 10-class classification \u2014 UTM Zone 35N / WGS 84",ha="center",color=W,fontsize=10)

# legend panel (right, on black)
lx=0.675; y=0.90
fig.text(lx,0.925,"Land cover classes (CORINE Level-3)",color=W,fontsize=11.5,fontweight="bold",va="top")
figH=fig.get_size_inches()[1]; lh=8.6*1.4/72/figH
for c in range(K):
    t=f"{c+1}. {NAMES[c+1]} \u2014 {EXPL[c+1]}  ({pct[c]:.1f}%)"
    wr=textwrap.fill(t,width=42); nl=wr.count("\n")+1
    fig.add_artist(Rectangle((lx,y-lh*0.95),0.016,lh*0.9,transform=fig.transFigure,
                   facecolor=tuple(np.clip(colors[c],0,1)),edgecolor=W,lw=0.5))
    fig.text(lx+0.022,y,wr,color=W,fontsize=8.6,va="top",linespacing=1.35)
    y-=nl*lh+0.010

# metadata (bottom-left black)
meta=("Scenes (Collection-2 L1, 2025):  West 181/032 L8 2025-07-22 (LC81810322025203LGN00);\n"
      "East 180/032 L9 2025-07-07 (LC91800322025188LGN00).  Projection UTM 35N / WGS 84, 250 m.\n"
      "Method: nearest-prototype to 2025 reference (GRASS i.cluster / i.maxlik equivalent). Source: authors.")
fig.text(0.037,0.055,meta,ha="left",va="top",color=W,fontsize=8.2,linespacing=1.5)

# scale bar
km=50.0; frac=km*1000/(Emax-Emin)*0.62; x0,y0=0.50,0.045
fig.add_artist(mpl.lines.Line2D([x0,x0+frac],[y0,y0],color=W,lw=2.2,transform=fig.transFigure))
for xx,lb in [(x0,"0"),(x0+frac/2,"25"),(x0+frac,"50 km")]:
    fig.add_artist(mpl.lines.Line2D([xx,xx],[y0,y0+0.008],color=W,lw=1.2,transform=fig.transFigure))
    fig.text(xx,y0-0.016,lb,ha="center",va="top",color=W,fontsize=8)
# north arrow
ax.annotate("N",xy=(0.965,0.86),xytext=(0.965,0.79),xycoords="axes fraction",ha="center",va="center",
            color=W,fontsize=13,fontweight="bold",arrowprops=dict(arrowstyle="-|>",color=W,lw=1.8))

fig.savefig("marmara_landcover_utm35n.png",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=0.05)
print("saved marmara_landcover_utm35n.png")
for c in range(K): print(f"  {c+1:2d} {NAMES[c+1]:20s} {pct[c]:5.1f}%")
