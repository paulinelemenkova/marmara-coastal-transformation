#!/usr/bin/env python3
"""fig10_comparison -- comparative sub-region CTI contributions + dominant land-cover
transitions (2015->2025).  (a) normalised reclamation (R'), urban/natural loss (U') and
mucilage-susceptibility (M') contributions per coastal sub-region; (b) Sankey of the
three super-class land-cover transitions over the coastal zone. Screening product
(see fig07-fig09 caveats)."""
import glob, numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.ticker import AutoMinorLocator
from scipy import ndimage

for f in glob.glob("/usr/share/fonts/truetype/liberation/LiberationSans-*.ttf"):
    mpl.font_manager.fontManager.addfont(f)
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":0.8})
FS_PANEL, FS_MID, FS_SM = 11, 9, 8

# ---- data from fig09 computation ----
ns={}; exec(open('/tmp/fig09_cti.py').read().split('# ===')[0].split('fig=plt.figure')[0], ns)
names=ns['names']; Rp,Up,Mp,CTI=ns['Rp'],ns['Up'],ns['Mp'],ns['CTI']
L15,L25,sea,coastland,pkm2=ns['L15'],ns['L25'],ns['sea'],ns['coastland'],ns['pkm2']
order=np.argsort(-CTI)
names_o=[names[k] for k in order]; R_o,U_o,M_o=Rp[order],Up[order],Mp[order]

# transition matrix (super-classes) over coastal zone + near-shore water
dist_land=ndimage.distance_transform_edt(~(L25>0))*166.8/1000.0
zone=coastland|(sea&(dist_land<5))
def sc(L):
    s=np.full(L.shape,-1,np.int8); s[L==0]=0; s[L==6]=1; s[(L>0)&(L!=6)]=2; return s
S15,S25=sc(L15),sc(L25); m=zone&(S15>=0)&(S25>=0)
T=np.array([[((S15==i)&(S25==j)&m).sum()*pkm2 for j in range(3)] for i in range(3)])
SLAB=["Water","Built-up","Other land"]

# =============================================================== figure =====
fig=plt.figure(figsize=(7.4,3.85))
gs=fig.add_gridspec(1,2,width_ratios=[1.15,1.0],wspace=0.28,left=0.09,right=0.985,top=0.91,bottom=0.34)
axa=fig.add_subplot(gs[0,0]); axb=fig.add_subplot(gs[0,1])

# ---- (a) grouped bars: R', U', M' per sub-region ----
COMP=[("Reclamation ($R'$)","#c81d33",R_o),("Urban/natural loss ($U'$)","#3b6ea5",U_o),
      ("Mucilage ($M'$)","#e2a24a",M_o)]
x=np.arange(len(names_o)); bw=0.26
for j,(lab,col,vals) in enumerate(COMP):
    axa.bar(x+(j-1)*bw,vals,bw,color=col,label=lab,zorder=3,edgecolor="white",linewidth=0.4)
axa.set_xticks(x); axa.set_xticklabels([n.replace(" ","\n",1) for n in names_o],fontsize=FS_SM)
axa.set_ylabel("Normalised component (0\u20131)",fontsize=FS_MID)
axa.set_ylim(0,1.08); axa.yaxis.set_minor_locator(AutoMinorLocator(2))
axa.tick_params(which="both",length=3,width=0.7,labelsize=FS_SM); axa.tick_params(which="minor",length=1.8)
axa.grid(axis="y",which="major",lw=0.4,color="0.85",zorder=0); axa.set_axisbelow(True)
for sp in ["top","right"]: axa.spines[sp].set_visible(False)
axa.text(-0.02,1.06,"(a)",transform=axa.transAxes,fontsize=FS_PANEL,fontweight="bold",va="bottom",ha="left")
axa.legend(fontsize=FS_SM,frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.14),ncol=3,
           handlelength=1.1,columnspacing=1.3,handletextpad=0.4)

# ---- (b) Sankey of super-class transitions ----
axb.set_xlim(0,1); axb.set_ylim(0,1); axb.axis("off")
axb.text(-0.02,1.06,"(b)",transform=axb.transAxes,fontsize=FS_PANEL,fontweight="bold",va="bottom",ha="left")
NCOL={"Water":"#5aa0d0","Built-up":"#c76b8e","Other land":"#9aa86a"}
tot=T.sum(); gap=0.04; xL,xR,nodew=0.06,0.94,0.055
src_tot=T.sum(1); tgt_tot=T.sum(0)
def stack(tots):
    h=(1-2*gap-(len(tots)-1)*gap)/1.0
    hs=(1-2*gap)*tots/tots.sum()-0  # proportional heights minus gaps handled below
    hs=(1-2*gap-(len(tots)-1)*gap)*tots/tots.sum()
    y=1-gap; out=[]
    for hh in hs: out.append((y-hh,y)); y-=hh+gap
    return out
Ly=stack(src_tot); Ry=stack(tgt_tot)
# draw nodes
for i,(y0,y1) in enumerate(Ly):
    axb.add_patch(Rectangle((xL,y0),nodew,y1-y0,fc=NCOL[SLAB[i]],ec="none",zorder=4))
    axb.text(xL-0.015,(y0+y1)/2,f"{SLAB[i]}",ha="right",va="center",fontsize=FS_SM,zorder=5)
for j,(y0,y1) in enumerate(Ry):
    axb.add_patch(Rectangle((xR-nodew,y0),nodew,y1-y0,fc=NCOL[SLAB[j]],ec="none",zorder=4))
    axb.text(xR+0.015,(y0+y1)/2,f"{SLAB[j]}",ha="left",va="center",fontsize=FS_SM,zorder=5)
# ribbons
Loff=[y1 for (y0,y1) in Ly]; Roff=[y1 for (y0,y1) in Ry]
for i in range(3):
    for j in range(3):
        f=T[i,j]
        if f<=0: continue
        h=(1-2*gap-2*gap)*f/tot*0  # placeholder
        hL=(Ly[i][1]-Ly[i][0])*f/src_tot[i]; hR=(Ry[j][1]-Ry[j][0])*f/tgt_tot[j]
        yL1=Loff[i]; yL0=yL1-hL; Loff[i]=yL0
        yR1=Roff[j]; yR0=yR1-hR; Roff[j]=yR0
        x0=xL+nodew; x1=xR-nodew; xm=(x0+x1)/2
        verts=[(x0,yL1),(xm,yL1),(xm,yR1),(x1,yR1),(x1,yR0),(xm,yR0),(xm,yL0),(x0,yL0),(x0,yL1)]
        codes=[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY]
        stable=(i==j)
        axb.add_patch(PathPatch(Path(verts,codes),fc=NCOL[SLAB[i]],ec="none",
                      alpha=0.75 if stable else 0.5,zorder=2))
axb.text(xL+nodew/2,1.0,"2015",ha="center",va="bottom",fontsize=FS_MID,fontweight="bold")
axb.text(xR-nodew/2,1.0,"2025",ha="center",va="bottom",fontsize=FS_MID,fontweight="bold")

fig.text(0.5,0.02,
   "Screening comparison (167 m grid). (b) coastal-zone super-class transitions; the reclamation signal is the "
   "Water$\\rightarrow$land flow,\nwhile Built-up$\\leftrightarrow$Other-land exchange lies within land-cover classification uncertainty.",
   ha="center",va="bottom",fontsize=FS_SM,bbox=dict(boxstyle="round,pad=0.4",fc="#f6f4ef",ec="0.7",lw=0.6))

fig.savefig("/tmp/fig10_comparison.pdf",bbox_inches="tight",pad_inches=0.02)
fig.savefig("/tmp/fig10_comparison.png",dpi=600,bbox_inches="tight",pad_inches=0.02)
print("Water->land (reclamation+prograd): %.0f km2 | Other->Built-up: %.0f | stable W/O/B: %.0f/%.0f/%.0f"%(
    T[0,1]+T[0,2],T[2,1],T[0,0],T[2,2],T[1,1]))
print("saved")
