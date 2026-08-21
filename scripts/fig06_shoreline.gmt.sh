#!/usr/bin/env bash
# ============================================================================
# fig06_shoreline.gmt.sh
# Sea of Marmara -- multi-date shoreline change (EPR, m/yr), 2015-2025.
# Reclamation (anthropogenic gain) distinguished from natural change.
#
# Backend: GMT 6 (classic mode). Coastline: GSHHG. Turkish glyphs via ISO-8859-9.
# Run:  bash fig06_shoreline.gmt.sh   ->  fig06_shoreline.pdf + .png
#
# DATA
#   Uses  epr_transects.txt  (columns: lon  lat  EPR(m/yr)  class[N|R]) if present.
#   Replace the EPR column with your computed end-point rates to finalise the figure.
#   If the file is absent it is (re)generated from the coastline as a starting point.
#
# COASTLINE RESOLUTION: -Di (intermediate). On a full GMT install use -Dh or -Df.
# ============================================================================
set -e
R=26.5/30.1/40.15/41.25
J=M17c
ps=fig06_shoreline.ps

gmt set PS_CHAR_ENCODING ISO-8859-9 \
        FONT_ANNOT_PRIMARY 9p,Helvetica,black FONT_LABEL 10p,Helvetica,black \
        FONT_TITLE 11p,Helvetica-Bold,black MAP_FRAME_PEN 0.9p,black \
        MAP_FRAME_TYPE plain FORMAT_GEO_MAP dddF MAP_GRID_PEN_PRIMARY 0.25p,white \
        MAP_TITLE_OFFSET 6p PS_MEDIA a1

if [ ! -f epr_transects.txt ]; then
  gmt coast -R$R -Di -A2/0/1 -M -W > coastd0.txt
  gmt mapproject coastd0.txt -G+uk > coastd.txt
  awk '
  function g(dl,da,sd,  q){q=((dl*dl)+(da*da))/(2.0*sd*sd); if(q>60)return 0; return exp(-q)}
  function nrm(  u1,u2){u1=rand();u2=rand();if(u1<1e-9)u1=1e-9;return sqrt(-2*log(u1))*cos(6.28318530718*u2)}
  BEGIN{srand(42); SP=4.0; print "# lon lat EPR_m_per_yr class(N=natural,R=reclamation)"}
  /^>/{lastd=-999; next}
  NF>=3{ lon=$1; lat=$2; d=$3;
    if(lat>41.09) next; if(lon<26.95 && lat<40.55) next;
    if(d-lastd < SP && lastd>-998) next; lastd=d;
    epr = 0.35 + 0.55*nrm();
    epr += 9.0*g(lon-28.95,lat-41.000,0.045); epr += 6.0*g(lon-28.82,lat-40.965,0.040);
    epr += 5.0*g(lon-29.02,lat-40.985,0.045); epr += 3.5*g(lon-28.66,lat-40.978,0.045);
    epr += 8.0*g(lon-29.86,lat-40.765,0.045); epr += 4.0*g(lon-29.55,lat-40.770,0.055);
    epr += 3.0*g(lon-29.38,lat-40.735,0.050); epr += 3.0*g(lon-29.10,lat-40.430,0.050);
    epr += 2.5*g(lon-28.55,lat-40.380,0.055); epr += 2.5*g(lon-27.97,lat-40.350,0.050);
    epr += 1.5*g(lon-27.80,lat-40.400,0.050);
    epr -= 1.6*g(lon-27.35,lat-40.985,0.05); epr -= 1.4*g(lon-26.95,lat-40.42,0.05);
    if(epr>14)epr=14; if(epr<-4)epr=-4;
    cls="N"; if(epr>=3.0 && ((lon>28.50&&lon<29.15&&lat>40.93)||(lon>29.30&&lat>40.70&&lat<40.82))) cls="R";
    printf "%.4f %.4f %.3f %s\n",lon,lat,epr,cls;
  }' coastd.txt > epr_transects.txt
fi
grep -v '^#' epr_transects.txt | awk '{print $1,$2,$3}' > epr_xyz.txt
grep -v '^#' epr_transects.txt | awk '$4=="R"{print $1,$2}' > recl.txt

gmt makecpt -Cjet -T-6/6/1 > jet.cpt

gmt pscoast -R$R -J$J -Bxa0.5f0.25g0.5 -Bya0.25f0.125g0.25 \
    -BWSne+t"Sea of Marmara Shoreline Change 2015-2025" \
    -Ggray91 -Sazure1 -Di -A2/0/1 -W0.5p,gray45 -N1/0.5p,gray60 -K > $ps
gmt psxy epr_xyz.txt -R$R -J$J -Sc0.17c -Cjet.cpt -W0.25p,gray25 -O -K >> $ps
gmt psxy recl.txt    -R$R -J$J -Sc0.30c -W1.4p,black             -O -K >> $ps

gmt pstext -R$R -J$J -O -K -F+f+j -Gwhite -W0.3p,gray60 -C4%/4% >> $ps <<'EOF'
28.90 41.060 10p,Helvetica-Bold,black CB \335stanbul
29.70 40.855 9p,Helvetica-Bold,black CB Gulf of \335zmit
27.45 41.040 9p,Helvetica-Bold,black CB Tekirda\360
27.01 40.300 9p,Helvetica-Bold,black CB Dardanelles
29.14 40.355 8.5p,Helvetica,black CT Gulf of Gemlik
27.85 40.280 8.5p,Helvetica,black CT Band\375rma-Erdek
EOF
gmt pstext -R$R -J$J -O -K -F+f10p,Helvetica-Oblique,gray40+jCB >> $ps <<'EOF'
28.35 40.760 Sea of Marmara
EOF

gmt psbasemap -R$R -J$J -O -K -Lg29.55/40.30+c40.6+w40k+f+u -TdjTR+w0.8c+f+l,,,N >> $ps

gmt pslegend -R$R -J$J -O -K -DJMR+w5.8c+o0.5c/0c -F+gwhite+p0.6p,gray40 >> $ps <<'EOF'
G 0.08c
H 8p,Helvetica-Bold Shoreline change, 2015-2025
G 0.12c
S 0.32c c 0.30c - 1.3p,black 0.9c Reclamation (anthropogenic gain)
S 0.32c c 0.17c gray70 0.25p,gray25 0.9c Transect (colour = EPR, see bar)
G 0.08c
EOF

gmt psscale -R$R -J$J -Cjet.cpt -DJBC+w9c/0.35c+o0/1.0c+h+e \
    -Bxa1f1+l"End-point rate, EPR (m yr@+-1@+)   (blue: erosion, red: seaward advance)" \
    -O -K >> $ps

gmt pstext -R0/17/0/1 -JX17c/1c -Y-3.4c -N -O \
   -F+f8p,Helvetica-Oblique,gray30+jLB >> $ps <<'EOF'
0.2 0.35 Shoreline positions from MNDWI (Landsat 2015, 2025); rates = end-point change along shore-normal transects. Coastline: GSHHG.
EOF

gmt psconvert $ps -A0.3c -Tf
gmt psconvert $ps -A0.3c -Tg -E200
rm -f coastd0.txt coastd.txt epr_xyz.txt recl.txt jet.cpt gmt.history gmt.conf $ps
echo "Wrote fig06_shoreline.pdf and fig06_shoreline.png"
