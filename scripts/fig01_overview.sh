#!/usr/bin/env bash
# fig01_overview.sh -- Sea of Marmara study-area / locator map (GMT 6 modern mode)
# Relief: GMT earth_relief_15s (SRTM15+ v2.7 / GEBCO family), same extent as the
# uploaded GEBCO 2026 tile (25-31E, 39-42N). Run: bash fig01_overview.sh
set -e
REG=25/31/39/42 ; PROJ=M16c
gmt set PROJ_LENGTH_UNIT cm
# exact plot size of the locator map -> inset frame sized to match (no empty white band)
DIMS=$(gmt mapproject -R23/45/35/43 -JM3.5c -W </dev/null)
IW=$(echo $DIMS | awk '{print $1}'); IH=$(echo $DIMS | awk '{print $2}')
gmt grdcut @earth_relief_15s -R$REG -Gmarmara.nc      # or: gmt grdconvert your_gebco.tif marmara.nc
gmt set FONT_ANNOT_PRIMARY 9p,Helvetica FONT_LABEL 10p,Helvetica FONT_TITLE 12p,Helvetica-Bold \
        MAP_TITLE_OFFSET 8p MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3.5p MAP_FRAME_PEN 1p \
        FORMAT_GEO_MAP dddF MAP_GRID_PEN_PRIMARY 0.25p,white PS_CHAR_ENCODING ISOLatin1+
gmt grdgradient marmara.nc -A300 -Ne0.6 -Gint.nc
gmt makecpt -Cabyss -T-1400/0/20 -Z > sea.cpt
gmt makecpt -C150/150/150,250/250/250 -T0/2600 > land.cpt   # greyscale land relief
gmt begin fig01_overview png,pdf
  # --- sea: shaded bathymetry ---
  gmt grdimage marmara.nc -R$REG -J$PROJ -Csea.cpt -Iint.nc -Baf \
      -BWSne+t"Sea of Marmara: study area, coastal sub-regions and Landsat WRS-2 coverage"
  # --- land: greyscale hillshade derived from the DEM, clipped to dry areas ---
  gmt coast -Gc -Dh -A15
    gmt grdimage marmara.nc -Cland.cpt -Iint.nc
  gmt coast -Q
  gmt basemap -Bg1                                          # white 1-deg grid lines
  gmt coast -Wthin,gray25 -N1/0.5p,gray45 -Dh -A15 -I1/0.4p,lightblue
  gmt plot -Sr+s -W1p,yellow,6_4:0 <<< $'27.80 40.10 30.00 41.75\n25.90 40.10 28.10 41.75\n27.60 39.00 29.80 40.15'
  gmt text -F+f7.5p,Helvetica-Bold,yellow+jLM <<< $'27.90 41.62 180/32\n26.00 41.62 181/32\n29.35 39.25 180/33'
  gmt text -F+f11p,Helvetica-BoldOblique,#89c3eb+jLM <<< '28.85 41.88 Black Sea'
  gmt text -F+f11p,Helvetica-BoldOblique,navy+jLM <<< '27.85 40.71 Sea of Marmara'
  gmt text -F+f11p,Helvetica-BoldOblique,navy+jCM <<< '25.55 40.58 Aegean Sea'
  gmt text -F+f8.5p,Helvetica-Oblique,gray15+jCM <<< $'30.02 40.93 Gulf of Izmit\n28.98 40.42 Gulf of Gemlik'
  gmt text -F+f8.5p,Helvetica-Oblique,gray15+jLM <<< '28.05 40.51 Kapidag Pen.'
  gmt plot -W0.7p,gray25 <<< $'>\n29.55 41.32\n29.04 41.14\n>\n25.55 39.96\n26.20 40.15\n>\n29.90 40.88\n29.65 40.73'
  gmt text -F+f8.5p,Helvetica-BoldOblique,gray15+jLM -N <<< '29.58 41.32 Bosphorus'
  gmt text -F+f8.5p,Helvetica-BoldOblique,gray15+jRM -N <<< '25.50 39.95 Dardanelles'
  gmt plot -Sc0.11c -Gblack -W0.3p,white <<< $'28.98 41.02\n27.51 40.98\n26.41 40.15\n27.97 40.36'
  gmt text -F+f8p,Helvetica,black+jBL -D0.10c/0.06c -Gwhite@40 <<< $'28.98 41.02 Istanbul\n27.51 40.98 Tekirdag'
  gmt text -F+f8p,Helvetica,black+jTL -D0.10c/-0.05c -Gwhite@40 <<< '27.97 40.36 Bandirma'
  gmt text -F+f8p,Helvetica,black+jTR -D-0.10c/-0.05c -Gwhite@40 <<< '26.41 40.15 Canakkale'
  gmt plot -Sr+s -W1.7p,red3 <<< $'28.45 40.78 29.40 41.25\n29.42 40.60 30.15 40.86\n27.45 40.17 29.20 40.55\n26.95 40.82 27.95 41.10\n26.05 39.98 26.88 40.40'
  gmt plot -Ss0.46c -Gred3 -W0.5p,white <<< $'28.53 41.18\n29.50 40.82\n27.55 40.49\n27.03 41.045\n26.13 40.34'
  gmt text -F+f9p,Helvetica-Bold,white+jCM <<< $'28.53 41.18 1\n29.50 40.82 2\n27.55 40.49 3\n27.03 41.045 4\n26.13 40.34 5'
  gmt basemap -TdjTL+w1.1c+f2+l,,,N+o0.35c/0.35c
  gmt basemap -LjBR+w100k+f+u+c40.5+o0.35c/0.45c
  gmt colorbar -Csea.cpt -DJMR+w7c/0.35c+o0.7c/0c -Bxa250f125+l"Depth" -By+l"m"
  printf '%s\n' 'H 9p,Helvetica-Bold Study sub-regions (red boxes)' \
    'L 8p,Helvetica L 1  Istanbul metropolitan coast' 'L 8p,Helvetica L 2  Gulf of Izmit' \
    'L 8p,Helvetica L 3  S. Marmara: Bandirma-Erdek-Kapidag-Gemlik' 'L 8p,Helvetica L 4  Tekirdag' \
    'L 8p,Helvetica L 5  \347anakkale Strait sector' 'L 7.5p,Helvetica-Oblique L Dashed yellow: Landsat WRS-2' > legend.txt
  gmt legend legend.txt -DjBL+w6.6c+o0.25c/0.25c -F+gwhite@10+p0.6p+c0.10c
  # locator inset flush in the top-right corner, frame == map size (no white box above)
  gmt inset begin -DjTR+w${IW}c/${IH}c+o0.15c/0.15c -F+gwhite+p0.8p
    gmt coast -R23/45/35/43 -JM3.5c -Ggray72 -Slightsteelblue1 -N1/0.4p,gray50 -Wfaint -Da -A300
    gmt plot -Sr+s -W1.3p,red3 <<< '25 39 31 42'
    gmt text -F+f6.5p,Helvetica-Oblique,gray20+jCM <<< $'36.0 38.6 Turkiye\n27.3 37.2 Aegean'
  gmt inset end
  gmt text -R$REG -J$PROJ -F+f8p,Helvetica-Oblique,gray30+jBL -N -Gwhite@30 \
    <<< '25.02 38.74 Data: SRTM15+ v2.7 / GEBCO family (GMT earth_relief 15 arcsec). Projection: Mercator.'
gmt end
