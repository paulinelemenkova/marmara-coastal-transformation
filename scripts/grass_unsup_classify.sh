#!/usr/bin/env bash
# =============================================================================
# grass_unsup_classify.sh
# Unsupervised 10-class land-cover classification of a Landsat scene in GRASS GIS
# 8.4, using i.cluster (spectral clustering) + i.maxlik (maximum-likelihood).
# Optionally relabels the 10 clusters against a CORINE reference by majority vote.
#
# Scene: Landsat 9 OLI, path 180 / row 032 (southern Sea of Marmara), 2025-07-07
#        e.g. LC09_L1TP_180032_20250707_20250709_02_T1  (Level-1) or the L2SP.
#
# Run (sandbox 'grass' runtime, or any GRASS 8.4 install):
#   bash grass_unsup_classify.sh
# =============================================================================
set -euo pipefail

# ----------------------------- CONFIG (edit) --------------------------------
export DATA="${DATA:-/sandbox/input}"          # folder with the untarred scene
export OUT="${OUT:-/sandbox/output}"
export GRASSDB="${GRASSDB:-/tmp/grassdb}"
export CRS="${CRS:-EPSG:32635}"                # UTM 35N (native for 180/032)
export SCENE="${SCENE:-LC09_L1TP_180032_20250707_20250709_02_T1}"
# Scheme: class 1 = water (MNDWI), classes 2..10 = 9 unsupervised LAND classes.
# Band suffixes to cluster on. Level-1: B2..B7. Level-2 SR: SR_B2..SR_B7.
# (must include B3 and B6 for the MNDWI water mask)
export BANDS="${BANDS:-B2 B3 B4 B5 B6 B7}"
# Optional CORINE raster (any categorical land-cover raster) for relabelling the
# clusters. Leave absent to keep generic Cluster 1..N labels.
export CORINE="${CORINE:-$DATA/corine_marmara.tif}"
# ----------------------------------------------------------------------------

mkdir -p "$OUT"; rm -rf "$GRASSDB"; mkdir -p "$GRASSDB"

# 'jet' colour table: class 1 = water (dark blue), 2..10 land -> red (r.colors rules=)
cat > "$OUT/cls_colors.txt" <<'COL'
1 0:0:128
2 0:0:255
3 0:96:255
4 0:212:255
5 77:255:170
6 170:255:77
7 255:230:0
8 255:122:0
9 255:19:0
10 128:0:0
nv 255:255:255
COL

cat > /tmp/_grass_body.sh <<'BODY'
#!/usr/bin/env bash
set -euo pipefail
export GRASS_OVERWRITE=1

echo "== 1. Import bands =="
inputs=""
for b in $BANDS; do
  f="$DATA/$SCENE/${SCENE}_${b}.TIF"
  r.import input="$f" output="b_${b}" --quiet
  inputs="${inputs:+$inputs,}b_${b}"
done
first=$(echo "$BANDS" | awk '{print $1}')
g.region raster="b_${first}" -p

echo "== 2. Mask scene border (drop black fill) =="
# valid where at least one band is non-zero
sumexpr=$(for b in $BANDS; do printf "b_%s+" "$b"; done | sed 's/+$//')
r.mapcalc "valid = if(($sumexpr) > 0, 1, null())" --quiet
r.mask raster=valid --quiet

echo "== 3. Band group =="
i.group group=scene subgroup=scene input="$inputs" --quiet

echo "== 4. Separate water with MNDWI (green=B3, swir1=B6) =="
# MNDWI > 0 = water; merges deep + shallow water into one class
r.mapcalc "mndwi = float(b_B3 - b_B6) / (b_B3 + b_B6)" --quiet
r.mapcalc "water = if(mndwi > 0.0, 1, null())" --quiet
r.mask -r --quiet
r.mapcalc "land = if(!isnull(valid) && isnull(water), 1, null())" --quiet

echo "== 5. Unsupervised classification of LAND into 9 (i.cluster + i.maxlik) =="
r.mask raster=land --quiet
i.cluster group=scene subgroup=scene signaturefile=sig classes=9 \
  reportfile="$OUT/cluster_report.txt" --quiet
i.maxlik  group=scene subgroup=scene signaturefile=sig \
  output=land_cls reject=land_reject --quiet
r.mask -r --quiet

echo "== 6. Combine: class 1 = water, 2..10 = land; jet colours + labels =="
r.mapcalc "lc_class = if(!isnull(water), 1, land_cls + 1)" --quiet
r.colors map=lc_class rules="$OUT/cls_colors.txt" --quiet
printf '1:Water\n' > /tmp/lbl.txt
for c in $(seq 2 10); do echo "$c:Class $c" >> /tmp/lbl.txt; done
r.category map=lc_class rules=/tmp/lbl.txt separator=":" --quiet

echo "== 6b. (optional) Relabel classes from CORINE by majority vote =="
if [ -f "$CORINE" ]; then
  r.import input="$CORINE" output=corine resample=nearest resolution=region --quiet
  r.stats -c -n input=lc_class,corine separator=space > "$OUT/class_vs_corine.txt"
  awk '{if($3>best[$1]){best[$1]=$3; maj[$1]=$2}}
       END{for(c in maj) print c" = "maj[c]}' \
      "$OUT/class_vs_corine.txt" | sort -n > /tmp/relabel.txt
  r.reclass input=lc_class output=lc_corine rules=/tmp/relabel.txt --quiet
  r.out.gdal input=lc_corine output="$OUT/lc10_corine.tif" type=Byte \
    createopt="COMPRESS=LZW" --quiet
  echo "   wrote CORINE-relabelled map: lc10_corine.tif"
fi

echo "== 7. Render map PNG (cairo) =="
export GRASS_RENDER_IMMEDIATE=cairo
export GRASS_RENDER_FILE="$OUT/LC09_180032_20250707_lc10.png"
export GRASS_RENDER_FILE_WIDTH=1500 GRASS_RENDER_FILE_HEIGHT=1200
export GRASS_RENDER_FILE_READ=TRUE
d.erase
d.rast map=lc_class
d.grid size=0.5 -g color=120:120:120 text_color=black fontsize=11
d.legend raster=lc_class -c -t at=6,58,86,89 fontsize=13 title="Class"
d.text text="Landsat 9 OLI - 2025-07-07 - Path 180 / Row 032 (southern Marmara)" \
  color=black size=2.6 at=2,97
d.text text="Data: USGS Landsat 9 C-2, ID: LC91800322025188LGN00 (2025/07/07). Method: k-means clustering (GRASS GIS i.cluster, i.maxlik). Source: authors." \
  color=50:50:50 size=1.5 at=2,3

echo "== 8. Export classified GeoTIFF + areas =="
r.out.gdal input=lc_class output="$OUT/lc10_clusters.tif" type=Byte \
  createopt="COMPRESS=LZW" --quiet
r.stats -a -n -l input=lc_class separator=tab > "$OUT/class_areas.txt"
echo "Done. Map + rasters + reports in $OUT"
BODY

echo "Creating GRASS project in $CRS and running ..."
grass -c "$CRS" "$GRASSDB/scene" --exec bash /tmp/_grass_body.sh
