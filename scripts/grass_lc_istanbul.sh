#!/usr/bin/env bash
# =============================================================================
# grass_lc_istanbul.sh
# Random-Forest land-cover classification of the Istanbul metropolitan
# Sea-of-Marmara coast from real Landsat Collection-2 Level-2 scenes, 2015 & 2025,
# in GRASS GIS 8.4 (produces the classified rasters behind fig04_lc_istanbul).
#
# SCENES (WRS-2 path 180 / row 31, Collection-2 Level-2 surface reflectance):
#   2015 : Landsat 8  OLI/TIRS  (summer, low cloud)
#   2025 : Landsat 9  OLI/TIRS  (summer, low cloud)   [Landsat 8 also fine]
# Download them first with search_download_landsat.py (companion script), then
# point DATA at the folder holding the two UNTARRED scene directories.
#
# Pipeline: import -> QA_PIXEL cloud/shadow mask -> SR scaling -> spectral
# indices -> band group -> Random Forest (r.learn.ml2) train+predict per epoch
# -> accuracy (r.kappa) -> change detection -> per-class areas -> GeoTIFF export.
#
# Run (sandbox 'grass' runtime, or any GRASS 8.4 install):
#   bash grass_lc_istanbul.sh
# =============================================================================
set -euo pipefail

# ----------------------------- CONFIG (edit) --------------------------------
export DATA="${DATA:-/sandbox/input}"       # holds the two untarred scene dirs
export OUT="${OUT:-/sandbox/output}"        # results go here
export GRASSDB="${GRASSDB:-/tmp/grassdb}"   # scratch GRASS database
export CRS="${CRS:-EPSG:32635}"             # UTM 35N — native CRS of path 180/31
export RES="${RES:-30}"                      # 30 m Landsat pixels

# Collection-2 Level-2 Product IDs of the two downloaded scenes.
# (edit to the exact IDs the search step returns for your chosen dates)
export S2015="${S2015:-LC08_L2SP_180031_20150724_20200909_02_T1}"
export S2025="${S2025:-LC09_L2SP_180031_20250728_20250730_02_T1}"

# Optional study-area clip box in CRS units (UTM 35N metres): W S E N.
# Leave empty ("") to keep the full scene. Values below crop to the Istanbul
# metropolitan Marmara shore (~Buyukcekmece -> Gebze).
export CLIP_W="${CLIP_W:-580000}"
export CLIP_S="${CLIP_S:-4530000}"
export CLIP_E="${CLIP_E:-700000}"
export CLIP_N="${CLIP_N:-4560000}"

# Training data: polygons/points per epoch with an INTEGER column 'class' (1..8).
# Build these from ESA WorldCover / CORINE / high-res imagery for each date.
# If you have only one set, point both at it (a temporal-transfer approximation).
export TRAIN2015="${TRAIN2015:-$DATA/train_2015.gpkg}"
export TRAIN2025="${TRAIN2025:-$DATA/train_2025.gpkg}"

# Optional independent validation points (column 'class') for accuracy.
export VALID="${VALID:-$DATA/validation.gpkg}"
# ----------------------------------------------------------------------------

mkdir -p "$OUT"
rm -rf "$GRASSDB"; mkdir -p "$GRASSDB"

# The 8-class coastal scheme (code -> label), matches the paper / fig04.
cat > "$OUT/classes.txt" <<'CLS'
1 Water
2 Built-up / impervious
3 Bare / reclaimed
4 Vegetation
5 Agriculture
6 Wetland
7 Port / industrial
8 Beach / sand
CLS

# GRASS colour rules for the 8 classes (r.colors rules= format).
cat > "$OUT/lc_colors.txt" <<'COL'
1 58:127:196
2 192:57:43
3 224:160:48
4 46:139:87
5 166:217:106
6 106:90:205
7 77:77:77
8 242:230:184
nv 255:255:255
COL

# ---- The GRASS job (runs inside the project; inherits the exported vars) ----
cat > /tmp/_grass_body.sh <<'BODY'
#!/usr/bin/env bash
set -euo pipefail
export GRASS_OVERWRITE=1

# helper: import the 6 optical SR bands + QA of one scene under a tag (e15/e25)
import_scene () {
  local pid="$1" tag="$2"
  local base="$DATA/$pid/$pid"
  local b
  for b in 2 3 4 5 6 7; do
    r.import input="${base}_SR_B${b}.TIF" output="${tag}_dn${b}" --quiet
  done
  r.import input="${base}_QA_PIXEL.TIF" output="${tag}_qa" --quiet
}

# helper: cloud mask (C2 QA_PIXEL bits) + SR scaling + spectral indices + group
prepare_epoch () {
  local tag="$1"
  # clear = not Fill(bit0) & not DilatedCloud(bit1) & not Cloud(bit3) & not Shadow(bit4)
  r.mapcalc "${tag}_clear = if( ((${tag}_qa & 1)==0) && ((${tag}_qa & 2)==0) && \
                                 ((${tag}_qa & 8)==0) && ((${tag}_qa & 16)==0), 1, null())" --quiet
  # C2 L2 surface-reflectance scaling: DN*0.0000275 - 0.2, masked to clear sky
  local b
  for b in 2 3 4 5 6 7; do
    r.mapcalc "${tag}_sr${b} = (${tag}_dn${b} * 0.0000275 - 0.2) * ${tag}_clear" --quiet
  done
  # indices: NDVI (veg), MNDWI (water/shoreline), NDBI (built-up)
  r.mapcalc "${tag}_ndvi  = float(${tag}_sr5 - ${tag}_sr4) / (${tag}_sr5 + ${tag}_sr4)" --quiet
  r.mapcalc "${tag}_mndwi = float(${tag}_sr3 - ${tag}_sr6) / (${tag}_sr3 + ${tag}_sr6)" --quiet
  r.mapcalc "${tag}_ndbi  = float(${tag}_sr6 - ${tag}_sr5) / (${tag}_sr6 + ${tag}_sr5)" --quiet
  # feature stack (identical band order for both epochs)
  i.group group="${tag}_stack" subgroup="${tag}_stack" \
    input="${tag}_sr2,${tag}_sr3,${tag}_sr4,${tag}_sr5,${tag}_sr6,${tag}_sr7,${tag}_ndvi,${tag}_mndwi,${tag}_ndbi" --quiet
}

# helper: rasterise a training vector to a label map aligned to the region
rasterise_training () {
  local src="$1" name="$2"
  v.import input="$src" output="${name}_v" --quiet
  v.to.rast input="${name}_v" output="$name" use=attr attribute_column=class --quiet
}

echo "== 1. Import scenes =="
import_scene "$S2015" e15
import_scene "$S2025" e25

echo "== 2. Region: full 2015 extent at ${RES} m, optional clip =="
g.region raster=e15_dn4 res="$RES" -a
if [ -n "${CLIP_W:-}" ]; then
  g.region n="$CLIP_N" s="$CLIP_S" e="$CLIP_E" w="$CLIP_W" res="$RES" -a
fi
g.region -p

echo "== 3. Cloud mask + SR scaling + indices + band group (both epochs) =="
prepare_epoch e15
prepare_epoch e25

echo "== 4. Random Forest (r.learn.ml2) =="
g.extension extension=r.learn.ml2 --quiet || true   # install addon (network on)
rasterise_training "$TRAIN2015" train15_lbl
r.learn.train  group=e15_stack training_map=train15_lbl \
  model_name=RandomForestClassifier n_estimators=500 \
  save_model="$OUT/rf_2015.gz" --quiet
r.learn.predict group=e15_stack load_model="$OUT/rf_2015.gz" output=lc_2015 --quiet

if [ -f "$TRAIN2025" ]; then TR25="$TRAIN2025"; else TR25="$TRAIN2015"; fi
rasterise_training "$TR25" train25_lbl
r.learn.train  group=e25_stack training_map=train25_lbl \
  model_name=RandomForestClassifier n_estimators=500 \
  save_model="$OUT/rf_2025.gz" --quiet
r.learn.predict group=e25_stack load_model="$OUT/rf_2025.gz" output=lc_2025 --quiet

# labels + colours
r.category map=lc_2015 rules="$OUT/classes.txt" separator=space
r.category map=lc_2025 rules="$OUT/classes.txt" separator=space
r.colors map=lc_2015 rules="$OUT/lc_colors.txt"
r.colors map=lc_2025 rules="$OUT/lc_colors.txt"

echo "== 5. Accuracy (independent validation, if provided) =="
if [ -f "$VALID" ]; then
  rasterise_training "$VALID" valid_lbl
  r.kappa -w classification=lc_2025 reference=valid_lbl output="$OUT/accuracy_2025.txt" || true
  cat "$OUT/accuracy_2025.txt" || true
fi

echo "== 6. Change detection + per-class areas =="
# from->to coded map (only where class changed)
r.mapcalc "lc_change = if(lc_2015 != lc_2025, lc_2015*10 + lc_2025, null())" --quiet
r.stats -a -n -l input=lc_2015 separator=tab > "$OUT/area_2015.txt"
r.stats -a -n -l input=lc_2025 separator=tab > "$OUT/area_2025.txt"
r.report -n map=lc_2015,lc_2025 units=k,p > "$OUT/lc_report.txt" || true
echo "-- 2015 areas (m2) --"; cat "$OUT/area_2015.txt"
echo "-- 2025 areas (m2) --"; cat "$OUT/area_2025.txt"

echo "== 7. Export GeoTIFFs =="
r.out.gdal input=lc_2015 output="$OUT/lc_istanbul_2015.tif" type=Byte \
  createopt="COMPRESS=LZW,TILED=YES" --quiet
r.out.gdal input=lc_2025 output="$OUT/lc_istanbul_2025.tif" type=Byte \
  createopt="COMPRESS=LZW,TILED=YES" --quiet
r.out.gdal input=lc_change output="$OUT/lc_istanbul_change.tif" type=Int16 \
  createopt="COMPRESS=LZW,TILED=YES" --quiet
echo "Done. Classified rasters + change map + area tables in $OUT"
BODY

echo "Creating GRASS project in $CRS and running the pipeline ..."
grass -c "$CRS" "$GRASSDB/marmara" --exec bash /tmp/_grass_body.sh
