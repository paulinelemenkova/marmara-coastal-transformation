#!/bin/bash
###############################################################################
# rf_seeded_from_isoclust.sh
# Cluster-seeded Random Forest land cover, Sea of Marmara 2015 & 2025.
#
# Idea: use the existing ISOCLUST 10-class maps as TRAINING LABELS for a
# supervised Random Forest run on the same feature stack. RF regularises the
# per-pixel clustering (uses the full stack context, gives class probabilities)
# and yields a defensible "supervised RF" map with the SAME 10 classes.
#
# Run inside YOUR GRASS mapset (the one that already holds the ISOCLUST result
# and the feature-stack bands), e.g.:
#   grass /path/to/grassdata/marmara/PERMANENT --exec bash rf_seeded_from_isoclust.sh
#
# GRASS >= 8.x. Needs the r.learn.ml2 addon (installed below). Network on.
###############################################################################
set -euo pipefail
export GRASS_OVERWRITE=1
export GRASS_MESSAGE_FORMAT=plain

# ----------------------------------------------------------------------------
# 0. EDIT THESE NAMES to match your mapset  (g.list type=raster  to check)
# ----------------------------------------------------------------------------
# ISOCLUST 10-class label rasters (integer 1..10), one per epoch:
ISO_2015="lc_iso_2015"
ISO_2025="lc_iso_2025"

# Feature-stack bands per epoch (the ~10 layers from your Preprocessing section:
# 6 reflective bands + NDVI + MNDWI + NBR + BSI). List them in the SAME order
# for both epochs so one model transfers cleanly.
STACK_2015="b2_2015,b3_2015,b4_2015,b5_2015,b6_2015,b7_2015,ndvi_2015,mndwi_2015,nbr_2015,bsi_2015"
STACK_2025="b2_2025,b3_2025,b4_2025,b5_2025,b6_2025,b7_2025,ndvi_2025,mndwi_2025,nbr_2025,bsi_2025"

CLASSES="1 2 3 4 5 6 7 8 9 10"   # ISOCLUST class codes
NPTS=600                          # training+val points PER class (stratified)
TRAIN_FRAC=0.70                   # 70% train / 30% held-out
NTREES=500                        # RF trees (report this in the paper)
SEED=42
OUT="/sandbox/output"; mkdir -p "$OUT"   # change to a normal dir if not in sandbox

# ----------------------------------------------------------------------------
# 1. Region + addon
# ----------------------------------------------------------------------------
g.region raster="$ISO_2025" -p          # both epochs share this grid
g.extension extension=r.learn.ml2 || true

# ----------------------------------------------------------------------------
# 2. Clean the ISOCLUST labels (mode filter) so RF learns clean class cores,
#    not salt-and-pepper edge pixels.  Comment out to train on raw clusters.
# ----------------------------------------------------------------------------
r.neighbors input="$ISO_2025" output=iso_clean method=mode size=3
# r.neighbors input="$ISO_2015" output=iso_clean_2015 method=mode size=3  # per-epoch variant

# ----------------------------------------------------------------------------
# 3. Stratified sample: NPTS random pixels per class -> sparse label maps,
#    then split into disjoint train / validation label rasters.
# ----------------------------------------------------------------------------
samp_list=""
for c in $CLASSES; do
  r.mapcalc "cls_${c} = if(iso_clean == ${c}, ${c}, null())"
  # r.random keeps the INPUT cell value at sampled cells (no attribute-column guesswork)
  r.random input="cls_${c}" npoints="$NPTS" raster="samp_${c}"
  samp_list="${samp_list}${samp_list:+,}samp_${c}"
done
r.patch input="$samp_list" output=all_samples

r.mapcalc "rnd = rand(0.0, 1.0)" seed="$SEED"
r.mapcalc "train_labels = if(rnd <= ${TRAIN_FRAC}, all_samples, null())"
r.mapcalc "val_labels   = if(rnd >  ${TRAIN_FRAC}, all_samples, null())"

# ----------------------------------------------------------------------------
# 4. Imagery groups (feature stacks)
# ----------------------------------------------------------------------------
i.group group=stack_2015 subgroup=stack_2015 input="$STACK_2015"
i.group group=stack_2025 subgroup=stack_2025 input="$STACK_2025"

# ----------------------------------------------------------------------------
# 5. Train ONE RF on 2025 seed labels, predict BOTH epochs.
#    Default = single model, so any 2015->2025 difference is driven by the
#    imagery, not by two different classifiers (cleanest for change detection).
#    NOTE: 2015 = Landsat 8, 2025 = Landsat 8/9 + different atmosphere, so a
#    2025-trained model is transferred across sensors — acceptable here because
#    your ISOCLUST already fixes prototypes on 2025, but state it in Methods.
# ----------------------------------------------------------------------------
r.learn.train group=stack_2025 training_map=train_labels \
  model_name=RandomForestClassifier n_estimators="$NTREES" \
  n_jobs=-1 save_model="$OUT/rf_marmara.gz"

r.learn.predict group=stack_2025 load_model="$OUT/rf_marmara.gz" output=rf_2025
r.learn.predict group=stack_2015 load_model="$OUT/rf_marmara.gz" output=rf_2015

# --- PER-EPOCH ALTERNATIVE (uncomment to train a separate model each year) ----
# r.neighbors input="$ISO_2015" output=iso_clean_2015 method=mode size=3
# ... rebuild all_samples_2015 / train_labels_2015 from iso_clean_2015 ...
# r.learn.train group=stack_2015 training_map=train_labels_2015 \
#   model_name=RandomForestClassifier n_estimators="$NTREES" save_model="$OUT/rf_2015.gz"
# r.learn.predict group=stack_2015 load_model="$OUT/rf_2015.gz" output=rf_2015
# ------------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# 6. Carry ISOCLUST colours + class labels onto the RF maps so the figure
#    styling is identical to your existing land-cover figure.
# ----------------------------------------------------------------------------
r.colors map=rf_2025 raster="$ISO_2025"
r.colors map=rf_2015 raster="$ISO_2025"
r.category map="$ISO_2025" > "$OUT/classes.txt"
r.category map=rf_2025 rules="$OUT/classes.txt"
r.category map=rf_2015 rules="$OUT/classes.txt"

# ----------------------------------------------------------------------------
# 7. FIDELITY check vs held-out ISOCLUST pixels.
#    *** This is NOT thematic accuracy *** — it only quantifies how well RF
#    reproduces the clustering. For the OA/kappa you REPORT, swap val_labels
#    for an independent reference raster (WorldCover/CORINE or interpreted
#    points), e.g.:  r.kappa -w classification=rf_2025 reference=ref_val ...
# ----------------------------------------------------------------------------
r.kappa -w classification=rf_2025 reference=val_labels \
  output="$OUT/rf2025_fidelity_to_isoclust.txt" || true

# ----------------------------------------------------------------------------
# 8. Export GeoTIFFs (drop straight into your existing 2-panel plot script)
# ----------------------------------------------------------------------------
r.out.gdal input=rf_2015 output="$OUT/rf_landcover_2015.tif" type=Byte \
  createopt="COMPRESS=DEFLATE" nodata=0
r.out.gdal input=rf_2025 output="$OUT/rf_landcover_2025.tif" type=Byte \
  createopt="COMPRESS=DEFLATE" nodata=0

# ----------------------------------------------------------------------------
# 9. (optional) quick-look PNG per epoch via headless cairo
# ----------------------------------------------------------------------------
for yr in 2015 2025; do
  export GRASS_RENDER_IMMEDIATE=cairo GRASS_RENDER_FILE_READ=TRUE
  export GRASS_RENDER_FILE="$OUT/rf_${yr}_quicklook.png"
  export GRASS_RENDER_FILE_WIDTH=1600 GRASS_RENDER_FILE_HEIGHT=1400
  d.erase
  d.rast "rf_${yr}"
  d.legend raster="rf_${yr}" -c title="Land cover ${yr}" at=5,55,2,5
  d.grid 0.5 text_color=white
  d.barscale ; d.northarrow at=92,15
  d.text text="RF (seeded from ISOCLUST) ${yr}" color=black size=3
done

echo "Done. Outputs in $OUT: rf_landcover_2015.tif, rf_landcover_2025.tif, quicklooks, fidelity report."
