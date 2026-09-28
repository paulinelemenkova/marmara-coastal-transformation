# Sea of Marmara coastal transformation — analysis scripts

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22048481.svg)](https://doi.org/10.5281/zenodo.22048481)

Open, reproducible workflow for mapping **land-cover change, shoreline reclamation and
mucilage susceptibility in the Sea of Marmara, Türkiye (2015–2025)** from open Landsat
imagery using **GRASS GIS**, **Python** and **GMT**.

These scripts reproduce every figure, table and derived layer in the article:

> Lemenkova, P., (2026). Land-cover change, shoreline reclamation and
> mucilage susceptibility in the Sea of Marmara, Türkiye. Manuscript in preparation.
> Code archived at https://doi.org/10.5281/zenodo.22048481

All input data are open access; no proprietary data or software are required.

## Requirements

- **Python 3.12** with the packages in [`requirements.txt`](requirements.txt)
  (`pip install -r requirements.txt`).
- **GRASS GIS ≥ 8.4** with the `r.learn.ml2` add-on
  (`g.extension extension=r.learn.ml2`) for the cluster-seeded Random Forest.
- **GMT ≥ 6** (Generic Mapping Tools) for the overview and shoreline maps.

## Data

Landsat Collection 2 Level-2 surface reflectance (WRS-2 paths 180 & 181, row 032; July
2015 and July 2025) is searched and downloaded by `scripts/search_download_landsat.py`
(via the Microsoft Planetary Computer STAC API). Ancillary open layers used in the paper:
Copernicus GLO-30 DEM, EMODnet/GEBCO bathymetry, ESA WorldCover and CORINE. Place the
downloaded scenes in a local `data/` folder (not tracked in git) before running the
GRASS workflow.

## How to reproduce

Run in this order (paths and mapset names are set at the top of each script):

1. **Download imagery** — `scripts/search_download_landsat.py`
2. **GRASS preprocessing + indices** — `scripts/grass_lc_istanbul.sh`
   (import, west–east mosaic, cloud mask, NDVI/MNDWI/NDBI, feature stack)
3. **Unsupervised seed** — `scripts/grass_unsup_classify.sh` (ISOCLUST → 10 classes)
4. **Cluster-seeded Random Forest** — `scripts/rf_seeded_from_isoclust.sh`
   (also `scripts/classify_scene.py`, `scripts/mosaic_classify_map_2015.py`,
   `scripts/mosaic_classify_map_2025.py`)
5. **Figures** — the `fig*` / `make_fig*` scripts (see index below)
6. **Results tables & statistics** — `scripts/fill_results.py`

## Script index

**Data & GRASS workflow**
- `search_download_landsat.py` — search/download Landsat C2 L2 scenes
- `grass_lc_istanbul.sh` — main GRASS pipeline: mosaic, cloud mask, spectral indices, stack
- `grass_unsup_classify.sh` — unsupervised ISOCLUST classification (training seed)
- `rf_seeded_from_isoclust.sh` — cluster-seeded Random Forest via `r.learn.ml2`
- `classify_scene.py`, `mosaic_classify_map_2015.py`, `mosaic_classify_map_2025.py` — per-epoch classification and mapping

**True-colour mosaics** (Figures 2–3)
- `mosaic_map.py` (2025), `mosaic_map_2015.py` (2015)

**Spectral indices** (Figures 7–10)
- NDVI: `ndvi_mosaic.py`, `ndvi_2panel_2015_2025_h.py`, `ndvi_2panel_2015_2025_v.py`
- MNDWI: `mndwi_mosaic.py`, `mndwi_2panel_2015_2025.py`, `mndwi_2panel_2015_2025_v.py`
- NDWI: `ndwi_mosaic.py`, `ndwi_mosaic_2015.py`, `ndwi_mosaic_2025.py`, `ndwi_2panel_2015_2025.py`
- NDBI: `ndbi_mosaic_seamless.py`, `ndbi_mosaic_2015.py`, `ndbi_2panel_2015_2025.py`
- `landcover_2panel.py` — unsupervised land-cover two-panel (Figure 6)

**Maps, figures & analysis**
- `fig01_overview.sh` — study-area / overview map, GMT (Figure 1)
- `fig02_pipeline.py` — processing pipeline diagram (Figure 4)
- `fig03_grassworkflow.py` — GRASS module workflow (Figure 5)
- `fig06_shoreline.gmt.sh` — shoreline change (EPR) map, GMT (Figure 12)
- `make_fig07_reclamation.py` — reclaimed-land map & statistics (Figure 13)
- `fig08_mucilage.py` — mucilage-susceptibility map (Figure 14)
- `fig09_cti.py` — composite coastal-transformation index (Figure 15)
- `fig10_comparison.py` — sub-region comparison (Figure 16)
- `make_fig11_drivers.py` — driver matrix (Figure 17)
- `fill_results.py` — land-cover composition and sub-region statistics (Tables 6–7)
- `figstyle.py` — shared Matplotlib style (imported by the plotting scripts)
- `epr_transects.txt` — shore-normal transect end-point rates (input to the shoreline map)

## Citing this repository

If you use these scripts, please cite both the article and this archived code (see
[`CITATION.cff`](CITATION.cff) and the Zenodo DOI above).

## License

Code is released under the [MIT License](LICENSE). Derived data files (e.g.
`epr_transects.txt`) are released under CC BY 4.0. Input datasets remain under their
original providers' licenses (USGS/NASA, Copernicus, ESA).

## Acknowledgement

Supported by the Scientific and Technological Research Council of Türkiye (TÜBİTAK),
BİDEB 2221 programme (reference B.14.2.TBT.0.06.01.02-220-859260).
