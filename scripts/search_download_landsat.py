#!/usr/bin/env python3
"""
search_download_landsat.py
--------------------------
Find and download the REAL Landsat Collection-2 Level-2 scenes used by
grass_lc_istanbul.sh: WRS-2 path 180 / row 31 (Istanbul Marmara coast), one
low-cloud summer scene for 2015 and one for 2025.

Source: Microsoft Planetary Computer STAC (collection 'landsat-c2-l2').
No login/API key required; assets are signed automatically.

    pip install pystac-client planetary-computer requests
    python search_download_landsat.py

For each epoch it prints the scenes ranked by cloud cover, picks the clearest,
and downloads the six optical SR bands (B2-B7) + QA_PIXEL into
    ./landsat/<PRODUCT_ID>/<PRODUCT_ID>_SR_B*.TIF
which is exactly the layout grass_lc_istanbul.sh expects (set DATA=./landsat).

Alternative source (needs a free USGS EROS login):
    pip install landsatxplore
    landsatxplore search --dataset landsat_ot_c2_l2 \
        --location 41.01 28.98 --clouds 10 \
        --start 2015-06-01 --end 2015-09-30
"""

import os
import sys
import requests

try:
    from pystac_client import Client
    import planetary_computer as pc
except ImportError:
    sys.exit("pip install pystac-client planetary-computer requests")

# --------------------------- search parameters ------------------------------ #
STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLLECTION = "landsat-c2-l2"
PATH, ROW = 180, 31                       # WRS-2 path/row for Istanbul
MAX_CLOUD = 10                            # percent
EPOCHS = {
    "2015": "2015-06-01/2015-09-30",     # summer window, leaf-on, low sun-glint
    "2025": "2025-06-01/2025-09-30",
}
ASSETS = ["blue", "green", "red", "nir08", "swir16", "swir22", "qa_pixel"]
# Planetary Computer asset keys -> Landsat band-file suffixes (what GRASS expects)
SUFFIX = {"blue": "SR_B2", "green": "SR_B3", "red": "SR_B4", "nir08": "SR_B5",
          "swir16": "SR_B6", "swir22": "SR_B7", "qa_pixel": "QA_PIXEL"}
OUTDIR = "landsat"
# ---------------------------------------------------------------------------- #


def download(url, dest):
    if os.path.exists(dest):
        print(f"      exists  {os.path.basename(dest)}")
        return
    with requests.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    print(f"      saved   {os.path.basename(dest)}")


def main():
    client = Client.open(STAC, modifier=pc.sign_inplace)
    os.makedirs(OUTDIR, exist_ok=True)

    for epoch, dates in EPOCHS.items():
        print(f"\n=== {epoch}  (path {PATH} / row {ROW}, cloud < {MAX_CLOUD}%) ===")
        search = client.search(
            collections=[COLLECTION],
            datetime=dates,
            query={
                "landsat:wrs_path": {"eq": f"{PATH:03d}"},
                "landsat:wrs_row": {"eq": f"{ROW:03d}"},
                "eo:cloud_cover": {"lt": MAX_CLOUD},
            },
        )
        items = sorted(search.items(), key=lambda it: it.properties["eo:cloud_cover"])
        if not items:
            print("   no scene under the cloud threshold — widen the window/threshold")
            continue
        for it in items[:5]:
            print(f"   {it.id}   cloud={it.properties['eo:cloud_cover']:.1f}%  "
                  f"{it.properties['datetime'][:10]}")

        best = items[0]
        print(f"   -> downloading clearest: {best.id}")
        dst = os.path.join(OUTDIR, best.id)
        os.makedirs(dst, exist_ok=True)
        for key in ASSETS:
            if key not in best.assets:
                print(f"      MISSING asset {key}")
                continue
            fname = f"{best.id}_{SUFFIX[key]}.TIF"
            download(best.assets[key].href, os.path.join(dst, fname))

    print(f"\nDone. Set DATA={OUTDIR} and the S2015/S2025 product IDs "
          f"in grass_lc_istanbul.sh to the IDs printed above.")


if __name__ == "__main__":
    main()
