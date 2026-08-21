#!/usr/bin/env bash
# Run this INSIDE the folder that contains your unzipped "scripts/" directory.
# It removes macOS junk, the accidental duplicate, and moves the stray sample PNG.
set -e
find . -name '.DS_Store' -delete
find . -name '._*' -delete
rm -rf __MACOSX
rm -f "scripts/grass_lc_istanbul copy.sh"
mkdir -p examples
[ -f "scripts/fig_ndvi_2panel_2015_2025.png" ] && mv "scripts/fig_ndvi_2panel_2015_2025.png" examples/
echo "Cleanup done. Remaining scripts:"
ls -1 scripts/ | grep -v '^\.'
