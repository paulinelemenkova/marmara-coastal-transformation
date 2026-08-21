#!/usr/bin/env python3
r"""
fill_results.py -- compute the Results numbers for the Sea of Marmara paper from
standard GRASS outputs and fill the LaTeX \ph{} placeholders in the two
subsections "Classification Accuracy" and "Land Cover Change and Urbanization
(2015--2025)" plus Table~\ref{tab:accuracy}.

It does NOT invent anything: every number is derived from data you export from
GRASS. If you don't pass the data, it fills nothing.

--------------------------------------------------------------------------------
INPUTS (CSV, easiest to export straight from GRASS; comma-separated)
--------------------------------------------------------------------------------
1) confusion.csv  -- one confusion count per row, per sub-region:
       predicted,reference,subregion,count
   Export with (labels via -l keep it readable):
       r.stats -cn -l input=classified,reference,subregions separator=comma \
               output=confusion.csv
   (3-column predicted,reference,count is accepted too, but then the per-region
    accuracy table cannot be filled.)

2) areas_2015.csv and areas_2025.csv -- class area per sub-region, one per epoch:
       class,subregion,area
   Export each epoch with:
       r.stats -an -l input=class_2015,subregions separator=comma \
               output=areas_2015.csv
       r.stats -an -l input=class_2025,subregions separator=comma \
               output=areas_2025.csv
   r.stats -a reports area in square metres -> use --area-units m2 (default).
   (2-column class,area is accepted; then only basin-wide stats are computed.)

Header rows are auto-detected and skipped.

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python3 fill_results.py \
      --tex article_marmara_coasts.tex \
      --confusion confusion.csv \
      --areas2015 areas_2015.csv --areas2025 areas_2025.csv \
      --emit-blocks results_filled.tex \
      --write-tex article_marmara_coasts.filled.tex     # optional full copy

  python3 fill_results.py --selftest    # runs the maths on SYNTHETIC numbers
                                        # (prints to screen only; writes nothing)

Add --turkish to render sub-region names with Turkish diacritics
(\.{I}stanbul, Gulf of \.{I}zmit, Tekirda\u{g}/\c{C}anakkale).
--------------------------------------------------------------------------------
EDIT THE CONFIG BLOCK BELOW so the class ids/labels and sub-region ids in your
CSVs map to the right semantic groups. The script prints the grouping it used so
you can check it.
"""
import argparse
import csv
import math
import re
import sys
from collections import defaultdict

# ============================ CONFIG (edit me) ==============================
YEAR1, YEAR2 = 2015, 2025

# Sub-region ids EXACTLY as they appear in the CSV `subregion` column, in the
# order the accuracy table lists them. If you exported with -l these are labels;
# otherwise they are category numbers (put them as strings).
SUBREGION_ORDER = ["Istanbul coast", "Gulf of Izmit",
                   "Southern Marmara", "Tekirdag/Canakkale"]

# Which sub-region id plays which role in the prose sentences.
IST_SUB   = "Istanbul coast"        # "The Istanbul metropolitan coast ..."
SOUTH_SUB = "Southern Marmara"      # "In the southern Marmara ..."

# Display names printed into the LaTeX (plain to match the current file).
SUBREGION_DISPLAY = {
    "Istanbul coast":     "Istanbul coast",
    "Gulf of Izmit":      "Gulf of Izmit",
    "Southern Marmara":   "Southern Marmara",
    "Tekirdag/Canakkale": "Tekirdag/Canakkale",
}
# Turkish (ASCII-escaped) variant, used with --turkish.
SUBREGION_DISPLAY_TR = {
    "Istanbul coast":     r"\.{I}stanbul coast",
    "Gulf of Izmit":      r"Gulf of \.{I}zmit",
    "Southern Marmara":   "Southern Marmara",
    "Tekirdag/Canakkale": r"Tekirda\u{g}/\c{C}anakkale",
}

# Map each land-cover class id/label (as it appears in the CSV `class` and in the
# confusion `predicted`/`reference` columns) to a FINE group tag. Defaults follow
# the paper's ten CORINE Level-3 classes; replace the keys with YOUR ids/labels.
CLASS_TO_GROUP = {
    "water":                         "water",
    "forest":                        "forest",
    "dense vegetation":              "dense_veg",
    "shrub/transitional":            "shrub",
    "irrigated crops":               "irrig_crops",
    "cropland":                      "cropland",
    "urban and built-up surfaces":   "built_up",
    "agricultural mosaic":           "agri_mosaic",
    "dry grassland/fallow":          "grassland",
    "bare or sparsely vegetated":    "bare",
}
# Human labels for the "spectral overlap between the X and Y classes" sentence.
# Defaults to identity (i.e. the CSV already carries readable labels).
CLASS_LABELS = {}

# Super-groups built from the fine tags above (edit membership if needed).
G_BUILT_UP           = {"built_up"}
G_WATER              = {"water"}
G_NATURAL_SEMINAT    = {"forest", "dense_veg", "shrub", "grassland", "bare"}
G_VEGETATION_OPEN    = {"forest", "dense_veg", "shrub", "grassland", "bare"}
G_AGRI_WETLAND       = {"irrig_crops", "cropland", "agri_mosaic", "wetland"}
# ===========================================================================

BASIN = "__basin__"


# ------------------------------- loaders -----------------------------------
def _looks_like_header(row, numeric_col):
    try:
        float(row[numeric_col]); return False
    except (ValueError, IndexError):
        return True


def load_confusion(path):
    """Return conf[sub][pred][ref] = count. sub='__all__' if no subregion col."""
    conf = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r and any(c.strip() for c in r)]
    if not rows:
        raise ValueError(f"{path}: empty")
    ncol = len(rows[0])
    if _looks_like_header(rows[0], ncol - 1):
        rows = rows[1:]
    for r in rows:
        r = [c.strip() for c in r]
        if ncol >= 4:
            pred, ref, sub, cnt = r[0], r[1], r[2], r[3]
        elif ncol == 3:
            pred, ref, sub, cnt = r[0], r[1], "__all__", r[2]
        else:
            raise ValueError(f"{path}: need >=3 columns, got {ncol}")
        conf[sub][pred][ref] += int(round(float(cnt)))
    return conf


_UNIT = {"m2": 1e-6, "ha": 1e-2, "km2": 1.0}  # -> km^2


def load_areas(path, units):
    """Return A[sub][class] = area_km2. sub='__all__' if no subregion col."""
    A = defaultdict(lambda: defaultdict(float))
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r and any(c.strip() for c in r)]
    if not rows:
        raise ValueError(f"{path}: empty")
    ncol = len(rows[0])
    if _looks_like_header(rows[0], ncol - 1):
        rows = rows[1:]
    f = _UNIT[units]
    for r in rows:
        r = [c.strip() for c in r]
        if ncol >= 3:
            cls, sub, area = r[0], r[1], r[2]
        elif ncol == 2:
            cls, sub, area = r[0], "__all__", r[1]
        else:
            raise ValueError(f"{path}: need >=2 columns, got {ncol}")
        A[sub][cls] += float(area) * f
    return A


# ---------------------------- accuracy maths -------------------------------
def confusion_metrics(mat):
    """mat[pred][ref]=count -> dict with oa, kappa, weighted f1, per-class pa/ua,
    and the largest off-diagonal (pred,ref) pair."""
    classes = sorted(set(mat) | {r for p in mat for r in mat[p]})
    n = {(p, r): mat.get(p, {}).get(r, 0) for p in classes for r in classes}
    N = sum(n.values())
    if N == 0:
        raise ValueError("empty confusion matrix")
    diag = sum(n[(c, c)] for c in classes)
    row = {p: sum(n[(p, r)] for r in classes) for p in classes}   # predicted tot
    col = {r: sum(n[(p, r)] for p in classes) for r in classes}   # reference tot
    oa = 100.0 * diag / N
    pe = sum(row[c] * col[c] for c in classes) / (N * N)
    kappa = (oa / 100.0 - pe) / (1 - pe) if (1 - pe) else float("nan")
    pa, ua, f1 = {}, {}, {}
    for c in classes:
        pa[c] = 100.0 * n[(c, c)] / col[c] if col[c] else float("nan")   # producer
        ua[c] = 100.0 * n[(c, c)] / row[c] if row[c] else float("nan")   # user
        s = pa[c] + ua[c]
        f1[c] = (2 * pa[c] * ua[c] / s) if s else float("nan")
    support = col                                       # weight by true support
    wsum = sum(support[c] for c in classes if not math.isnan(f1[c]))
    wf1 = (sum(support[c] * f1[c] for c in classes if not math.isnan(f1[c])) / wsum
           / 100.0) if wsum else float("nan")           # 0..1 scale
    off = [(n[(p, r)], p, r) for p in classes for r in classes if p != r]
    off.sort(reverse=True)
    worst_pair = (off[0][2], off[0][1]) if off and off[0][0] > 0 else (None, None)
    return dict(oa=oa, kappa=kappa, wf1=wf1, pa=pa, ua=ua, f1=f1,
                worst_pair=worst_pair)


# ------------------------------ formatting ---------------------------------
def f_oa(x):    return f"{x:.1f}"
def f_kap(x):   return f"{x:.2f}"
def f_pct(x):   return f"{x:.0f}"
def f_km(x):    return f"{x:.0f}" if abs(x) >= 10 else f"{x:.1f}"
def f_rate(x):  return f"{x:.1f}"


def group_area(A, sub, fine_set):
    tot = 0.0
    for cls, a in A[sub].items():
        if CLASS_TO_GROUP.get(cls) in fine_set:
            tot += a
    return tot


def basin_area(A, fine_set, subs):
    return sum(group_area(A, s, fine_set) for s in subs)


# ------------------------------ assembly -----------------------------------
def build_values(conf, A1, A2, disp, area_subs):
    per_sub_ok = "__all__" not in conf
    if not per_sub_ok:
        raise SystemExit("ERROR: confusion.csv has no subregion column; the "
                         "per-region accuracy table cannot be filled. Re-export "
                         "with  r.stats -cn -l input=classified,reference,"
                         "subregions ...")
    metrics = {s: confusion_metrics(conf[s]) for s in conf}
    for s in SUBREGION_ORDER:
        if s not in metrics:
            raise SystemExit(f"ERROR: sub-region '{s}' not found in confusion.csv "
                             f"(present: {sorted(metrics)}). Fix SUBREGION_ORDER.")

    oas = {s: metrics[s]["oa"] for s in SUBREGION_ORDER}
    kaps = {s: metrics[s]["kappa"] for s in SUBREGION_ORDER}
    f1s = {s: metrics[s]["wf1"] for s in SUBREGION_ORDER}
    best = max(SUBREGION_ORDER, key=lambda s: oas[s])
    worst = min(SUBREGION_ORDER, key=lambda s: oas[s])
    a_ref, b_ref = metrics[worst]["worst_pair"]
    lab = lambda c: CLASS_LABELS.get(c, c) if c is not None else "n/a"

    # water & built-up producer/user accuracy floor across all sub-regions
    wb = []
    for s in SUBREGION_ORDER:
        m = metrics[s]
        for c in m["pa"]:
            g = CLASS_TO_GROUP.get(c)
            if g in (G_WATER | G_BUILT_UP):
                wb += [m["pa"][c], m["ua"][c]]
    wb = [x for x in wb if not math.isnan(x)]
    pu_threshold = int(math.floor(min(wb))) if wb else 0

    # area-based change statistics
    bu1 = basin_area(A1, G_BUILT_UP, area_subs)
    bu2 = basin_area(A2, G_BUILT_UP, area_subs)
    nat1 = basin_area(A1, G_NATURAL_SEMINAT, area_subs)
    nat2 = basin_area(A2, G_NATURAL_SEMINAT, area_subs)
    if bu1 <= 0 or nat1 <= 0:
        raise SystemExit("ERROR: zero baseline built-up or natural area at 2015 "
                         "(check CLASS_TO_GROUP mapping and area units).")
    builtup_exp = 100.0 * (bu2 - bu1) / bu1
    natural_contr = 100.0 * (nat1 - nat2) / nat1

    dt = YEAR2 - YEAR1
    ib1 = group_area(A1, IST_SUB, G_BUILT_UP)
    ib2 = group_area(A2, IST_SUB, G_BUILT_UP)
    if ib1 <= 0:
        raise SystemExit(f"ERROR: zero Istanbul built-up area at {YEAR1}.")
    ist_rate = 100.0 * (ib2 - ib1) / (ib1 * dt)
    ist_vo1 = group_area(A1, IST_SUB, G_VEGETATION_OPEN)
    ist_vo2 = group_area(A2, IST_SUB, G_VEGETATION_OPEN)
    ist_vegopen_loss = ist_vo1 - ist_vo2

    sa1 = group_area(A1, SOUTH_SUB, G_AGRI_WETLAND)
    sa2 = group_area(A2, SOUTH_SUB, G_AGRI_WETLAND)
    if sa1 <= 0:
        raise SystemExit(f"ERROR: zero Southern-Marmara agri/wetland area at {YEAR1}.")
    south_contr = 100.0 * (sa1 - sa2) / sa1

    # ---- ordered (token, value) lists matching the file's \ph{} sequence ----
    acc = []
    acc += [("XX.X", f_oa(min(oas.values()))), ("XX.X", f_oa(max(oas.values())))]
    acc += [("0.XX", f_kap(min(kaps.values()))), ("0.XX", f_kap(max(kaps.values())))]
    acc += [("0.XX", f_kap(min(f1s.values()))), ("0.XX", f_kap(max(f1s.values())))]
    acc += [("XX", disp[best]), ("XX.X", f_oa(oas[best])), ("0.XX", f_kap(kaps[best]))]
    acc += [("XX", disp[worst]), ("XX.X", f_oa(oas[worst]))]
    acc += [("XX", lab(a_ref)), ("XX", lab(b_ref))]
    acc += [("XX", str(pu_threshold))]
    for s in SUBREGION_ORDER:                         # table rows
        acc += [("XX.X", f_oa(oas[s])), ("0.XX", f_kap(kaps[s])), ("0.XX", f_kap(f1s[s]))]

    lc = [("XX", f_pct(builtup_exp)), ("XX", f_pct(natural_contr)),
          ("XX", f_km(ib1)), ("XX", f_km(ib2)), ("XX.X", f_rate(ist_rate)),
          ("XX", f_km(ist_vegopen_loss)), ("XX", f_pct(south_contr))]

    summary = dict(best=disp[best], worst=disp[worst],
                   oa_range=(f_oa(min(oas.values())), f_oa(max(oas.values()))),
                   confused=(lab(a_ref), lab(b_ref)), pu=pu_threshold,
                   builtup_exp=builtup_exp, natural_contr=natural_contr,
                   ist_bu=(ib1, ib2), ist_rate=ist_rate,
                   ist_vegopen_loss=ist_vegopen_loss, south_contr=south_contr,
                   per_sub={s: (oas[s], kaps[s], f1s[s]) for s in SUBREGION_ORDER})
    return acc, lc, summary


# ------------------------------ tex filling --------------------------------
def fill_block(block, expected, tag):
    it = iter(expected)
    used = [0]

    def repl(m):
        tok = m.group(1)
        try:
            exp_tok, val = next(it)
        except StopIteration:
            raise SystemExit(f"[{tag}] more \\ph{{}} in the file than computed "
                             "values -- the section wording changed; re-check "
                             "the placeholder order.")
        if tok != exp_tok:
            raise SystemExit(f"[{tag}] placeholder mismatch: file has \\ph{{{tok}}} "
                             f"where value '{val}' (expected token '{exp_tok}') was "
                             "due. The section's placeholder order changed.")
        used[0] += 1
        return val

    out, _ = re.subn(r"\\ph\{([^}]*)\}", repl, block)
    leftover = list(it)
    if leftover:
        raise SystemExit(f"[{tag}] {len(leftover)} computed values had no "
                         "placeholder to go into -- wording changed.")
    return out, used[0]


def locate(tex, start_marker, end_marker):
    i = tex.index(start_marker)
    j = tex.index(end_marker, i + len(start_marker))
    return i, j


def main():
    ap = argparse.ArgumentParser(description="Fill Results placeholders from GRASS outputs.")
    ap.add_argument("--tex", default="article_marmara_coasts.tex")
    ap.add_argument("--confusion")
    ap.add_argument("--areas2015")
    ap.add_argument("--areas2025")
    ap.add_argument("--area-units", choices=list(_UNIT), default="m2")
    ap.add_argument("--emit-blocks", default="results_filled.tex")
    ap.add_argument("--write-tex", default=None)
    ap.add_argument("--turkish", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    disp = SUBREGION_DISPLAY_TR if args.turkish else SUBREGION_DISPLAY

    if args.selftest:
        run_selftest(disp)
        return

    if not (args.confusion and args.areas2015 and args.areas2025):
        ap.error("need --confusion, --areas2015 and --areas2025 (or --selftest)")

    conf = load_confusion(args.confusion)
    A1 = load_areas(args.areas2015, args.area_units)
    A2 = load_areas(args.areas2025, args.area_units)
    area_subs = [s for s in A1 if s != "__all__"] or ["__all__"]

    acc, lc, summ = build_values(conf, A1, A2, disp, area_subs)

    tex = open(args.tex, encoding="utf-8").read()
    a0, a1 = locate(tex, r"\subsection{Classification Accuracy}",
                    r"\subsection{Land Cover Change")
    l0 = tex.index(r"\subsection{Land Cover Change and Urbanization (2015--2025)}")
    l1 = tex.index(r"\begin{figure}", l0)

    acc_block, na = fill_block(tex[a0:a1], acc, "Classification Accuracy")
    lc_block, nl = fill_block(tex[l0:l1], lc, "Land Cover Change")

    with open(args.emit_blocks, "w", encoding="utf-8") as fh:
        fh.write(acc_block.rstrip() + "\n\n" + lc_block.rstrip() + "\n")
    print(f"Filled {na} placeholders in Classification Accuracy and {nl} in Land "
          f"Cover Change -> {args.emit_blocks}")

    if args.write_tex:
        new = tex[:a0] + acc_block + tex[a1:l0] + lc_block + tex[l1:]
        open(args.write_tex, "w", encoding="utf-8").write(new)
        print(f"Wrote full filled manuscript -> {args.write_tex}")

    print_summary(summ)


def print_summary(s):
    print("\n--- computed values ---")
    print(f"  OA range           : {s['oa_range'][0]}-{s['oa_range'][1]} %")
    print(f"  best / worst region: {s['best']}  /  {s['worst']}")
    print(f"  most-confused pair : {s['confused'][0]} <-> {s['confused'][1]}")
    print(f"  water+built-up P/U : above {s['pu']} %")
    print(f"  built-up expansion : {s['builtup_exp']:.1f} %")
    print(f"  natural contraction: {s['natural_contr']:.1f} %")
    print(f"  Istanbul built-up  : {s['ist_bu'][0]:.1f} -> {s['ist_bu'][1]:.1f} km^2 "
          f"({s['ist_rate']:.2f} %/yr)")
    print(f"  Istanbul veg+open  : -{s['ist_vegopen_loss']:.1f} km^2")
    print(f"  S.Marmara agri+wet : -{s['south_contr']:.1f} %")
    print("  per sub-region (OA%, kappa, F1):")
    for k, (o, ka, f) in s["per_sub"].items():
        print(f"    {k:20s} {o:5.1f}  {ka:.2f}  {f:.2f}")


# ------------------------------- self test ---------------------------------
def run_selftest(disp):
    print("=" * 70)
    print("SELF-TEST -- SYNTHETIC numbers to check the arithmetic only.")
    print("These are NOT results and are never written to any .tex file.")
    print("=" * 70)
    import random
    random.seed(0)
    classes = list(CLASS_TO_GROUP)
    conf = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    A1 = defaultdict(lambda: defaultdict(float))
    A2 = defaultdict(lambda: defaultdict(float))
    for s in SUBREGION_ORDER:
        for c in classes:
            for r in classes:
                conf[s][c][r] = random.randint(40, 80) if c == r else random.randint(0, 5)
            A1[s][c] = random.uniform(20, 200)
            A2[s][c] = A1[s][c] * (1.15 if CLASS_TO_GROUP[c] == "built_up"
                                   else random.uniform(0.85, 1.0))
    acc, lc, summ = build_values(conf, A1, A2, disp,
                                 [s for s in SUBREGION_ORDER])
    print_summary(summ)
    # exercise the token machinery against the real file if present
    for texname in ("article_marmara_coasts.tex", "working.tex"):
        try:
            tex = open(texname, encoding="utf-8").read()
        except FileNotFoundError:
            continue
        a0, a1 = locate(tex, r"\subsection{Classification Accuracy}",
                        r"\subsection{Land Cover Change")
        l0 = tex.index(r"\subsection{Land Cover Change and Urbanization (2015--2025)}")
        l1 = tex.index(r"\begin{figure}", l0)
        _, na = fill_block(tex[a0:a1], acc, "acc")
        _, nl = fill_block(tex[l0:l1], lc, "lc")
        print(f"\nToken alignment vs {texname}: OK "
              f"({na} accuracy + {nl} land-cover placeholders matched).")
        break


if __name__ == "__main__":
    main()
