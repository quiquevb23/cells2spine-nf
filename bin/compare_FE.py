#!/usr/bin/env python3
"""
compare_FE.py

Compare Functional Enrichment (Reactome CSVs only) between Cell2location and CSIDE
for each area and cell type, for both GSEA and ORA results.

Outputs:
  - per-celltype Venn diagrams for "up" and "down" terms (if direction available)
  - per-celltype summary CSV (overlap counts + jaccard)
  - area-level UpSet plot combining all celltypes/methods
"""

import argparse
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from upsetplot import from_contents, plot as upset_plot

# -------------------------
def find_term_and_pcols(df):
    cols = [c.strip() for c in df.columns]
    term_col = next((c for c in cols if re.search(r"term|name|description", c, re.I)), None)
    pcol = next((c for c in cols if re.search(r"adj.*p|p.adjust|fdr|padj|adj.p|qvalue", c, re.I)), None)
    nes_col = next((c for c in cols if re.search(r"nes|normalizedenrich|enrichmentscore", c, re.I)), None)
    return term_col, pcol, nes_col

def read_reactome_terms(path, pval_cutoff=0.1):
    """
    Read a Reactome CSV and return dict with sets: {'up', 'down', 'both'}.
    - Uses NES (or enrichmentScore) sign for up/down if present.
    - Filters by adjusted p-value (if present) using pval_cutoff.
    - If no NES found, returns 'both' (no direction).
    """
    df = pd.read_csv(path)
    if df.empty:
        return {"up": set(), "down": set(), "both": set()}
    term_col, pcol, nes_col = find_term_and_pcols(df)
    if term_col is None:
        # give up: return empty
        return {"up": set(), "down": set(), "both": set()}

    # filter by p-value if available
    if pcol is not None:
        try:
            df = df[pd.to_numeric(df[pcol], errors="coerce") <= pval_cutoff]
        except Exception:
            pass

    terms = df[term_col].astype(str).tolist()

    # if NES-like column exists, split by sign
    if nes_col is not None:
        try:
            nes_vals = pd.to_numeric(df[nes_col], errors="coerce")
            up_terms = set(df.loc[nes_vals > 0, term_col].astype(str))
            down_terms = set(df.loc[nes_vals < 0, term_col].astype(str))
            both = set(terms)  # overall
            return {"up": up_terms, "down": down_terms, "both": both}
        except Exception:
            # fallback to undirected
            return {"up": set(), "down": set(), "both": set(terms)}
    else:
        # No direction info -> return everything under 'both'
        return {"up": set(), "down": set(), "both": set(terms)}

def extract_celltype_from_filename(filename, area):
    """
    Given filename like DEA_dorsal_gm_Astrocytes_limma_GSEA_Reactome.csv
    and area='dorsal_gm', extract 'Astrocytes' robustly (handles underscores in ct).
    """
    prefix = f"DEA_{area}_"
    if not filename.startswith(prefix):
        return None
    m = re.match(rf"^{re.escape(prefix)}(.+?)_.*Reactome\.csv$", filename)
    #m = re.match(rf"^{re.escape(prefix)}(.+?)_.*GO_BP\.csv$", filename)
    if m:
        return m.group(1)
    # fallback: try removing known suffix patterns
    s = filename[len(prefix):]
    s = re.sub(r"_limma_.*GO_BP\.csv$", "", s)
    s = re.sub(r"_.*GO_BP\.csv$", "", s)
    #s = re.sub(r"_limma_.*Reactome\.csv$", "", s)
    #s = re.sub(r"_.*Reactome\.csv$", "", s)
    return s or None

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

# -------------------------
def compare_area_analysis(area, analysis, cell2loc_base, rctd_base, out_base):
    """
    Compare a single analysis type (GSEA or ORA) for one area.
    Only Reactome CSVs considered.
    """
    whole_sample = (area == "WHOLE_SAMPLE")
    c2l_dir  = os.path.join(cell2loc_base, analysis) if whole_sample else os.path.join(cell2loc_base, analysis, area)
    rctd_dir = os.path.join(rctd_base, analysis)     if whole_sample else os.path.join(rctd_base, analysis, area)
    if not os.path.isdir(c2l_dir):
        print(f"[WARN] Missing Cell2location FE dir: {c2l_dir}")
        return []
    if not os.path.isdir(rctd_dir):
        print(f"[WARN] Missing CSIDE FE dir: {rctd_dir}")
        return []

    # gather cell types present in both (based on Reactome csv names)
    rctd_files = [f for f in os.listdir(rctd_dir) if f.endswith("GO_BP.csv") and
                  f.startswith(f"DEA_WHOLE_SAMPLE_" if whole_sample else f"DEA_{area}_")]
    cts = []
    for f in rctd_files:
        ct = extract_celltype_from_filename(f, area)
        if ct:
            # check matching file in cell2loc
            # search for a file in cell2loc that contains the ct and endswith Reactome.csv
            matches = [g for g in os.listdir(c2l_dir) if g.endswith("GO_BP.csv") and g.startswith(f"DEA_{area}_") and re.search(rf"_{re.escape(ct)}(_|$)", g)]
            if len(matches) > 0:
                cts.append(ct)
    cts = sorted(set(cts))

    summary_rows = []

    # per-ct comparisons
    for ct in cts:
        # find matching file paths
        rctd_match = next((os.path.join(rctd_dir, f) for f in os.listdir(rctd_dir)
                          if f.endswith("GO_BP.csv") and f.startswith(f"DEA_{area}_") and re.search(rf"_{re.escape(ct)}(_|$)", f)), None)
        c2l_match  = next((os.path.join(c2l_dir, f) for f in os.listdir(c2l_dir)
                          if f.endswith("GO_BP.csv") and f.startswith(f"DEA_{area}_") and re.search(rf"_{re.escape(ct)}(_|$)", f)), None)
        if not rctd_match or not c2l_match:
            continue

        sets_rctd = read_reactome_terms(rctd_match)
        sets_c2l  = read_reactome_terms(c2l_match)

        # prepare output dir per area/analysis/ct
        ct_out = os.path.join(out_base, analysis, area)
        ensure_dir(ct_out)

        # try to produce up/down venns if both methods provide directed info (NES present)
        have_direction = (len(sets_rctd["up"]) + len(sets_rctd["down"]) > 0) or (len(sets_c2l["up"]) + len(sets_c2l["down"]) > 0)

        if have_direction:
            # UP
            s1_up = sets_c2l["up"]
            s2_up = sets_rctd["up"]
            plt.figure(figsize=(4,4))
            venn2([s1_up, s2_up], set_labels=("Cell2location", "CSIDE"))
            plt.title(f"{analysis} | {area} | {ct} (GO_BP) — UP")
            plt.tight_layout()
            plt.savefig(os.path.join(ct_out, f"{ct}_{analysis}_GO_BP_UP_venn.png"), dpi=300)
            plt.close()

            # DOWN
            s1_dn = sets_c2l["down"]
            s2_dn = sets_rctd["down"]
            plt.figure(figsize=(4,4))
            venn2([s1_dn, s2_dn], set_labels=("Cell2location", "CSIDE"))
            plt.title(f"{analysis} | {area} | {ct} (GO_BP) — DOWN")
            plt.tight_layout()
            plt.savefig(os.path.join(ct_out, f"{ct}_{analysis}_GO_BP_DOWN_venn.png"), dpi=300)
            plt.close()
        else:
            # undirected / combined
            s1 = sets_c2l["both"]
            s2 = sets_rctd["both"]
            plt.figure(figsize=(4,4))
            venn2([s1, s2], set_labels=("Cell2location", "CSIDE"))
            plt.title(f"{analysis} | {area} | {ct} (GO_BP) — ALL")
            plt.tight_layout()
            plt.savefig(os.path.join(ct_out, f"{ct}_{analysis}_GO_BP_ALL_venn.png"), dpi=300)
            plt.close()

        # compute overlaps and jaccard for up/down/both
        def jaccard(a,b):
            if len(a|b)==0:
                return 0.0
            return len(a&b)/len(a|b)

        summary = {
            "analysis": analysis,
            "area": area,
            "celltype": ct,
            "n_c2l_up": len(sets_c2l["up"]),
            "n_c2l_down": len(sets_c2l["down"]),
            "n_c2l_all": len(sets_c2l["both"]),
            "n_rctd_up": len(sets_rctd["up"]),
            "n_rctd_down": len(sets_rctd["down"]),
            "n_rctd_all": len(sets_rctd["both"]),
            "up_overlap": len(sets_c2l["up"] & sets_rctd["up"]),
            "down_overlap": len(sets_c2l["down"] & sets_rctd["down"]),
            "all_overlap": len(sets_c2l["both"] & sets_rctd["both"]),
            "up_jaccard": jaccard(sets_c2l["up"], sets_rctd["up"]),
            "down_jaccard": jaccard(sets_c2l["down"], sets_rctd["down"]),
            "all_jaccard": jaccard(sets_c2l["both"], sets_rctd["both"]),
            "c2l_file": c2l_match,
            "rctd_file": rctd_match
        }
        summary_rows.append(summary)

    # area-level UpSet (combine all ct-method sets)
    if len(summary_rows) > 0:
        content = {}
        for row in summary_rows:
            area_dir = os.path.join(out_base, analysis, area)
            ct = row["celltype"]
            # read again to ensure we capture current sets (or reuse earlier if you refactor)
            # Here re-read
            rctd_match = row["rctd_file"]
            c2l_match  = row["c2l_file"]
            sets_r = read_reactome_terms(rctd_match)
            sets_c = read_reactome_terms(c2l_match)
            # include both up/down and all as separate items to show overlaps
            if len(sets_c["up"])>0 or len(sets_r["up"])>0:
                content[f"{ct}_C2L_up"] = sets_c["up"]
                content[f"{ct}_RCTD_up"] = sets_r["up"]
                content[f"{ct}_C2L_down"] = sets_c["down"]
                content[f"{ct}_RCTD_down"] = sets_r["down"]
            else:
                content[f"{ct}_C2L_all"] = sets_c["both"]
                content[f"{ct}_RCTD_all"] = sets_r["both"]

        if len(content) > 0:
            upset_data = from_contents(content)
            if upset_data.empty:
                print(f"[WARN] Skipping UpSet for {area} ({analysis}) — no data to plot.")
            else:
                ensure_dir(os.path.join(out_base, analysis, area))
                plt.figure(figsize=(10,5))
                upset_plot(upset_data)
                plt.title(f"{analysis} | {area} | GO_BP FE UpSet")
                plt.tight_layout()
                plt.savefig(os.path.join(out_base, analysis, area, f"{area}_{analysis}_GO_BP_upset.png"), dpi=300)
                plt.close()
        else:
            print(f"[WARN] Skipping UpSet for {area} ({analysis}) — empty content.")

    # return summary rows for aggregation
    return summary_rows

# -------------------------
def main(args):
    area = args.area
    cell2loc_base = args.cell2loc_dir
    rctd_base = args.rctd_dir
    out_base = args.output_dir
    ensure_dir(out_base)

    analyses = ["GSEA", "ORA"]

    all_summary = []
    for analysis in analyses:
        res = compare_area_analysis(area, analysis, cell2loc_base, rctd_base, out_base)
        if res:
            all_summary.extend(res)

    if all_summary:
        df = pd.DataFrame(all_summary)
        out_csv = os.path.join(out_base, f"FE_GO_BP_overlap_summary_{area}.csv")
        df.to_csv(out_csv, index=False)
        print(f"[INFO] Wrote FE summary -> {out_csv}")
    else:
        print(f"[WARN] No FE GO_BP comparisons produced for area {area}")

# -------------------------
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Compare Functional Enrichment GO_BP outputs (GSEA/ORA) between Cell2location and CSIDE")
    p.add_argument("--area", required=True)
    p.add_argument("--cell2loc_dir", required=True, help="base FE dir (contains GSEA/ and ORA/ subfolders)")
    p.add_argument("--rctd_dir", required=True, help="base FE dir for CSIDE")
    p.add_argument("--output_dir", required=True)
    args = p.parse_args()
    main(args)
