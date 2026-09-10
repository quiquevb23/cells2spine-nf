#!/usr/bin/env python3
"""Add stub: blocks to Nextflow modules that lack them.
Run from the cells2spine-nf directory:
    python3 add_stubs.py
"""
import sys
from pathlib import Path

STUBS = {
    "modules/cell2loc_owndata.nf": """\
    stub:
    \"\"\"
    mkdir -p cell2location_map/stub
    touch cell2location_map/stub/placeholder
    \"\"\"
""",
    "modules/compare_cell2loc_rctd.nf": """\
    stub:
    \"\"\"
    mkdir -p DA Gene_Expr/WHOLE_SAMPLE DEGs/WHOLE_SAMPLE FEs/WHOLE_SAMPLE FEs/global_heatmaps
    touch DA/placeholder
    \"\"\"
""",
    "modules/cside.nf": """\
    stub:
    \"\"\"
    mkdir -p RCTD/stub CSIDE/stub
    touch RCTD/stub/placeholder CSIDE/stub/placeholder
    \"\"\"
""",
    "modules/ct_gene_expr_enrichment.nf": """\
    stub:
    \"\"\"
    mkdir -p FunctionalEnrichment/GSEA FunctionalEnrichment/ORA
    touch FunctionalEnrichment/GSEA/placeholder FunctionalEnrichment/ORA/placeholder
    \"\"\"
""",
    "modules/ct_gene_expr_percelltype.nf": """\
    stub:
    \"\"\"
    mkdir -p CT_Gene_expr/stub
    touch CT_Gene_expr/stub/placeholder
    \"\"\"
""",
    "modules/deg_area_enrichment.nf": """\
    stub:
    \"\"\"
    mkdir -p FunctionalEnrichment/GSEA FunctionalEnrichment/ORA
    touch FunctionalEnrichment/GSEA/placeholder FunctionalEnrichment/ORA/placeholder
    \"\"\"
""",
    "modules/deg_area_pipeline.nf": """\
    stub:
    \"\"\"
    mkdir -p DEA_areas_good/stub
    touch DEA_areas_good/stub/placeholder
    \"\"\"
""",
    "modules/extract_spatial_inputs.nf": """\
    stub:
    \"\"\"
    mkdir -p counts coords
    touch counts/placeholder coords/placeholder
    \"\"\"
""",
    "modules/masigpro.nf": """\
    stub:
    \"\"\"
    mkdir -p masigpro_results/stub
    touch masigpro_results/stub/placeholder
    \"\"\"
""",
}


def patch(path_str, stub_block):
    p = Path(path_str)
    if not p.exists():
        print(f"[MISSING] {path_str}")
        return

    text = p.read_text()

    if "stub:" in text:
        print(f"[SKIP]    {path_str}  (already has stub:)")
        return

    # Insert stub block before the final closing brace of the process
    lines = text.rstrip("\n").split("\n")
    insert_at = None
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].strip() == "}":
            insert_at = i
            break

    if insert_at is None:
        print(f"[ERROR]   {path_str}  — could not find closing }}", file=sys.stderr)
        return

    stub_lines = stub_block.rstrip("\n").split("\n")
    lines[insert_at:insert_at] = [""] + stub_lines
    p.write_text("\n".join(lines) + "\n")
    print(f"[OK]      {path_str}")


if __name__ == "__main__":
    for path, stub in STUBS.items():
        patch(path, stub)
    print("\nDone.")

