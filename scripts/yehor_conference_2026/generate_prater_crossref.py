#!/usr/bin/env python3
"""
Generate cross-reference TSV of truly sex-specific genes vs Prater 2021 DEGs.

Usage:
    python3 scripts/yehor_conference_2026/generate_prater_crossref.py
"""

import csv
import openpyxl

BASE = "articles/yehor_conference_2026/data/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
SEX_DIR = f"{BASE}/sex_stratified"
VENN_DIR = f"{SEX_DIR}/venn_m_vs_f"
PRATER_PATH = "articles/yehor_conference_2026/references/Prater_2021_RNA-Seq_first-second_trimester_transition_supp_tables.xlsx"
OUTPUT = "articles/yehor_conference_2026/data/sex_specific_vs_prater.tsv"


def load_prater():
    wb = openpyxl.load_workbook(PRATER_PATH, read_only=True)
    ws = wb["T1 DEGs_results_table_l2fc1"]
    rows = list(ws.iter_rows(values_only=True))
    result = {}
    for row in rows[1:]:
        ensembl, baseMean, l2fc, lfcSE, stat, pval, padj, gene_name, entrez, desc = row
        if entrez is not None and str(entrez) != "NA":
            try:
                result[str(int(entrez))] = {
                    "symbol": gene_name,
                    "log2FC": float(l2fc),
                    "padj": float(padj) if padj and str(padj) != "NA" else None,
                }
            except (ValueError, TypeError):
                pass
    wb.close()
    return result


def load_truly_sex_specific():
    path = f"{VENN_DIR}/group_female_only_annotated.tsv"
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return [r for r in reader if r["near_miss_m"] == "FALSE"]


def main():
    prater = load_prater()
    truly = load_truly_sex_specific()

    with open(OUTPUT, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "ENTREZID", "SYMBOL", "logFC_f", "logFC_m", "adjP_f",
            "direction", "in_prater", "prater_log2FC", "prater_padj",
        ])
        in_count = 0
        out_count = 0
        for g in truly:
            eid = g["ENTREZID"]
            lf = float(g["logFC_f"])
            direction = "up" if lf > 0 else "down"
            in_prater = eid in prater
            row = [
                eid,
                g["SYMBOL"],
                f"{lf:.4f}",
                f"{float(g['logFC_m']):.4f}",
                f"{float(g['adjP_f']):.2e}",
                direction,
                "TRUE" if in_prater else "FALSE",
                f"{prater[eid]['log2FC']:.4f}" if in_prater else "",
                f"{prater[eid]['padj']:.2e}" if in_prater and prater[eid]["padj"] is not None else "",
            ]
            writer.writerow(row)
            if in_prater:
                in_count += 1
            else:
                out_count += 1

    print(f"Saved: {OUTPUT}")
    print(f"Total truly sex-specific: {len(truly)}")
    print(f"In Prater: {in_count}, NOT in Prater: {out_count}")


if __name__ == "__main__":
    main()
