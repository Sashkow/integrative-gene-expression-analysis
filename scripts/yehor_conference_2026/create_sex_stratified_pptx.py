#!/usr/bin/env python3
"""
Build PPTX presentation for sex-stratified DE analysis (6ds run).

Usage:
    python3 scripts/yehor_conference_2026/create_sex_stratified_pptx.py
"""

import csv
import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR_TYPE

BASE = "articles/yehor_conference_2026/data/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
SEX_DIR = os.path.join(BASE, "sex_stratified")
VENN3_DIR = os.path.join(SEX_DIR, "venn_analysis")
VENN2_DIR = os.path.join(SEX_DIR, "venn_m_vs_f")
PLOTS_DIR = os.path.join(SEX_DIR, "plots")
DATA_DIR = "articles/yehor_conference_2026/data"
CROSSREF_PATH = os.path.join(DATA_DIR, "sex_specific_vs_prater.tsv")
OUTPUT = os.path.join(
    "articles/yehor_conference_2026",
    "sex_stratified_presentation.pptx",
)

DARK_BLUE = RGBColor(0x1A, 0x23, 0x7E)
MED_BLUE = RGBColor(0x30, 0x4F, 0xFE)
DARK_GREY = RGBColor(0x37, 0x47, 0x4F)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
ACCENT_GREEN = RGBColor(0x2E, 0x7D, 0x32)
ACCENT_RED = RGBColor(0xC6, 0x28, 0x28)
ACCENT_ORANGE = RGBColor(0xE6, 0x51, 0x00)
LIGHT_BG = RGBColor(0xF5, 0xF5, 0xF5)
LIGHT_GREY = RGBColor(0x9E, 0x9E, 0x9E)


def set_slide_bg(slide, color):
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_textbox(slide, left, top, width, height, text, font_size=14,
                bold=False, color=DARK_GREY, alignment=PP_ALIGN.LEFT,
                font_name="Calibri"):
    txbox = slide.shapes.add_textbox(left, top, width, height)
    tf = txbox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = bold
    p.font.color.rgb = color
    p.font.name = font_name
    p.alignment = alignment
    return txbox


def add_multiline_textbox(slide, left, top, width, height, lines,
                          font_name="Calibri"):
    """lines: list of (text, font_size, bold, color, alignment)"""
    txbox = slide.shapes.add_textbox(left, top, width, height)
    tf = txbox.text_frame
    tf.word_wrap = True
    for i, (text, font_size, bold, color, alignment) in enumerate(lines):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.text = text
        p.font.size = Pt(font_size)
        p.font.bold = bold
        p.font.color.rgb = color
        p.font.name = font_name
        p.alignment = alignment
        p.space_after = Pt(2)
    return txbox


def add_image_scaled(slide, img_path, left, top, max_width, max_height):
    from PIL import Image
    im = Image.open(img_path)
    w_px, h_px = im.size
    aspect = w_px / h_px
    target_aspect = max_width / max_height
    if aspect > target_aspect:
        width = max_width
        height = int(max_width / aspect)
    else:
        height = max_height
        width = int(max_height * aspect)
    slide.shapes.add_picture(img_path, left, top, width, height)


def add_table_cell(table, row, col, text, font_size=11, bold=False,
                   color=DARK_GREY, alignment=PP_ALIGN.CENTER,
                   font_name="Calibri"):
    cell = table.cell(row, col)
    cell.text = str(text)
    p = cell.text_frame.paragraphs[0]
    p.font.size = Pt(font_size)
    p.font.bold = bold
    p.font.color.rgb = color
    p.font.name = font_name
    p.alignment = alignment


# --------------- data loaders ---------------

def load_truly_sex_specific_female():
    path = os.path.join(VENN2_DIR, "group_female_only_annotated.tsv")
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    return [r for r in rows if r["near_miss_m"] == "FALSE"]


def load_truly_sex_specific_male():
    path = os.path.join(VENN2_DIR, "group_male_only_annotated.tsv")
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    return [r for r in rows if r["near_miss_f"] == "FALSE"]


def load_go_up():
    path = os.path.join(VENN2_DIR, "go_bp_female_truly_sex_specific_up.csv")
    with open(path) as f:
        reader = csv.DictReader(f)
        return list(reader)


def load_novel_candidates():
    with open(CROSSREF_PATH) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return [r for r in reader if r["in_prater"] == "FALSE"]


# --------------- slides ---------------

def make_title_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.8), Inches(1.2), Inches(11.7), Inches(1.5),
        "Improving Integrative Analysis of Human Placental Tissue "
        "Microarray Data Through Missing Value Imputation "
        "and Fetal Sex Prediction",
        font_size=28, bold=True, color=DARK_BLUE,
        alignment=PP_ALIGN.CENTER,
    )
    add_textbox(
        slide,
        Inches(0.8), Inches(3.0), Inches(11.7), Inches(0.6),
        "Yehor Poliakov, Oleksandr Lykhenko, Maria Obolenska",
        font_size=18, color=DARK_GREY,
        alignment=PP_ALIGN.CENTER,
    )
    add_textbox(
        slide,
        Inches(0.8), Inches(3.8), Inches(11.7), Inches(0.5),
        "Placental transcriptome: 1st vs 2nd trimester",
        font_size=16, color=LIGHT_GREY,
        alignment=PP_ALIGN.CENTER,
    )
    add_textbox(
        slide,
        Inches(0.8), Inches(6.5), Inches(11.7), Inches(0.5),
        "GSE100051 • GSE122214 • GSE28551 • GSE37901 • GSE93520 • GSE9984",
        font_size=12, color=LIGHT_GREY,
        alignment=PP_ALIGN.CENTER,
    )


def make_intro_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
        "Integrative Analysis",
        font_size=28, bold=True, color=DARK_BLUE,
        alignment=PP_ALIGN.CENTER,
    )

    left_lines = [
        ("What is integrative analysis?", 16, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Standardized preprocessing of raw data "
         "from multiple studies",
         13, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Direct merging into one unified "
         "expression matrix",
         13, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Unified statistical modeling across "
         "all samples",
         13, False, DARK_GREY, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(0.5), Inches(1.2), Inches(6.0), Inches(3.0),
        left_lines,
    )

    right_lines = [
        ("Why better than meta-analysis?", 16, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Meta-analysis:", 13, True, ACCENT_ORANGE, PP_ALIGN.LEFT),
        ("  combines summary statistics (p-values, effect sizes) "
         "from separate studies",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Integrative analysis:", 13, True, ACCENT_GREEN, PP_ALIGN.LEFT),
        ("  combines raw data → more statistical power",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("  can model sample-level covariates (sex, "
         "gestational age) across studies",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(6.8), Inches(1.2), Inches(6.0), Inches(3.5),
        right_lines,
    )

    challenge_lines = [
        ("Challenges we address:", 16, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Different microarray platforms measure different gene sets — "
         "intersection loses >50% of genes",
         13, False, ACCENT_ORANGE, PP_ALIGN.LEFT),
        ("  → softImpute imputation: 8,260 → 17,531 genes (2.1× more)",
         13, True, ACCENT_GREEN, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Fetal sex is absent from most GEO metadata",
         13, False, ACCENT_ORANGE, PP_ALIGN.LEFT),
        ("  → massiR prediction from Y-chromosome probe intensity",
         13, True, ACCENT_GREEN, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(0.5), Inches(4.5), Inches(12.0), Inches(3.0),
        challenge_lines,
    )


def make_dataset_table_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
        "6 GEO Microarray Datasets of Human Placenta",
        font_size=28, bold=True, color=DARK_BLUE,
        alignment=PP_ALIGN.CENTER,
    )

    rows_count = 8  # header + 6 datasets + total
    cols_count = 5
    tbl_left = Inches(0.5)
    tbl_top = Inches(1.2)
    tbl_width = Inches(12.3)
    tbl_height = Inches(3.5)

    shape = slide.shapes.add_table(
        rows_count, cols_count, tbl_left, tbl_top, tbl_width, tbl_height,
    )
    table = shape.table

    table.columns[0].width = Inches(1.8)
    table.columns[1].width = Inches(3.0)
    table.columns[2].width = Inches(2.8)
    table.columns[3].width = Inches(2.5)
    table.columns[4].width = Inches(2.2)

    headers = ["GEO ID", "Platform", "1T samples (M/F)",
               "2T samples (M/F)", "Genes"]
    for c, hdr in enumerate(headers):
        add_table_cell(table, 0, c, hdr, font_size=12, bold=True,
                       color=WHITE)
        cell = table.cell(0, c)
        cell.fill.solid()
        cell.fill.fore_color.rgb = DARK_BLUE

    dataset_rows = [
        ("GSE100051", "Illumina (GPL10558)",  "42 (16M, 26F)", "7 (2M, 5F)",   "11,761"),
        ("GSE122214", "Affymetrix (GPL570)",  "4 (2M, 2F)",   "—",            "16,925"),
        ("GSE28551",  "ABI (GPL2986)",        "16 (9M, 7F)",  "—",            "11,013"),
        ("GSE37901",  "Affymetrix (GPL570)",  "—",            "4 (3M, 1F)",   "16,925"),
        ("GSE93520",  "Agilent (GPL6480)",    "36 (17M, 19F)","—",            "12,055"),
        ("GSE9984",   "Affymetrix (GPL570)",  "4 (2M, 2F)",   "4 (1M, 3F)",   "16,925"),
    ]
    for r, row_data in enumerate(dataset_rows, start=1):
        for c, val in enumerate(row_data):
            add_table_cell(table, r, c, val, font_size=11)

    total_row = 7
    add_table_cell(table, total_row, 0, "Total", font_size=12, bold=True,
                   color=DARK_BLUE)
    add_table_cell(table, total_row, 1, "4 platforms, 3 manufacturers",
                   font_size=11, bold=True, color=DARK_BLUE)
    add_table_cell(table, total_row, 2, "102 (46M, 56F)",
                   font_size=11, bold=True, color=DARK_BLUE)
    add_table_cell(table, total_row, 3, "15 (6M, 9F)",
                   font_size=11, bold=True, color=DARK_BLUE)
    add_table_cell(table, total_row, 4, "117 samples",
                   font_size=11, bold=True, color=DARK_BLUE)
    for c in range(cols_count):
        cell = table.cell(total_row, c)
        cell.fill.solid()
        cell.fill.fore_color.rgb = RGBColor(0xE8, 0xEA, 0xF6)

    comparison_lines = [
        ("Intersection of all platforms:  8,260 genes",
         16, True, ACCENT_ORANGE, PP_ALIGN.LEFT),
        ("After softImpute imputation:  17,531 genes  (2.1× more)",
         16, True, ACCENT_GREEN, PP_ALIGN.LEFT),
        ("", 6, False, DARK_GREY, PP_ALIGN.LEFT),
        ("All M/F counts predicted by massiR "
         "(0 known fetal sex in GEO metadata)",
         12, False, LIGHT_GREY, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(0.5), Inches(5.2), Inches(12.0), Inches(2.0),
        comparison_lines,
    )


def make_methods_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
        "Materials and Methods",
        font_size=28, bold=True, color=DARK_BLUE,
        alignment=PP_ALIGN.CENTER,
    )

    left_lines = [
        ("Preprocessing Pipeline", 16, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Background correction, between-sample normalization, "
         "gene-level summarization per dataset",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Missing value imputation (softImpute — "
         "low-rank matrix factorization)",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Fetal sex prediction (massiR — "
         "Y-chromosome probe intensity)",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Batch effect correction "
         "(ComBat, reference: GSE100051)",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(0.5), Inches(1.1), Inches(6.0), Inches(5.5),
        left_lines,
    )

    right_lines = [
        ("Differential Expression", 16, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• limma (linear models for microarrays)",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Thresholds: |logFC| ≥ 1, adjusted p-value < 0.05",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 6, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Sex-Stratified Comparisons", 16, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• 1T vs 2T combined (sex as covariate)",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• 1T vs 2T males only",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• 1T vs 2T females only",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• M vs F within each trimester",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 6, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Near-Miss Analysis", 16, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Genes significant in one sex but nearly significant "
         "in the other (FDR < 0.10, |logFC| > 0.8)",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Distinguishes power artifacts from "
         "true sex-specific biology",
         12, False, ACCENT_ORANGE, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(6.8), Inches(1.1), Inches(6.0), Inches(5.5),
        right_lines,
    )


def make_volcano_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.5), Inches(0.15), Inches(12), Inches(0.5),
        "Volcano Plots: Sex-Stratified Comparisons",
        font_size=24, bold=True, color=DARK_BLUE,
    )

    volcanos = [
        ("volcano_males_1t_vs_2t.png", "Males: 1T vs 2T (337 DEGs)"),
        ("volcano_females_1t_vs_2t.png", "Females: 1T vs 2T (504 DEGs)"),
        ("volcano_1t_male_vs_female.png", "1T: Male vs Female (3 DEGs)"),
        ("volcano_2t_male_vs_female.png", "2T: Male vs Female (1 DEG)"),
    ]

    positions = [
        (Inches(0.15), Inches(0.7)),
        (Inches(6.35), Inches(0.7)),
        (Inches(0.15), Inches(3.9)),
        (Inches(6.35), Inches(3.9)),
    ]

    img_w = Inches(6.2)
    img_h = Inches(3.6)

    for (fname, label), (lf, tp) in zip(volcanos, positions):
        img_path = os.path.join(PLOTS_DIR, fname)
        add_image_scaled(slide, img_path, lf, tp, img_w, img_h)

    ref_lines = [
        ("Published M vs F comparisons:",
         10, True, ACCENT_GREEN, PP_ALIGN.LEFT),
        ("• Gonzalez 2018 (PMID 29335024): M vs F late T1, RNA-seq, "
         "58 sex-DEGs (25 X, 15 Y, 18 autosomal)",
         8, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Braun 2021 (PMID 33150487): M vs F at 11-16 wk (~T2), "
         "RNA-seq, 322 sex-DE transcripts",
         8, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• Flowers/Gonzalez 2024 (PMID 38537412): M vs F at T1 and T3, "
         "RNA-seq, 94 sex-DEGs at T1, 26 at T3",
         8, False, DARK_GREY, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(0.15), Inches(6.85), Inches(12.5), Inches(0.65),
        ref_lines,
    )


def _make_single_volcano_slide(prs, img_filename, title):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)
    img_path = os.path.join(PLOTS_DIR, img_filename)
    add_image_scaled(
        slide, img_path,
        Inches(0.3), Inches(0.15),
        Inches(12.7), Inches(7.2),
    )


def make_volcano_combined_slide(prs):
    _make_single_volcano_slide(
        prs, "volcano_combined_1t_vs_2t.png",
        "Combined: 1T vs 2T")


def make_volcano_males_slide(prs):
    _make_single_volcano_slide(
        prs, "volcano_males_1t_vs_2t_top10.png",
        "Males: 1T vs 2T")


def make_volcano_females_slide(prs):
    _make_single_volcano_slide(
        prs, "volcano_females_1t_vs_2t_top10.png",
        "Females: 1T vs 2T")


def make_stratification_advantage_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
        "Sex Stratification Reveals More DEGs",
        font_size=24, bold=True, color=DARK_BLUE,
    )

    img_path = os.path.join(PLOTS_DIR, "venn_1_2_sex_stratified.png")
    add_image_scaled(
        slide, img_path,
        Inches(0.5), Inches(1.0), Inches(8.0), Inches(6.0),
    )

    lines = [
        ("1T vs 2T DEG counts:", 14, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Combined (sex as covariate):", 13, False, DARK_GREY, PP_ALIGN.LEFT),
        ("  447 DEGs", 14, True, DARK_GREY, PP_ALIGN.LEFT),
        ("", 6, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Males only:", 13, False, DARK_GREY, PP_ALIGN.LEFT),
        ("  337 DEGs", 14, True, DARK_GREY, PP_ALIGN.LEFT),
        ("", 6, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Females only:", 13, False, DARK_GREY, PP_ALIGN.LEFT),
        ("  504 DEGs", 18, True, ACCENT_RED, PP_ALIGN.LEFT),
        ("", 8, False, DARK_GREY, PP_ALIGN.LEFT),
        ("504 > 447: sex stratification",
         13, True, ACCENT_GREEN, PP_ALIGN.LEFT),
        ("reveals additional genes",
         13, True, ACCENT_GREEN, PP_ALIGN.LEFT),
        ("", 8, False, DARK_GREY, PP_ALIGN.LEFT),
        ("317 core genes shared",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("across all three analyses",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(8.8), Inches(1.0), Inches(4.0), Inches(6.0),
        lines,
    )


def make_near_miss_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.3), Inches(0.15), Inches(12.5), Inches(0.5),
        "Near-Miss Analysis: Power Artifacts vs True Biology",
        font_size=24, bold=True, color=DARK_BLUE,
    )

    img_path = os.path.join(VENN2_DIR, "venn_m_vs_f.png")
    add_image_scaled(
        slide, img_path,
        Inches(0.3), Inches(0.8), Inches(7.0), Inches(3.5),
    )

    filter_lines = [
        ("Filtering results:", 14, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• 317 shared DEGs (100% direction concordant)",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• 20 male-only: 17 near-miss → ~3 truly sex-specific",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("• 187 female-only: 117 near-miss →",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("  ~70 truly sex-specific",
         13, True, ACCENT_RED, PP_ALIGN.LEFT),
        ("", 6, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Power context:", 12, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("Male 2T: n=6  |  Female 2T: n=9",
         11, False, DARK_GREY, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(7.5), Inches(0.8), Inches(5.3), Inches(3.5),
        filter_lines,
    )

    # Gene tables at bottom
    truly_f = load_truly_sex_specific_female()
    up_f = [g for g in truly_f if float(g["logFC_f"]) > 0]
    down_f = [g for g in truly_f if float(g["logFC_f"]) < 0]
    truly_m = load_truly_sex_specific_male()

    summary_lines = [
        (f"Female-only candidate sex-specific: 70 genes  "
         f"(↑ {len(up_f)} up  |  ↓ {len(down_f)} down in 2T)",
         12, True, ACCENT_RED, PP_ALIGN.LEFT),
        (f"Male-only candidate sex-specific: {len(truly_m)} genes",
         12, True, MED_BLUE, PP_ALIGN.LEFT),
        ("99.5% direction concordant — "
         "quantitative differences, not qualitative",
         12, True, ACCENT_ORANGE, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(0.3), Inches(4.4), Inches(12.5), Inches(0.9),
        summary_lines,
    )

    top_f = sorted(truly_f, key=lambda r: abs(float(r["logFC_f"])),
                   reverse=True)[:6]
    hdr_f = "Female-only  logFC♀   logFC♂   adjP♀       Dir"
    rows_f = [hdr_f]
    for g in top_f:
        lf = float(g["logFC_f"])
        lm = float(g["logFC_m"])
        ap = float(g["adjP_f"])
        arrow = "↑" if lf > 0 else "↓"
        rows_f.append(
            f"{g['SYMBOL']:12s} {lf:+.2f}    {lm:+.2f}    "
            f"{ap:.1e}   {arrow}")
    add_textbox(
        slide,
        Inches(0.3), Inches(5.3), Inches(6.0), Inches(2.0),
        "\n".join(rows_f),
        font_size=9, color=DARK_GREY, font_name="Consolas",
    )

    hdr_m = "Male-only    logFC♂   logFC♀   adjP♂       Dir"
    rows_m = [hdr_m]
    for g in truly_m:
        lm = float(g["logFC_m"])
        lf = float(g["logFC_f"])
        ap = float(g["adjP_m"])
        arrow = "↑" if lm > 0 else "↓"
        rows_m.append(
            f"{g['SYMBOL']:12s} {lm:+.2f}    {lf:+.2f}    "
            f"{ap:.1e}   {arrow}")
    add_textbox(
        slide,
        Inches(6.5), Inches(5.3), Inches(6.0), Inches(1.5),
        "\n".join(rows_m),
        font_size=9, color=MED_BLUE, font_name="Consolas",
    )


def make_go_dotplot_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.5), Inches(0.2), Inches(12), Inches(0.5),
        "GO Enrichment: Candidate Sex-Specific Genes (Upregulated)",
        font_size=24, bold=True, color=DARK_BLUE,
    )

    img_path = os.path.join(VENN2_DIR,
                            "go_bp_female_truly_sex_specific_up_dotplot.png")
    add_image_scaled(
        slide, img_path,
        Inches(0.5), Inches(0.8), Inches(8.5), Inches(6.5),
    )

    lines = [
        ("35 significant GO-BP terms", 13, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("(BH-corrected p < 0.05)", 10, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 6, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Dominant theme:", 12, True, ACCENT_RED, PP_ALIGN.LEFT),
        ("Vascular remodeling", 14, True, ACCENT_RED, PP_ALIGN.LEFT),
        ("+ oxygen sensing", 14, True, ACCENT_RED, PP_ALIGN.LEFT),
        ("", 6, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Universe: 17,531 genes", 10, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Test set: 36 up genes", 10, False, DARK_GREY, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(9.3), Inches(1.0), Inches(3.2), Inches(5.5),
        lines,
    )


def make_prater_comparison_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
        "Validation Against Prater 2021 RNA-seq",
        font_size=24, bold=True, color=DARK_BLUE,
    )

    img_path = os.path.join(DATA_DIR, "scatterplot_logfc_ours_vs_prater.png")
    add_image_scaled(
        slide, img_path,
        Inches(0.3), Inches(1.0), Inches(7.5), Inches(6.0),
    )

    lines = [
        ("Prater 2021 (PMID 34100896):", 14, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("RNA-seq, 7–8 wk vs 13–14 wk",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("n = 14, sex as covariate (not stratified)",
         12, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 6, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Shared DEGs with our combined analysis:",
         13, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("  333 DEGs, 100% direction concordant",
         13, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 8, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Of our 70 candidate sex-specific genes:",
         13, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("  21/70 also DE in Prater",
         13, False, DARK_GREY, PP_ALIGN.LEFT),
        ("  (detectable even when pooled)",
         11, False, LIGHT_GREY, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("  49/70 only detectable",
         14, True, ACCENT_RED, PP_ALIGN.LEFT),
        ("  with sex stratification",
         14, True, ACCENT_RED, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(8.0), Inches(1.0), Inches(4.8), Inches(6.0),
        lines,
    )


def make_novel_candidates_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.3), Inches(0.15), Inches(12.5), Inches(0.5),
        "49 Sex-Specific Candidates: Only Detectable with Stratification",
        font_size=24, bold=True, color=DARK_BLUE,
    )

    novel = load_novel_candidates()

    intro_lines = [
        ("These 49 genes are significant in female-only 1T vs 2T analysis",
         13, False, DARK_GREY, PP_ALIGN.LEFT),
        ("but NOT detected by Prater 2021 (combined, sex as covariate).",
         13, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 4, False, DARK_GREY, PP_ALIGN.LEFT),
        ("Without massiR sex prediction + sex-stratified DE,",
         13, True, ACCENT_ORANGE, PP_ALIGN.LEFT),
        ("these genes would be invisible.",
         13, True, ACCENT_ORANGE, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(0.3), Inches(0.7), Inches(12.5), Inches(1.6),
        intro_lines,
    )

    top_up = sorted(
        [g for g in novel if g["direction"] == "up"],
        key=lambda r: abs(float(r["logFC_f"])),
        reverse=True,
    )[:8]
    top_down = sorted(
        [g for g in novel if g["direction"] == "down"],
        key=lambda r: abs(float(r["logFC_f"])),
        reverse=True,
    )[:6]

    hdr = "Gene         logFC♀   logFC♂   adjP♀       Dir"
    rows_up = [f"Upregulated novel candidates ({sum(1 for g in novel if g['direction'] == 'up')} genes):", ""]
    rows_up.append(hdr)
    for g in top_up:
        lf = float(g["logFC_f"])
        lm = float(g["logFC_m"])
        ap = float(g["adjP_f"])
        rows_up.append(
            f"{g['SYMBOL']:12s} {lf:+.2f}    {lm:+.2f}    "
            f"{ap}   ↑")
    add_textbox(
        slide,
        Inches(0.3), Inches(2.4), Inches(6.0), Inches(4.5),
        "\n".join(rows_up),
        font_size=9, color=DARK_GREY, font_name="Consolas",
    )

    rows_down = [f"Downregulated novel candidates ({sum(1 for g in novel if g['direction'] == 'down')} genes):", ""]
    rows_down.append(hdr)
    for g in top_down:
        lf = float(g["logFC_f"])
        lm = float(g["logFC_m"])
        ap = float(g["adjP_f"])
        rows_down.append(
            f"{g['SYMBOL']:12s} {lf:+.2f}    {lm:+.2f}    "
            f"{ap}   ↓")
    add_textbox(
        slide,
        Inches(6.5), Inches(2.4), Inches(6.0), Inches(4.5),
        "\n".join(rows_down),
        font_size=9, color=DARK_GREY, font_name="Consolas",
    )

    go_note = [
        ("GO theme: many of these genes relate to vascular remodeling "
         "and oxygen sensing — consistent with slide 8 enrichment results",
         12, True, ACCENT_GREEN, PP_ALIGN.LEFT),
    ]
    add_multiline_textbox(
        slide,
        Inches(0.3), Inches(6.7), Inches(12.5), Inches(0.6),
        go_note,
    )


def make_conclusion_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
        "Conclusions",
        font_size=28, bold=True, color=DARK_BLUE,
        alignment=PP_ALIGN.CENTER,
    )

    conclusions = [
        ("1. Better annotation:",
         14, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("   softImpute expanded gene coverage from 8,260 "
         "(intersection) to 17,531 genes (2.1×)",
         14, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 8, False, DARK_GREY, PP_ALIGN.LEFT),
        ("2. Better integration:",
         14, True, DARK_BLUE, PP_ALIGN.LEFT),
        ("   massiR sex prediction enabled sex-stratified analysis, "
         "revealing 504 female DEGs vs 447 combined",
         14, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 8, False, DARK_GREY, PP_ALIGN.LEFT),
        ("3. 70 candidate sex-specific genes identified,",
         14, True, ACCENT_RED, PP_ALIGN.LEFT),
        ("   enriched in vascular remodeling and "
         "oxygen sensing pathways",
         14, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 8, False, DARK_GREY, PP_ALIGN.LEFT),
        ("4. 49/70 are novel — undetectable without "
         "sex stratification,",
         14, True, ACCENT_ORANGE, PP_ALIGN.LEFT),
        ("   absent from Prater 2021 RNA-seq comparison",
         14, False, DARK_GREY, PP_ALIGN.LEFT),
        ("", 8, False, DARK_GREY, PP_ALIGN.LEFT),
        ("5. Sex differences are quantitative "
         "(99.5% direction concordant), not qualitative —",
         14, False, DARK_GREY, PP_ALIGN.LEFT),
        ("   female placentas show stronger vascular response "
         "during 1T→2T transition",
         14, False, DARK_GREY, PP_ALIGN.LEFT),
    ]

    add_multiline_textbox(
        slide,
        Inches(0.8), Inches(1.2), Inches(11), Inches(5.5),
        conclusions,
    )


def _add_pipeline_box(slide, x, y, w, h, title, subtitle_lines, color,
                      title_size=16, sub_size=13):
    shape = slide.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE, x, y, w, h,
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    shape.shadow.inherit = False

    tf = shape.text_frame
    tf.word_wrap = True
    tf.auto_size = None

    p_title = tf.paragraphs[0]
    p_title.text = title
    p_title.font.size = Pt(title_size)
    p_title.font.bold = True
    p_title.font.color.rgb = WHITE
    p_title.font.name = "Calibri"
    p_title.alignment = PP_ALIGN.CENTER
    p_title.space_before = Pt(0)
    p_title.space_after = Pt(8)

    for line in subtitle_lines:
        p_sub = tf.add_paragraph()
        p_sub.text = line
        p_sub.font.size = Pt(sub_size)
        p_sub.font.bold = False
        p_sub.font.color.rgb = WHITE
        p_sub.font.name = "Calibri"
        p_sub.alignment = PP_ALIGN.CENTER
        p_sub.space_before = Pt(0)
        p_sub.space_after = Pt(0)

    return shape


def _add_pipeline_arrow(slide, x, y, w, h):
    arrow = slide.shapes.add_shape(
        MSO_SHAPE.NOTCHED_RIGHT_ARROW, x, y, w, h,
    )
    arrow.fill.solid()
    arrow.fill.fore_color.rgb = LIGHT_GREY
    arrow.line.fill.background()
    return arrow


def _add_down_arrow(slide, x, y, w, h):
    arrow = slide.shapes.add_shape(
        MSO_SHAPE.DOWN_ARROW, x, y, w, h,
    )
    arrow.fill.solid()
    arrow.fill.fore_color.rgb = LIGHT_GREY
    arrow.line.fill.background()
    return arrow


def _add_diagonal_arrow_connector(slide, x1, y1, x2, y2):
    from lxml import etree
    cnx = slide.shapes.add_connector(
        MSO_CONNECTOR_TYPE.STRAIGHT, x1, y1, x2, y2,
    )
    cnx.line.color.rgb = LIGHT_GREY
    cnx.line.width = Pt(2)
    ln = cnx._element.find(
        './/{http://schemas.openxmlformats.org/drawingml/2006/main}ln')
    if ln is not None:
        tail = etree.SubElement(
            ln,
            '{http://schemas.openxmlformats.org/drawingml/2006/main}tailEnd')
        tail.set('type', 'triangle')
        tail.set('w', 'lg')
        tail.set('len', 'lg')
    return cnx


def make_pipeline_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide,
        Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
        "Analysis Pipeline",
        font_size=28, bold=True, color=DARK_BLUE,
        alignment=PP_ALIGN.CENTER,
    )

    steps = [
        ("6 GEO Datasets", ["4 platforms", "117 samples"], DARK_BLUE),
        ("Incomplete Matrix", ["17,531 genes", "29% values missing"],
         DARK_GREY),
        ("Complete Matrix", ["softImpute", "17,531 × 117"], ACCENT_GREEN),
        ("Normalized Matrix", ["ComBat", "ref: GSE100051"], MED_BLUE),
    ]

    arrows = [
        "Merge by gene ID",
        "Imputation",
        "Batch correction",
        "limma\n(imputed excluded)",
    ]

    box_w = Inches(1.85)
    box_h = Inches(1.5)
    arrow_w = Inches(1.05)
    arrow_h = Inches(0.3)

    result_steps = [
        ("Combined DEGs", ["1T vs 2T: 447"], RGBColor(0x63, 0x49, 0x42)),
        ("Female DEGs", ["1T vs 2T ♀: 504"], ACCENT_RED),
        ("Male DEGs", ["1T vs 2T ♂: 337"], RGBColor(0x00, 0x69, 0x5C)),
    ]
    result_h = Inches(1.1)
    result_gap = Inches(0.15)
    result_block_h = 3 * result_h + 2 * result_gap

    n_main = len(steps)
    total_w = n_main * box_w + n_main * arrow_w + box_w
    slide_w = Inches(13.333)
    start_x = int((slide_w - total_w) / 2)
    row_center_y = Inches(3.75)
    box_y = int(row_center_y - box_h / 2)
    arrow_y = int(row_center_y - arrow_h / 2)

    caption_h = Inches(0.35)
    caption_y = box_y - Inches(0.85)
    caption_w = box_w + arrow_w

    for i, (title, subs, color) in enumerate(steps):
        x = start_x + i * (box_w + arrow_w)
        _add_pipeline_box(slide, x, box_y, box_w, box_h,
                          title, subs, color)

        if i < n_main - 1:
            ax = x + box_w
            _add_pipeline_arrow(slide, ax, arrow_y, arrow_w, arrow_h)
            cx = ax - box_w // 2
            add_textbox(slide, cx, caption_y, caption_w, caption_h,
                        arrows[i], font_size=16, color=DARK_GREY,
                        alignment=PP_ALIGN.CENTER)

    # --- 3 diagonal arrows from last box to each result box ---
    last_box_x = start_x + (n_main - 1) * (box_w + arrow_w)
    result_x = last_box_x + box_w + arrow_w
    result_top_y = int(row_center_y - result_block_h / 2)

    start_arrow_x = last_box_x + box_w
    start_arrow_y = int(box_y + box_h // 2)

    add_textbox(slide, last_box_x + box_w // 2, caption_y,
                caption_w, caption_h,
                arrows[-1], font_size=16, color=DARK_GREY,
                alignment=PP_ALIGN.CENTER)

    for i, (title, subs, color) in enumerate(result_steps):
        ry = result_top_y + i * (result_h + result_gap)
        end_arrow_y = int(ry + result_h // 2)
        _add_diagonal_arrow_connector(
            slide, start_arrow_x, start_arrow_y, result_x, end_arrow_y)
        _add_pipeline_box(slide, result_x, ry, box_w, result_h,
                          title, subs, color,
                          title_size=16, sub_size=16)


def make_pipeline_visual_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    set_slide_bg(slide, WHITE)

    add_textbox(
        slide, Inches(0.5), Inches(0.3), Inches(12), Inches(0.6),
        "Analysis Pipeline",
        font_size=28, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER,
    )

    # Per-dataset: (name, genes, samples) — sorted tallest to shortest
    ds_info = sorted([
        ("GSE122214", 16925, 4),
        ("GSE37901",  16925, 4),
        ("GSE9984",   16925, 8),
        ("GSE93520",  12055, 36),
        ("GSE100051", 11761, 49),
        ("GSE28551",  11013, 16),
    ], key=lambda x: -x[1])
    total_genes = 17531
    total_samples = 117
    gene_fracs = [g / total_genes for _, g, _ in ds_info]
    tallest_frac = max(gene_fracs)
    tallest_h_func = lambda mh: int(mh * tallest_frac)

    ds_colors = [
        RGBColor(0xE6, 0x9F, 0x00),
        RGBColor(0x56, 0xB4, 0xE9),
        RGBColor(0x00, 0x9E, 0x73),
        RGBColor(0xCC, 0x79, 0xA7),
        RGBColor(0x00, 0x72, 0xB2),
        RGBColor(0xD5, 0x5E, 0x00),
    ]
    IMPUTED_CLR = RGBColor(0xFF, 0xCC, 0x80)
    LINE_CLR = RGBColor(0xE0, 0xE0, 0xE0)

    # Column widths proportional to sample counts
    max_h = Inches(2.4)
    tallest_h = tallest_h_func(max_h)
    merged_total_w = Inches(1.8)
    sample_counts = [s for _, _, s in ds_info]
    col_widths = [int(merged_total_w * s / total_samples)
                  for s in sample_counts]

    sep_gap = Inches(0.04)
    block1_w = sum(col_widths) + 5 * sep_gap
    block_w = sum(col_widths)

    gap_w = Inches(0.70)
    arrow_shape_w = Inches(0.45)
    arrow_shape_h = Inches(0.28)

    result_h = Inches(1.0)
    result_v_gap = Inches(0.12)
    result_total_h = 3 * result_h + 2 * result_v_gap

    total_slide_content = (block1_w + 4 * gap_w + 3 * block_w
                           + block_w)
    slide_w = Inches(13.333)
    sx = int((slide_w - total_slide_content) / 2)
    top_y = Inches(1.7)
    mid_y = top_y + tallest_h // 2

    b = [sx]
    b.append(b[0] + block1_w + gap_w)
    b.append(b[1] + block_w + gap_w)
    b.append(b[2] + block_w + gap_w)
    rx = b[3] + block_w + gap_w

    def col_rect(cx, cy, w, h, color):
        s = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, cx, cy, w, h)
        s.fill.solid()
        s.fill.fore_color.rgb = color
        s.line.color.rgb = LINE_CLR
        s.line.width = Pt(0.5)

    # --- BLOCK 1: 6 separate datasets (proportional sizes) ---
    cx = b[0]
    for i in range(6):
        cw = col_widths[i]
        ch = int(max_h * gene_fracs[i])
        col_rect(cx, top_y, cw, ch, ds_colors[i])

        cx += cw + sep_gap

    # --- BLOCK 2: Merged incomplete (aligned top, bounding box) ---
    cx = b[1]
    for i in range(6):
        cw = col_widths[i]
        ch = int(max_h * gene_fracs[i])
        col_rect(cx, top_y, cw, ch, ds_colors[i])

        cx += cw
    bb = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, b[1], top_y, block_w, tallest_h)
    bb.fill.background()
    bb.line.color.rgb = DARK_GREY
    bb.line.width = Pt(1.5)

    # --- BLOCK 3: Complete (imputed areas filled) ---
    cx = b[2]
    for i in range(6):
        cw = col_widths[i]
        ch = int(max_h * gene_fracs[i])
        col_rect(cx, top_y, cw, ch, ds_colors[i])

        if ch < tallest_h:
            miss_y = top_y + ch
            miss_h = tallest_h - ch
            col_rect(cx, miss_y, cw, miss_h, IMPUTED_CLR)
        cx += cw

    # --- BLOCK 4: Normalized (original data MED_BLUE, imputed area pale blue) ---
    PALE_BLUE = RGBColor(0xBB, 0xDE, 0xFB)
    cx = b[3]
    for i in range(6):
        cw = col_widths[i]
        ch = int(max_h * gene_fracs[i])
        col_rect(cx, top_y, cw, ch, MED_BLUE)
        if ch < tallest_h:
            col_rect(cx, top_y + ch, cw, tallest_h - ch, PALE_BLUE)
        cx += cw

    # --- Arrows between blocks 1-2, 2-3, 3-4 + captions above arrows ---
    arrow_y = int(mid_y - arrow_shape_h // 2)
    block_ends = [b[0] + block1_w, b[1] + block_w, b[2] + block_w]
    arrow_labels = ["Merge by\ngene ID", "softImpute", "ComBat"]
    caption_h = Inches(0.35)
    caption_y = top_y - Inches(0.80)
    caption_extra_w = Inches(0.6)

    for j in range(3):
        ax = block_ends[j] + (gap_w - arrow_shape_w) // 2
        _add_pipeline_arrow(slide, ax, arrow_y,
                            arrow_shape_w, arrow_shape_h)
        cx = block_ends[j] - caption_extra_w // 2
        add_textbox(slide, cx, caption_y,
                    gap_w + caption_extra_w, caption_h,
                    arrow_labels[j],
                    font_size=16, color=DARK_GREY,
                    alignment=PP_ALIGN.CENTER)

    # --- 3 diagonal arrows from block 4 to DEG result boxes ---
    result_top = int(mid_y - result_total_h // 2)
    last_end = b[3] + block_w

    result_steps = [
        ("Combined DEGs", ["1T vs 2T: 447"], RGBColor(0x63, 0x49, 0x42)),
        ("Female DEGs", ["1T vs 2T ♀: 504"], ACCENT_RED),
        ("Male DEGs", ["1T vs 2T ♂: 337"], RGBColor(0x00, 0x69, 0x5C)),
    ]

    start_arrow_x = last_end
    start_arrow_y = int(mid_y)

    add_textbox(slide, last_end - caption_extra_w // 2, caption_y,
                gap_w + caption_extra_w, caption_h,
                "limma\n(imputed excl.)",
                font_size=16, color=DARK_GREY,
                alignment=PP_ALIGN.CENTER)

    for i, (title, subs, color) in enumerate(result_steps):
        ry = result_top + i * (result_h + result_v_gap)
        end_arrow_y = int(ry + result_h // 2)
        _add_diagonal_arrow_connector(
            slide, start_arrow_x, start_arrow_y, rx, end_arrow_y)
        _add_pipeline_box(slide, rx, ry, block_w, result_h,
                          title, subs, color,
                          title_size=14, sub_size=14)

    # --- Labels below blocks ---
    lbl_y = top_y + tallest_h + Inches(0.15)
    for lx, lw, txt in [
        (b[0], block1_w,
         "6 GEO Datasets\n4 platforms, 117 samples"),
        (b[1], block_w, "Incomplete Matrix\n17,531 genes"),
        (b[2], block_w, "Complete Matrix\n29% imputed"),
        (b[3], block_w, "Normalized Matrix"),
    ]:
        add_textbox(slide, lx, lbl_y, lw, Inches(0.6),
                    txt, font_size=12, color=DARK_GREY,
                    alignment=PP_ALIGN.CENTER)


def main():
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    make_title_slide(prs)                    # 1
    make_intro_slide(prs)                    # 2
    make_dataset_table_slide(prs)            # 3
    make_methods_slide(prs)                  # 4
    make_pipeline_slide(prs)                 # 5
    make_pipeline_visual_slide(prs)          # 6
    make_volcano_slide(prs)                  # 7
    make_volcano_combined_slide(prs)         # 8
    make_volcano_males_slide(prs)            # 9
    make_volcano_females_slide(prs)          # 10
    make_stratification_advantage_slide(prs) # 11
    make_near_miss_slide(prs)                # 12
    make_go_dotplot_slide(prs)               # 13
    make_prater_comparison_slide(prs)        # 14
    make_novel_candidates_slide(prs)         # 15
    make_conclusion_slide(prs)               # 12

    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    prs.save(OUTPUT)
    print(f"Saved: {OUTPUT}")
    print(f"Slides: {len(prs.slides)}")


if __name__ == "__main__":
    main()
