#!/usr/bin/env bash
# Build main.docx from main.tex via pandoc.
#
# pandoc handles natbib \citep, \thebibliography (rendered inline),
# booktabs tables and most math out of the box. The one thing it can't
# render is the TikZ flowchart in section 2 — we work around that by
# rasterising the flowchart's PDF page from main.pdf and substituting
# an \includegraphics for the tikzpicture block in a temp .tex.
#
# Requires: pandoc, pdftoppm (poppler-utils), pdftotext, awk.
#
# Usage:
#   ./build_docx.sh           # writes main.docx alongside main.tex
#
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

TEX=main.tex
PDF=main.pdf
DOCX=main.docx
TMPTEX=.main_docx_src.tex
FIG=flowchart_for_docx.png

command -v pandoc   >/dev/null || { echo "pandoc not installed"   >&2; exit 1; }
command -v pdftoppm >/dev/null || { echo "pdftoppm not installed" >&2; exit 1; }
command -v pdftotext >/dev/null || { echo "pdftotext not installed" >&2; exit 1; }

if [[ ! -f "$PDF" || "$TEX" -nt "$PDF" ]]; then
  echo "Compiling $TEX -> $PDF (pdflatex x2 for refs)..."
  pdflatex -interaction=nonstopmode "$TEX" >/dev/null
  pdflatex -interaction=nonstopmode "$TEX" >/dev/null
fi

# Locate the page that contains the flowchart caption.
FIG_PAGE=
for p in $(seq 1 "$(pdfinfo "$PDF" | awk '/^Pages:/ {print $2}')"); do
  if pdftotext -layout -f "$p" -l "$p" "$PDF" - 2>/dev/null \
       | grep -q "Figure 1: Direct-merge"; then
    FIG_PAGE=$p
    break
  fi
done
if [[ -z "$FIG_PAGE" ]]; then
  echo "Could not locate flowchart page in $PDF" >&2
  exit 1
fi
echo "Flowchart found on page $FIG_PAGE"

# Render that page as PNG. pdftoppm names files with zero-padded page idx.
rm -f flowchart_tmp-*.png
pdftoppm -png -r 220 -f "$FIG_PAGE" -l "$FIG_PAGE" "$PDF" flowchart_tmp >/dev/null
mv flowchart_tmp-*.png "$FIG"

# Build a temp .tex with the tikzpicture block replaced by \includegraphics.
awk -v fig="$FIG" '
  /\\begin{tikzpicture}/ {
    in_tikz = 1
    print "\\includegraphics[width=\\linewidth]{" fig "}"
    next
  }
  /\\end{tikzpicture}/  { in_tikz = 0; next }
  !in_tikz              { print }
' "$TEX" > "$TMPTEX"

pandoc "$TMPTEX" \
  --from=latex \
  --to=docx \
  --output="$DOCX" \
  --standalone \
  --resource-path=.

rm -f "$TMPTEX"

# Post-process: set body paragraph alignment to justified (full width).
python3 "$SCRIPT_DIR/justify_docx.py" "$DOCX"

echo "Wrote: $(pwd)/$DOCX"
