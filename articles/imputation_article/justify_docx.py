#!/usr/bin/env python3
# Post-process a pandoc-generated .docx in place:
#   1. Inject <w:jc w:val="both"/> into docDefaults so body paragraphs are
#      justified by default.
#   2. Fix the broken <w:tblW w:type="pct" w:w="0.0"/> that pandoc 2.9 emits
#      on every table. LibreOffice reads 0.0 as "no width" and renders the
#      table at an arbitrary huge width; replace with 5000 (= 100% in
#      fiftieths of a percent) so tables fit the page width in both
#      LibreOffice and Google Docs.
import os, re, sys, zipfile

JC = '<w:jc w:val="both"/>'

def patch_document(text: str) -> str:
    return re.sub(
        r'<w:tblW w:type="pct" w:w="0\.0" ?/>',
        '<w:tblW w:type="pct" w:w="5000"/>',
        text,
    )

def patch_styles(text: str) -> str:
    if 'w:jc w:val="both"' in text:
        return text  # already justified
    if re.search(r'<w:pPrDefault>\s*<w:pPr>', text):
        return re.sub(r'(<w:pPrDefault>\s*<w:pPr>)',
                      r'\1' + JC, text, count=1)
    if re.search(r'<w:pPrDefault\s*/>', text):
        return re.sub(r'<w:pPrDefault\s*/>',
                      f'<w:pPrDefault><w:pPr>{JC}</w:pPr></w:pPrDefault>',
                      text, count=1)
    if '<w:pPrDefault>' in text:
        return text.replace('<w:pPrDefault>',
                            f'<w:pPrDefault><w:pPr>{JC}</w:pPr>', 1)
    return text.replace('<w:docDefaults>',
        f'<w:docDefaults><w:pPrDefault><w:pPr>{JC}</w:pPr></w:pPrDefault>', 1)

def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: justify_docx.py FILE.docx")
    src = sys.argv[1]
    tmp = src + ".tmp"
    with zipfile.ZipFile(src, 'r') as zin:
        items = [(n, zin.read(n)) for n in zin.namelist()]
    with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zout:
        for name, data in items:
            if name == 'word/styles.xml':
                data = patch_styles(data.decode('utf-8')).encode('utf-8')
            elif name == 'word/document.xml':
                data = patch_document(data.decode('utf-8')).encode('utf-8')
            zout.writestr(name, data)
    os.replace(tmp, src)
    print(f"Justified body alignment + fixed table widths in: {src}")

if __name__ == '__main__':
    main()
