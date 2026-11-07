import sys
from pypdf import PdfReader
for p in sys.argv[1:]:
    try:
        r = PdfReader(p)
        t = "\n".join((pg.extract_text() or "") for pg in r.pages)
        open(p[:-4] + ".txt", "w").write(t)
        print(p, len(r.pages), "pages", len(t), "chars")
    except Exception as e:
        print(p, "ERR", e)
