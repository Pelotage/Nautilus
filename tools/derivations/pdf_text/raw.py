import re, zlib, sys
def txt(data):
    out=[]
    for m in re.finditer(rb'stream\r?\n', data):
        s=m.end()
        try:
            d=zlib.decompressobj().decompress(data[s:s+2000000])
        except Exception:
            continue
        if b'TJ' in d or b'Tj' in d:
            # crude text extraction
            parts=re.findall(rb'\[(.*?)\]\s*TJ|\((.*?)\)\s*Tj', d, re.S)
            line=[]
            for a,b in parts:
                seg=a or b
                line.append(b''.join(re.findall(rb'\(((?:\\.|[^\\)])*)\)', seg)) if a else seg)
            out.append(b' '.join(line))
    return b'\n'.join(out)
for p in sys.argv[1:]:
    t=txt(open(p,'rb').read())
    open(p+'.raw.txt','wb').write(t)
    print(p,len(t))
