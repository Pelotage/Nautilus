import re, zlib, sys
data = open(sys.argv[1],'rb').read()
out=[]
for m in re.finditer(rb'stream\r?\n', data):
    start=m.end()
    end=data.find(b'endstream', start)
    chunk=data[start:end]
    try:
        d=zlib.decompress(chunk)
    except Exception:
        try: d=zlib.decompressobj().decompress(chunk)
        except Exception: continue
    if b'BT' not in d: continue
    txt=[]
    for bt in re.findall(rb'BT(.*?)ET', d, re.S):
        for s in re.findall(rb'\((?:\\.|[^\\)])*\)|\[(?:.*?)\]\s*TJ|<[0-9A-Fa-f]+>', bt, re.S):
            if s.startswith(b'['):
                parts=re.findall(rb'\(((?:\\.|[^\\)])*)\)', s)
                txt.append(b''.join(parts))
            elif s.startswith(b'('):
                txt.append(s[1:-1])
        txt.append(b'\n')
    out.append(b''.join(txt))
sys.stdout.buffer.write(b'\n'.join(out))
