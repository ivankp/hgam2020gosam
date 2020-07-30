#!/usr/bin/env python3

import sys
if len(sys.argv)!=3:
    print("usage:",sys.argv[0],"in.cbor out.cbor")
    sys.exit(1)
print(sys.argv[1]+' -> '+sys.argv[2])

from cbor import cbor
import re
from collections import defaultdict

with open(sys.argv[1],'rb') as f:
    hf = cbor.load(f)

class scale_pdf:
    def __init__(self):
        self.scale = [ ]
        self.pdf = [ ]
    def __call__(self,i,pdf,ren,fac):
        if pdf==0:
            if (ren,fac)==(1,1):
                self.scale.insert(0,i)
                self.pdf.insert(0,i)
            else:
                self.scale.append(i)
        else:
            self.pdf.append(i)

sets = defaultdict(scale_pdf)
indep = [ ] # independent weights
weight_re = re.compile(r'([^:]+):(\S+)\s+ren:(\S+)\s+fac:(\S+)')
for i, weight in enumerate(hf["bins"][0][0]):
    m = weight_re.match(weight)
    if m is None:
        indep.append((weight,i))
        continue
    g = m.groups()
    # print(g)
    sets[g[0]](i,int(g[1]),float(g[2]),float(g[3]))

for name, s in sets.items():
    print(name)
    print(s.scale)
    print(s.pdf)

hf['bins'] = [ w for w,i in indep ] + list(sets.keys())

print()
for hname, h in hf["hists"].items():
    print(hname)
    for hbin in h['bins'][1]:
        ws1 = hbin[0]
        ws2 = [ ws1[i] for w,i in indep ]
        for s in sets.values():
            w = ws1[s.scale[0]]
            ws2.append([
                *w,
                [ min(ws1[i][0] for i in s.scale)-w[0],
                  max(ws1[i][0] for i in s.scale)-w[0] ],
                [ min(ws1[i][0] for i in s.pdf)-w[0],
                  max(ws1[i][0] for i in s.pdf)-w[0] ]
            ])
        hbin[0] = ws2

with open(sys.argv[2],'wb') as f:
    cbor.dump(hf,f)

