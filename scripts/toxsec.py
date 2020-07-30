#!/usr/bin/env python3

import sys
if len(sys.argv) not in (2,3):
    print("usage:",sys.argv[0],"in.cbor [out.cbor]")
    sys.exit(1)
print(sys.argv[1]+' -> '+sys.argv[-1])

from cbor import cbor

with open(sys.argv[1],'rb') as f:
    hf = cbor.load(f)

for hname, h in hf["hists"].items():
    print(hname)
    for ws,n,nent in h['bins'][1]:
        if n>1:
            for w in ws:
                w[1] = (w[1] - w[0]/n)/(n*(n-1)) # 1310.7439 p.16
                w[0] /= n

with open(sys.argv[-1],'wb') as f:
    cbor.dump(hf,f)

