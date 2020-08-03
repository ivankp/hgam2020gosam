#!/usr/bin/env python3

import sys
if len(sys.argv)<3:
    print("usage:",sys.argv[0],"in.cbor out.json [weights ...]")
    sys.exit(1)
print(sys.argv[1]+' -> '+sys.argv[2])

from cbor import cbor
import json, re, math
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

# for name, s in sets.items():
#     print(name)
#     print(s.scale)
#     print(s.pdf)

hf['bins'] = [ w for w,i in indep ] + list(sets.keys())

# print()
# make envelopes
for hname, h in hf["hists"].items():
    # print(hname)
    for hbin in h['bins'][1]:
        ws1 = hbin[0]
        ws2 = [ ws1[i] for w,i in indep ]
        for s in sets.values():
            w = ws1[s.scale[0]]
            ws2.append([
                w[0],
                math.sqrt(w[1]),
                [ min(ws1[i][0] for i in s.scale)-w[0],
                  max(ws1[i][0] for i in s.scale)-w[0] ],
                [ min(ws1[i][0] for i in s.pdf)-w[0],
                  max(ws1[i][0] for i in s.pdf)-w[0] ]
            ])
        hbin[0] = ws2

def subaxes(dim,n):
    for axis in dim:
        yield axis
    for i in range(n-len(dim)):
        yield dim[-1]

def overflow(axes):
    o = [ False ]
    for dim in axes:
        for i,axis in enumerate(subaxes(dim,len(o))):
            o[i] = [True] + [o[i]] * (len(axis)-1) + [True]
        o = [ item for sublist in o for item in sublist ]
    return o

def mask_list(xs,mask):
    for x,f in zip(xs,mask):
        if not f:
            yield x

# change structure & remove overflow bins
out = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for k,wname in enumerate(hf['bins']):
    if len(sys.argv)>3:
        if wname not in sys.argv[3:]:
            continue
    for hname, h in hf["hists"].items():
        hout = out[wname][hname]
        axes = hout['axes'] \
             = [ [ hf['axes'][i] for i in ii ] for ii in h['axes'] ]
        for y in mask_list( h['bins'][1], overflow(axes) ):
            y = y[0][k]
            # print(y)
            for y,yname in zip(y,('xsec','mc_unc','scale','pdf')):
                # print(yname,y)
                hout[yname].append(y)

with open(sys.argv[2],'w') as f:
    # json.dump(out,f,separators=(',',':'))
    jstr = json.dumps(out,indent=2)
    jstr = re.sub(r'(\d,|\[)\s+',r'\1',jstr)
    jstr = re.sub(r'\s+\]',r']',jstr)
    jstr = re.sub(r'\],\s+\[',r'],[',jstr)
    f.write(jstr)

