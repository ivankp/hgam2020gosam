#!/usr/bin/env python3

import sys
if len(sys.argv)!=3:
    print("usage:",sys.argv[0],"in.json out.pdf")
    sys.exit(1)
print(sys.argv[1]+' -> '+sys.argv[2])

import math, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

with open(sys.argv[1],'r') as f:
    hf = json.load(f)

rew = "HT1 CT14nlo"

with PdfPages(sys.argv[2]) as out:
    for name, hist in hf[rew].items():
        print(name)

        yv = hist['xsec']
        yu = hist['mc_unc']
        xv = [ i+0.5 for i in range(len(yv)) ]
        xu = [ 0.5 for x in xv ]

        plt.figure(figsize=(6,4))
        plt.margins(x=0)

        plt.errorbar(xv, yv, yu, xu, ' ', linewidth=2, label=rew)

        plt.title(name)
        out.savefig(bbox_inches='tight')
        plt.close()

