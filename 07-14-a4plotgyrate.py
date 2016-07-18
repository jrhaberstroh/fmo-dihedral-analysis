#!/home/jhaberstroh/anaconda/bin/python
from __future__ import division, print_function
import numpy as np
import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
from jh_loadxvg import jh_loadxvg
from rrid2resid_4BCL import rrid2resid_4BCL
import os
SRCDIR=os.path.dirname(os.path.realpath(__file__))

raise NotImplementedError("This is a direct copy of RMSF code, not yet r_gyration code.")

IN_FILE=SRCDIR+"/output/rmsf.txt.xvg"

rmsf = jh_loadxvg(IN_FILE)

backbone_index=np.arange(len(rmsf[:,1]))
rrid_index = backbone_index / 3
rrid_index = rrid_index.astype(int)
print(rrid_index)
resid_index = [rrid2resid_4BCL(i) for i in rrid_index]

plt.plot(resid_index,rmsf[:,1], 'o')
plt.xlim([-7,366])
plt.xlabel("Residue Index")
plt.ylabel("rmsf, nm")
plt.title("RMSF of backbone atoms of 4BCL monomer")

plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-07-15/backbone_rmsf.png", bbox_inches="tight")
plt.clf()
