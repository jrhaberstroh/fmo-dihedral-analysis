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


IN_FILE=SRCDIR+"/output/gyrate_4BCL.xvg"

r_gyration = jh_loadxvg(IN_FILE)


# backbone_index=np.arange(len(rmsf[:,1]))
# rrid_index = backbone_index / 3
# rrid_index = rrid_index.astype(int)
# print(rrid_index)
# resid_index = [rrid2resid_4BCL(i) for i in rrid_index]
# 
plt.plot(r_gyration[:,0], r_gyration[:,1])
plt.title("Radius of gyration for 104ns FMO")
plt.xlabel("Time, ns")
plt.ylabel("Radius of gyration, nm")
# 
outname="gyrate_4BCL.png"
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-07-15/" + outname, bbox_inches="tight")
plt.clf()
