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


IN_FILE=SRCDIR+"/output/rotamers_4BCL.xvg"
rotamers = jh_loadxvg(IN_FILE)


# backbone_index=np.arange(len(rmsf[:,1]))
# rrid_index = backbone_index / 3
# rrid_index = rrid_index.astype(int)
# print(rrid_index)
# resid_index = [rrid2resid_4BCL(i) for i in rrid_index]
#
print(rotamers.shape)
plt.plot(rotamers[:,0], rotamers[:,1], label="S2min")
plt.plot(rotamers[:,0], rotamers[:,2], label="S2max")
plt.title("Rotamers for 104ns FMO")
leg = plt.legend(loc=3, fancybox=True)
leg.get_frame().set_alpha(.5)
plt.xlabel("Time, ns")

outname="rotamersS2_4BCL.png"
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-07-15/" + outname, bbox_inches="tight")
plt.clf()

plt.plot(rotamers[:,0], rotamers[:,3], label="Phi", linewidth=3, alpha = .7)
plt.plot(rotamers[:,0], rotamers[:,4], label="Psi")
plt.plot(rotamers[:,0], rotamers[:,5], label="Omega")
plt.title("Rotamers for 104ns FMO")
leg = plt.legend(loc=3, fancybox=True)
leg.get_frame().set_alpha(.5)
plt.xlabel("Time, ns")

outname="rotamersPPO_4BCL.png"
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-07-15/" + outname, bbox_inches="tight")
plt.clf()

