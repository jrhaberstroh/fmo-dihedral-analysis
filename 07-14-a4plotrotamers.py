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
#
resid_index = [rrid2resid_4BCL(i) for i in (rotamers[:,0]-1)]
resid_index[-1] = 366

print(rotamers.shape)
plt.plot(resid_index, rotamers[:,1], label="S2min")
plt.plot(resid_index, rotamers[:,2], label="S2max")
plt.title("Rotamers for 104ns FMO")
leg = plt.legend(loc=3, fancybox=True)
leg.get_frame().set_alpha(.5)
plt.xlabel("Residue Index")
plt.ylabel("S2 order parameter")

outname="rotamersS2_4BCL"
plt.tight_layout()
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-07-15/" + outname + ".png")
fig = plt.gcf()
fig.set_size_inches(20, 6, forward=True)
plt.tight_layout()
plt.xlim([-50,400])
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-07-15/" + outname + "_wide" + ".png")
plt.close()
matplotlib.rcParams.update(matplotlib.rcParamsDefault)

plt.plot(resid_index, rotamers[:,3], label="Phi", linewidth=3, alpha = .7)
plt.plot(resid_index, rotamers[:,4], label="Psi")
plt.plot(resid_index, rotamers[:,5], label="Omega")
plt.title("Rotamers for 104ns FMO")
leg = plt.legend(loc=3, fancybox=True)
leg.get_frame().set_alpha(.5)
plt.xlabel("Residue Index")
plt.ylabel("Mean Rotamer Cosine")

outname="rotamersPPO_4BCL"

plt.tight_layout()
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-07-15/" + outname + ".png")
fig = plt.gcf()
fig.set_size_inches(20, 6, forward=True)
plt.tight_layout()
plt.xlim([-50,400])
plt.savefig("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-07-15/" + outname + "_wide" + ".png")
plt.clf()


