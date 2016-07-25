from __future__ import division
from __future__ import print_function
import numpy as np
import argparse
import string
import logging
from jh_hist_fit import jh_hist_fit
from rrid2resid_4BCL import rrid2resid_4BCL

parser = argparse.ArgumentParser(description="")
parser.add_argument("-f", required=True, type=str, help="csv file with spaces for delimiters for energy gap")
parser.add_argument("-T", action="store_true", help="display timeseries")
parser.add_argument("-H", action="store_true", help="display histograms")
parser.add_argument("-c", type=int, default=0, help="Chromophore number, 1-7")
parser.add_argument("-o", default=None, type=str, help="image output destination")
args   = parser.parse_args()
basename=string.split(args.f, "/")[-1]

colors = ("#7fc97f", "#beaed4", "#fdc086")
logging.basicConfig(level=logging.DEBUG)


# START: SAFESAVEPLOT
if not args.o is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt
def safesaveplot(out=None, suffix=None, transparent=False, clf=True):
    fig = plt.gcf()
    fig.set_size_inches(8,6, forward=True)
    plt.tight_layout()
    if out:
        plt.savefig(out+suffix, transparent=transparent)
        if clf:
            plt.clf()
    else:
        plt.show()
# END: SAFESAVEPLOT


timeseries = np.loadtxt(args.f)
t = np.arange(timeseries.shape[0]) + 1

if args.T:
    plt.plot(t, np.sum(timeseries, axis=1))
    plt.xlabel("Time, ns")
    plt.ylabel("Energy Gap shift, cm$^{-1}$")
    plt.title("timeseries of {}, t=100ns@100ps sampling".format(basename))
    safesaveplot(args.o+"_T", '.png')

if args.H:
    nbins=20
    x, y, yfit = jh_hist_fit(np.sum(timeseries, axis=1), bins=nbins)
    plt.plot(x,    y, '-k', label="empirical")
    plt.plot(x, yfit, '--k', alpha=.5, linewidth=2, label="gaussian fit")
    plt.ylabel("Counts, {} bins".format(nbins))
    plt.xlabel("Energy Gap shift (bin centers), cm$^{-1}$")
    plt.title("histogram of {}, t=100ns@100ps sampling".format(basename))
    safesaveplot(args.o+"_H", '.png')

resid = np.array([rrid2resid_4BCL(i) for i in xrange(357)])
rmsf = np.std(timeseries, axis=0)
# Remove solvent with a mask
mask = (resid != 0)
timeseries = timeseries[:, mask]
resid = resid[mask]
rmsf  =  rmsf[mask]

mean_i_cm1 = np.sum(timeseries, axis=0)
cov_ii_cm2 = np.cov( (timeseries-mean_i_cm1).T )
import numpy.linalg as LA
lv, pca_iv = LA.eigh(cov_ii_cm2)
lv     =     lv[  ::-1]
pca_iv = pca_iv[:,::-1]

plt.semilogy(np.sqrt(lv), "o")
plt.ylim([1E-2, 1E2])
plt.xlabel("Principal component vector")
plt.ylabel("Root Covariance Spectrum, cm$^{-1}$")
safesaveplot(args.o if args.o is None else args.o+"_spectrum", ".png")

rmsf_mask = (rmsf != 0)
loadv = np.sum( np.square(pca_iv[rmsf_mask, :]) * lv[np.newaxis, :] / rmsf[rmsf_mask, np.newaxis], axis=0 )
plt.semilogy(loadv, "ro", label="Loading")
plt.semilogy(np.sqrt(lv), "o", label="Covariance")
plt.ylim([1E-2, 1E2])
plt.xlabel("Principal component vector")
plt.ylabel("Root Covariance & Load Spectrum, cm$^{-1}$")
plt.legend()
safesaveplot(args.o if args.o is None else args.o+"_load", ".png")

for pca_num, pca_vec in enumerate( (pca_iv[:,0], pca_iv[:,1], pca_iv[:,2]) ):
    pc_timeseries = np.dot(pca_vec, timeseries.T)
    plt.plot(pc_timeseries, label="PC{}".format(pca_num+1), linewidth=2, color=colors[pca_num], alpha=.8)
plt.legend()
safesaveplot(args.o if args.o is None else args.o+"pcts", '.png')



plt.bar(resid-.5, rmsf, width=1.0, color="#2F4F4F", linewidth=.05)

for pca_num, pca_vec in enumerate( (pca_iv[:,0], pca_iv[:,1], pca_iv[:,2]) ):
    # Mask out all zero-valued sites
    zero_mask = (rmsf != 0)
    resid   = resid[zero_mask]
    pca_vec = pca_vec[zero_mask]
    rmsf    = rmsf[zero_mask]

    width  = 1. / 3.
    offset = -.5 + width * (pca_num)
    neg_mask = (pca_vec < 0)
    pca_vec = np.abs(pca_vec)
    # Use each eigenvector's loading (the fraction of total variance contributed)
    #   as the height of the bar
    pca_load_vec = ((np.square(pca_vec) * lv[pca_num]) / np.square(rmsf)) 
    ind_max = np.argmax(pca_load_vec)
    logging.debug("Maximum load for PC{} = {}@resid{}".format(pca_num+1, pca_load_vec[ind_max], resid[ind_max]))
    logging.debug("EigVal: {}, rmsf: {}".format(lv[pca_num], rmsf[ind_max]))
    logging.debug("Sum load for PC{} = {} (should be > 1? average sum load over all PCs is 1.)".format(pca_num+1, np.sum(pca_load_vec)))
    logging.debug("Total magnitued of PC{} = {} (must be 1!!)".format(pca_num+1, np.sum(np.square(pca_vec))))
    pca_plotter = pca_load_vec * rmsf
    
    ALPHA=1
    plt.bar(resid[neg_mask]+offset, pca_plotter[neg_mask], width=width, linewidth=0, alpha=ALPHA, color=colors[pca_num])
    pos_mask = np.logical_not(neg_mask)
    plt.bar(resid[pos_mask]+offset, pca_plotter[pos_mask], width=width, linewidth=0, alpha=ALPHA, color=colors[pca_num])

    pc_cutoff  = np.amax(pca_vec) / 3.0
    pc_cutmask = pca_vec >= pc_cutoff
    pc_cutlist = np.where(pc_cutmask)[0]
    print("c{}p{} {}".format(args.c, pca_num+1, " ".join([str(x) for x in pc_cutlist])))
    for index_in_pc in pc_cutlist:
        pos = resid[index_in_pc]
        value = pca_plotter[index_in_pc]
        plt.text(pos, value, "i{}//PC{}".format(pos, pca_num+1), fontsize=9, 
                horizontalalignment='center', verticalalignment='bottom',
                color=colors[pca_num])


plt.xlim([-7,366])
plt.xlabel("Residue index")
plt.ylabel("Gap shift MSF, cm$^{-1}$")
plt.title("Gap RMSF of {}, t=100ns@100ps sampling".format(basename))
safesaveplot(args.o, '.svg', clf=False)

safesaveplot(args.o, '.png', transparent=True)


