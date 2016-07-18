from __future__ import division
from __future__ import print_function
import numpy as np
import argparse
import string
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


# START: SAFESAVEPLOT
if not args.o is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt
def safesaveplot(out=None, suffix=None):
    if out:
        plt.savefig(out+suffix, bbox_inches="tight")
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
mask = (resid != 0)
resid = resid[mask]
rmsf  =  rmsf[mask]


plt.bar(resid-.5, rmsf, width=1.0, linewidth=0)    
plt.xlim([-7,366])
plt.xlabel("Residue index")
plt.ylabel("Gap shift MSF, cm$^{-1}$")
plt.title("Gap RMSF of {}, t=100ns@100ps sampling".format(basename))
safesaveplot(args.o, '.png')

