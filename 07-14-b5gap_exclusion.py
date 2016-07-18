#!/home/jhaberstroh/anaconda/bin/python
from __future__ import division
from __future__ import print_function
import numpy as np
import argparse
import string
from jh_hist_fit import jh_hist_fit
from rrid2resid_4BCL import rrid2resid_4BCL
import os
from argparse import Namespace
SRCDIR=os.path.dirname(os.path.realpath(__file__))

# file=$SRCDIR/output/cdc_c${f}_2016-07-14.txt
# python 07-14-b4gapRMSF.py -f $file -c $f -o $SUBGROUP_BASE/2016-07-15/c${f}_gapRMSF
FIRST=True
for chromo in (np.arange(7) + 1):
    this_file = SRCDIR+"/output/cdc_c{}_2016-07-14.txt".format(chromo)
    output = os.environ["SUBGROUP_BASE"]+"/2016-07-15/exclusion-c{}".format(chromo)
    
    #! Argparse commented out because it was used before chromo was looped
    #! The "args" Namespace object was created to maintain a compatible interface
    # parser = argparse.ArgumentParser(description="")
    # parser.add_argument("-f", required=True, type=str, help="csv file with spaces for delimiters for energy gap")
    # parser.add_argument("-T", action="store_true", help="display timeseries")
    # parser.add_argument("-H", action="store_true", help="display histograms")
    # parser.add_argument("-c", type=int, default=0, help="Chromophore number, 1-7")
    # parser.add_argument("-o", default=None, type=str, help="image output destination")
    # args   = parser.parse_args()

    #! Production args
    args = Namespace(f=this_file, T=True, H=True, c=chromo, o=output)
    #! Debug args
    #! args = Namespace(f=this_file, T=True, H=True, c=chromo, o=None)

    basename=string.split(args.f, "/")[-1]
    
    if FIRST:
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
        FIRST=False
    
    timeseries = np.loadtxt(args.f)
    t = (np.arange(timeseries.shape[0]) + 1) / 10

    resid = np.array([rrid2resid_4BCL(i) for i in xrange(357)])
    solvent_mask = (resid != 0)
    excl_mask = np.logical_not((resid >= 120) * (resid <= 140))
    excl_timeseries = timeseries[:, excl_mask]
    
    if args.T:
        plt.plot(t, np.sum(timeseries, axis=1), 'k', alpha=.7, label="Original")
        plt.plot(t, np.sum(excl_timeseries, axis=1), 'r', alpha=.7, label="Exclusion")
        plt.xlabel("Time, ns")
        plt.ylabel("Energy Gap shift, cm$^{-1}$")
        plt.title("timeseries of {}, t=100ns@100ps sampling".format(basename))
        safesaveplot(args.o if args.o is None else args.o+"_T", '.png')
    
    if args.H:
        nbins=20
        x, y, yfit = jh_hist_fit(np.sum(timeseries, axis=1), bins=nbins)
        excl_x, excl_y, excl_yfit = jh_hist_fit(np.sum(excl_timeseries, axis=1), bins=nbins)
        plt.plot(x,    y, '-k', label="empirical")
        plt.plot(x, yfit, '--k', alpha=.5, linewidth=2, label="gaussian fit")
        plt.plot(excl_x, excl_y,     '-r', label="empirical-exclusion")
        plt.plot(excl_x, excl_yfit, '--r', alpha=.5, linewidth=2, label="gaussian fit-exclusion")
        plt.ylabel("Counts, {} bins".format(nbins))
        plt.xlabel("Energy Gap shift (bin centers), cm$^{-1}$")
        plt.legend(loc=2, prop={'size':8})
        plt.title("histogram of {}, t=100ns@100ps sampling".format(basename))
        safesaveplot(args.o if args.o is None else args.o+"_H", '.png')
    
