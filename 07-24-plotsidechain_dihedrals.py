from __future__ import division, print_function
import numpy as np
import argparse
import rrid2resid_4BCL
from MDAnalysis.coordinates.xdrfile import libxdrfile2 as XDR
import logging
logging.basicConfig(level=logging.DEBUG)
import os
SRCDIR=os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser(description="")
parser.add_argument("-f", type=str, required=True, help="[input] .trr file of sidechains")
parser.add_argument("-resid", type=str, required=True, help="[input] list of residue ids for .trr file")
parser.add_argument("-cdc", type=str, default=None, help="[input] cdc trajectory input file")
parser.add_argument("-onlyplot", type=int, nargs="*", default=None, help="[input] List of dihedrals to plot (default=None, i.e. plot all dihedrals)")
parser.add_argument("-sc", action="store_true", help="[input] plot sin/cos?")
parser.add_argument("-chromo", type=str, default=None, help="[input] Chromophore number label for title")
parser.add_argument("-o", type=str, default=None, help="[output] image destination basename")
args   = parser.parse_args()
# START: SAFESAVEPLOT
if not args.o is None:
    import matplotlib
    matplotlib.use("agg")
import matplotlib.pyplot as plt
def safesaveplot(out=None, suffix=None, transparent=False, clf=True):
    fig = plt.gcf()
    if out:
        plt.savefig(out+suffix, transparent=transparent)
        if clf:
            plt.clf()
    else:
        plt.show()
# END: SAFESAVEPLOT

cdc = None
if not args.cdc is None:
    cdc = np.loadtxt(args.cdc)

resid = np.loadtxt(args.resid)
resid = resid.astype(int)


x_t = []

DIM=3
# read the number of "atoms", the number of dihedral angles
natoms = XDR.read_trr_natoms(args.f)
# allocate coordinate array of the right size and type
# (the type float32 is crucial to match the underlying C-code!!)
x = np.zeros((natoms, DIM), dtype=np.float32)
v = np.zeros((natoms, DIM), dtype=np.float32)
f = np.zeros((natoms, DIM), dtype=np.float32)
# allocate unit cell box
box = np.zeros((DIM, DIM), dtype=np.float32)
# open file
TRR = XDR.xdrfile_open(args.f, 'r')
status = XDR.exdrOK
while status == XDR.exdrOK:
    status,step,time, _, _, _, _ = XDR.read_trr(TRR, box, x, v, f)
    x_alloc = np.copy(x)
    x_t.append(x_alloc)
XDR.xdrfile_close(TRR)

x_t = np.array(x_t)

# Compute the number of dihedrals
padding = 0
even_shape=(x_t.shape[1] % 2 == 0)
if even_shape:
    ## Check if there are two zeros padded at the end.
    ## If there are, remove one from the bulk value.
    end_two = x_t[-1, -1, -2:]
    if all(end_two == 0.0):
        padding = -1
else:
    ## If the array shape[1] is odd, pad the array by one, and set the 
    ##   padding value to remove the last two dihedrals
    blank = np.zeros( (x_t.shape[0], 1, x_t.shape[2]) )
    x_t = np.concatenate(x_t, blank)
    padding = -2
N_dihedral_bulk = (x_t.shape[1]) // 2 * 3
N_dihedral = N_dihedral_bulk + padding
x_t = x_t.reshape( (x_t.shape[0], N_dihedral_bulk, 2) )
x_t = x_t[:, :N_dihedral, :]
logging.debug(x_t.shape)
logging.debug(resid.shape)

angle_t = np.arctan2(x_t[:, :, 0], x_t[:, :, 1]) * 180.0 / np.pi



if args.onlyplot is None:
    args.onlyplot = np.array(range(len(resid)))
else:
    logging.debug("Received onlyplot: {}".format(" ".join(map(str, args.onlyplot))))
for i in xrange(len(resid)):
    if i not in args.onlyplot:
        continue
    if args.sc:
        plt.plot(x_t[:, i, 1], label="sin")
        plt.plot(x_t[:, i, 0], label="cos")
        plt.xlabel("time/100ps")
        plt.ylabel("cos and sin of angle")
        plt.ylim([-1, 1])
        plt.title("timeseries for dihedral {}, resid {}".format(i, resid[i]))
        plt.legend()
        safesaveplot(args.o if args.o is None else args.o + "_sc_r{}i{}".format(resid[i], i), ".png")

    plt.plot(angle_t[:, i], 'ko', markersize=2)
    plt.xlabel("time/100ps")
    plt.ylabel("Angle, degrees")
    plt.ylim([-180.0, 180.0])
    if args.chromo is None:
        plt.title("timeseries for dihedral {}, resid {}".format(i, resid[i]))
    else:
        plt.title("timeseries for dihedral {}, r{}::c{}".format(i, resid[i], args.chromo))
    # Include cdc in the figure
    ax_left  = plt.gca()
    ax_right = ax_left.twinx()
    rrid = rrid2resid_4BCL.resid2rrid_4BCL(resid[i])
    ax_right.plot(cdc[:,rrid], "r.-", linewidth=.2, markersize=.5)
    ax_right.set_ylabel("cdc shift, cm$^{-1}$", color='r')
    mean_shift = np.mean(cdc[:,rrid])
    ax_right.set_ylim(mean_shift - 50, mean_shift + 50)
    for tl in ax_right.get_yticklabels():
        tl.set_color('r')
    safesaveplot(args.o if args.o is None else args.o + "_r{}i{}".format(resid[i], i), ".png")

