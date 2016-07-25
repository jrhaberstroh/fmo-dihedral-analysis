from __future__ import division, print_function
import numpy as np
import argparse
from MDAnalysis.coordinates.xdrfile import libxdrfile2 as XDR
import logging
logging.basicConfig(level=logging.DEBUG)
import os
SRCDIR=os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser(description="")
parser.add_argument("-f", type=str, required=True, help="[input] .trr file of sidechains")
parser.add_argument("-resid", type=str, required=True, help="[input] list of residue ids for .trr file")
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
    x_t.append(x)
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
    x_t = np.pad(x_t, ((0,0), (0,0), (0,1)), "constant", constant_values=0)
    padding = -2
N_dihedral_bulk = (x_t.shape[1]) // 2 * 3
N_dihedral = N_dihedral_bulk + padding
x_t = x_t.reshape( (x_t.shape[0], N_dihedral_bulk, 2) )
x_t = x_t[:, :N_dihedral, :]
print(x_t.shape)
print(resid.shape)

angle_t = np.arctan2(x_t[:, :, 0], x_t[:, :, 1]) * 180.0 / np.pi

for i in xrange(len(resid)):
    ##? Is this useful? Let's just look at the angles for now...
    ##? plt.plot(x_t[:, i, 1], label="sin")
    ##? plt.plot(x_t[:, i, 0], label="cos")
    ##? plt.xlabel("time/100ps")
    ##? plt.ylabel("cos and sin of angle")
    ##? plt.ylim([-1, 1])
    ##? plt.title("timeseries for dihedral {}, resid {}".format(i, resid[i]))
    ##? plt.legend()
    ##? safesaveplot(args.o if args.o is None else args.o + "_sc_r{}i{}".format(resid[i], i), ".png")

    plt.plot(angle_t[:, i])
    plt.xlabel("time/100ps")
    plt.ylabel("Angle, degrees")
    plt.ylim([-180.0, 180.0])
    plt.title("timeseries for dihedral {}, resid {}".format(i, resid[i]))
    plt.legend()
    safesaveplot(args.o if args.o is None else args.o + "_angle_r{}i{}".format(resid[i], i), ".png")

