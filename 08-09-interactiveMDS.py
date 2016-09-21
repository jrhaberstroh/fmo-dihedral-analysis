%pylab
import MDAnalysis.coordinates
import os
trr=os.environ["SCRATCH"] + "/2016-06-fmo500ns/dihedrals.trr"

def unravel_dihedrals(x_t):
    ## Unravel dihedral list from (x,y,z) to (sin,cos)
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
    return(x_t[:, :N_dihedral, :])


## Load the TRR data of dihedrals
TRR = MDAnalysis.coordinates.reader(trr)
x_t = np.array(TRR)
x_t = unravel_dihedrals(x_t)



## Compute covariance and eigenanalysis
dpca_t = x_t.reshape((1042,-1))
## Remove residues 120-140 that are in the loop
start=308*2
stop =394*2
dpca_t = np.hstack((dpca_t[:,:start], dpca_t[:,stop:] ))
## Remove resid 8-14 that unfold
stop =16*2
dpca_t = dpca_t[:,stop:]
dpcaCTR_t = dpca_t - np.mean(dpca_t, axis=0)
dpcaCovar = np.cov(dpcaCTR_t, rowvar=False)
from numpy.linalg import eigh
w,v = eigh(dpcaCovar)
w = w[::-1]
v = v[:, ::-1]
dpcaEIG_t = np.dot(dpcaCTR_t, v)

cumulative = zeros(dpcaEIG_t.shape[1])
for npc in xrange(7):
    vi = np.square(v[:,0])
    ind = (np.array(range(dpcaEIG_t.shape[1])) * .5) + .5
    plt.bar(ind, vi, width=.5, bottom=cumulative)
    cumulative = cumulative+vi


## Plot a simple 2d scatter of the first two PC's
my_cmap = plt.cm.get_cmap("winter")
c = my_cmap(np.linspace(0,1,dpcaEIG_t.shape[0]))
scatter = plt.scatter(dpcaEIG_t[:,0], dpcaEIG_t[:,1], color=c)

## Plot the cool 3d scatter of the first three PC's
from __future__ import print_function
"""
A very simple 'animation' of a 3D plot
"""
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import time
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
my_cmap = plt.cm.get_cmap("winter")
c = my_cmap(np.linspace(0,1,dpcaEIG_t.shape[0]))
scatter = ax.scatter(dpcaEIG_t[:,0], dpcaEIG_t[:,1], dpcaEIG_t[:,2], color=c)


## Do the MDS projection of principal components 1-6
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

data_ti = dpcaEIG_t[:,:6]
similarities = euclidean_distances(data_ti)
seed = np.random.RandomState(seed=3)
mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                   dissimilarity="precomputed", n_jobs=1)
pos = mds.fit(similarities).embedding_
my_cmap = plt.cm.get_cmap("winter")
c = my_cmap(np.linspace(0,1,data_ti.shape[0]))
plt.scatter(pos[:,0], pos[:,1], color=c)

