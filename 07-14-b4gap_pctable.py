from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import string
import numpy.random as RNG
SRCDIR=os.path.dirname(os.path.realpath(__file__))
SUBGROUP="/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-07-15"


parser = argparse.ArgumentParser(description="")
parser.add_argument("-f", type=str, help="file containing a list of principal components")
parser.add_argument("-wi", type=float, default=20)
parser.add_argument("-ht", type=float, default=6)
args   = parser.parse_args()

pc_elements = []
pc_labels  = []
with open(args.f) as f:
    for l in f:
        l_arr = l.split()
        print(l_arr[0]+": "+",".join(l_arr[1:]))
        pc_labels.append(l_arr[0])
        pc_elem = [int(x) for x in l_arr[1:]]
        pc_elements.append(pc_elem)

print(pc_elements)

flat_pc_elements = [item for sublist in pc_elements for item in sublist]
unique_elem = np.unique(flat_pc_elements)
has_unique_elem = []

for elem in unique_elem:
    has_this_unique_elem = []
    for pc_index, these_elems in enumerate(pc_elements):
        if elem in these_elems:
            has_this_unique_elem.append(pc_labels[pc_index])
    has_unique_elem.append(has_this_unique_elem)

for i in xrange(len(unique_elem)):
    print("{}: {}".format(unique_elem[i], has_unique_elem[i]))

scatter_pc = []
scatter_y = []
scatter_label = []
scatter_color = []
for elem, havers in zip(unique_elem, has_unique_elem):
    for haver in havers:
        scatter_pc.append(haver)
        scatter_y.append(elem)
        
scatter_x = [int(pc.split('c')[1].split('p')[0]) for pc in scatter_pc]
scatter_p = [int(pc.split('p')[1]) for pc in scatter_pc]
colors = ("#7fc97f", "#beaed4", "#fdc086")
scatter_c = [colors[p-1] for p in scatter_p]

print(scatter_x)
print(scatter_y)

jitter = .05
scatter_x = [x + RNG.normal(0, jitter) for x in scatter_x]
scatter_y = [y + RNG.normal(0, jitter) for y in scatter_y]

scatter_x = np.array(scatter_x)
scatter_y = np.array(scatter_y)
scatter_c = np.array(scatter_c)
for i in (2, 1, 0):
    mask = (scatter_c == colors[i])
    plt.scatter(scatter_y[mask], scatter_x[mask], color=colors[i], s=125, alpha=.8, label="PC{}".format(i+1))

plt.legend()

fig = plt.gcf()
fig.set_size_inches(args.wi, args.ht, forward=True)

plt.hlines(np.arange(7)+1, -7, 366, linewidths=2, linestyles='dashed', alpha=.1)


ax = plt.gca()
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.set_ylabel("chromophore number")
plt.tight_layout()
plt.savefig(SUBGROUP+"/pctable_overlay.png", transparent=True)


ax.invert_yaxis()
plt.tight_layout()
plt.savefig(SUBGROUP+"/pctable_overlay_inv.png", transparent=True)


ax.invert_yaxis()
plt.title("Significant Principal components ( > max / 3 )")
plt.xlabel("index")
plt.tight_layout()
plt.savefig(SUBGROUP+"/pctable.png")

plt.clf()


