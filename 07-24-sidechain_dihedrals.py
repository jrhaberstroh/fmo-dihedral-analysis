#!/home/jhaberstroh/anaconda/bin/python
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import argparse
import logging
# Uncomment for debugging messages
# logging.basicConfig(level=logging.DEBUG)
import os
SRCDIR=os.path.dirname(os.path.realpath(__file__))

GRO   ="/home/jhaberstroh/data/2015-09-FMO_conf/4BCL.gro"
AALIST=SRCDIR+"/input/aminoacid_dihedrals.txt"


def build_dihedrals(atom_list):
    all_dihedrals = []
    if len(atom_list) >= 4:
        for dih_count in xrange(len(atom_list) - 3):
            all_dihedrals.append(atom_list[dih_count:dih_count+4])
    return all_dihedrals


aa_dihedral = {}

with open(AALIST) as f:
    for l in f:
        l = l.split('#')[0]
        l = l.strip()
        if l != "":
            aa = l[:3]
            dihedrals = l[3:].split()
            aa_dihedral[aa] = dihedrals

amino_acids = aa_dihedral.keys()
for key in amino_acids:
    logging.debug(key + ": " + " ".join(aa_dihedral[key]))


IGNORE_TYPES=["BCL", "NA", "SOL"]

dihedral_list = []
dihedral_resid= []

with open(GRO) as f:
    _ = f.readline()
    n = f.readline()
    n = int(n)
    this_resid = 0
    this_dihedral_atomtype = []
    this_dihedral_atomnum  = []
    
    # Iterate through the number of lines specified by line 2 of .gro file
    for i in xrange(n):
        l = f.readline()
        label   = l[0:8]
        resid   = int(label[:5].strip())
        resname = label[-3:].strip()
        atomname= l[8:15].strip()
        atomnum = int(l[15:20].strip())

        if resname in IGNORE_TYPES:
            continue
        new_resid = (this_resid != resid)
        if new_resid:
            this_resid = resid
            logging.debug(resname.rstrip())
            badkey = (not resname in amino_acids)
            if badkey:
                raise ValueError("Amino acid {} found in GRO but not in aalist".format(resname))
            this_dihedral_atomtype = build_dihedrals(aa_dihedral[resname])
            this_dihedral_atomnum  = [ [-1, -1, -1, -1] for dihedral in this_dihedral_atomtype ]
            logging.debug(this_dihedral_atomtype)
        # this_dihedral_atomtype is set up, populate this_dihedral_atomnum
        for i, dihedral in enumerate(this_dihedral_atomtype):
            if atomname in dihedral:
                atom_index = dihedral.index(atomname)
                this_dihedral_atomnum[i][atom_index] = atomnum

        atomnum_flat = [item for sublist in this_dihedral_atomnum for item in sublist]
        atomnum_done = [atomnum != -1 for atomnum in atomnum_flat]
        if all(atomnum_done):
            for dihedral_atomnum in this_dihedral_atomnum:
                dihedral_list.append(dihedral_atomnum)
                dihedral_resid.append(this_resid)

for resid, dihedral in zip(dihedral_resid, dihedral_list):
    dihedral = ["{:4d}".format(atomnum) for atomnum in dihedral]
    print("{}#".format(resid) + " ".join(dihedral))
