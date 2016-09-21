#!/bin/bash
set -o nounset
set -o errexit
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

TRAJ=/home/jhaberstroh/data/2016-06-fmo500ns/2016-06-fmo500ns/output/md_170nsX10ps.xtc
OUTDIR=/home/jhaberstroh/data/2016-06-fmo500ns/2016-08-cdc_rerun

echo "Running ANGLE"
echo "[ SidechainDihedrals ]"                         > $SRCDIR/output/dihedral_list.ndx
$SRCDIR/07-24-sidechain_dihedrals.py | cut -d'#' -f2 >> $SRCDIR/output/dihedral_list.ndx
$SRCDIR/07-24-sidechain_dihedrals.py | cut -d'#' -f1  > $SRCDIR/output/dihedral_list_resid.ndx
g_angle -f $TRAJ -n $SRCDIR/output/dihedral_list.ndx -type dihedral -all \
        -or $OUTDIR/dihedrals_170nsX10ps.trr                      \
        -od $OUTDIR/dihedrals_dist_170nsX10ps
