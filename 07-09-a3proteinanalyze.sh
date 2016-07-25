#!/bin/bash
set -o errexit
set -o nounset
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

MYDIR="/mnt/pup1/jhaberstroh/data/2016-06-fmo500ns/2016-06-fmo500ns"
JOB=830203
FMO_CONF=$HOME/Jobs/2015-07-FMOconf
FIT_TPR="$MYDIR/output/__CENTER.tpr"
TRAJ="$MYDIR/output/ALL_md500-FIT.xtc"

if [ ! -e $SRCDIR/output ]; then
    mkdir $SRCDIR/output
else if [ ! -d $SRCDIR/output ]; then
    echo "ERROR: $SRCDIR/output exists, but is not a directory. Cannot continue."
    exit 888
fi; fi

RMSF=false
if [ "$RMSF" = "true" ]; then
    echo "Running RMSF"
    g_rmsf -f $TRAJ -s $FIT_TPR -o $SRCDIR/output/rmsf.txt
fi

ANGLE=true
if [ "$ANGLE" = "true" ]; then
    ##t echo "Running ANGLE"
    ##t echo "[ SidechainDihedrals ]"                         > $SRCDIR/output/dihedral_list.ndx
    ##t $SRCDIR/07-24-sidechain_dihedrals.py | cut -d'#' -f2 >> $SRCDIR/output/dihedral_list.ndx
    ##t $SRCDIR/07-24-sidechain_dihedrals.py | cut -d'#' -f1  > $SRCDIR/output/dihedral_list_resid.ndx
    ##t g_angle -f $TRAJ -n $SRCDIR/output/dihedral_list.ndx -type dihedral -all \
    ##t         -or $SCRATCH/2016-06-fmo500ns/dihedrals.trr                      \
    ##t         -od $SCRATCH/2016-06-fmo500ns/dihedral_dist
    python $SRCDIR/07-24-plotsidechain_dihedrals.py -f $SCRATCH/2016-06-fmo500ns/dihedrals.trr \
            -resid $SRCDIR/output/dihedral_list_resid.ndx -o $SUBGROUP_BASE/2016-07-29/dihedral
fi

RGYRO=false
if [ "$RGYRO" = "true" ]; then
    echo "Running RGYRO"
    g_gyrate -f $TRAJ -s $FIT_TPR -o $SRCDIR/output/gyrate_4BCL
fi

CHI=false
if [ "$CHI" = "true" ]; then
    echo "Running sidechain rotamers"
    g_chi -s $FIT_TPR -f $TRAJ -o $SRCDIR/output/rotamers_4BCL.xvg
fi
