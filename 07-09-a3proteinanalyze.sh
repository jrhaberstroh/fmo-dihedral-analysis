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

ANGLE=false
if [ "$ANGLE" = "true" ]; then
    echo "Running ANGLE"
    ##! NOTE: None of this works
    ##! USE CREATED INDEX
    ##! mk_angndx -s $FIT_TPR -n $SRCDIR/output/dihedrals.ndx
    ##! g_angle -f $TRAJ -n $SRCDIR/output/dihedrals.ndx -or $SRCDIR/output/angles_4BCL
    ##! USE MANUAL INDEX
    ##! g_angle -f $TRAJ -n $SRCDIR/07-12-4BCL_propers.ndx -or $SRCDIR/output/angles_4BCL
fi

RGYRO=false
if [ "$RGYRO" = "true" ]; then
    echo "Running RGYRO"
    g_gyrate -f $TRAJ -s $FIT_TPR -o $SRCDIR/output/gyrate_4BCL
fi

CHI=true
if [ "$CHI" = "true" ]; then
    echo "Running sidechain rotamers"
    g_chi -s $FIT_TPR -f $TRAJ -o $SRCDIR/output/rotamers_4BCL.xvg
fi
