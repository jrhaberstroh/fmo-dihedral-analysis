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
#     echo "Running ANGLE"
#     echo "[ SidechainDihedrals ]"                         > $SRCDIR/output/dihedral_list.ndx
#     $SRCDIR/07-24-sidechain_dihedrals.py | cut -d'#' -f2 >> $SRCDIR/output/dihedral_list.ndx
#     $SRCDIR/07-24-sidechain_dihedrals.py | cut -d'#' -f1  > $SRCDIR/output/dihedral_list_resid.ndx
#     g_angle -f $TRAJ -n $SRCDIR/output/dihedral_list.ndx -type dihedral -all \
#             -or $SCRATCH/2016-06-fmo500ns/dihedrals.trr                      \
#             -od $SCRATCH/2016-06-fmo500ns/dihedral_dist
    python $SRCDIR/07-24-plotsidechain_dihedrals.py -f $SCRATCH/2016-06-fmo500ns/dihedrals.trr \
            -resid $SRCDIR/output/dihedral_list_resid.ndx -o $SUBGROUP_BASE/2016-08-05/dihedral/c367 \
           -cdc $SRCDIR/output/cdc_c1_2016-07-14.txt -onlyplot 284 287 452 456 \
            -chromo 367
    python $SRCDIR/07-24-plotsidechain_dihedrals.py -f $SCRATCH/2016-06-fmo500ns/dihedrals.trr \
            -resid $SRCDIR/output/dihedral_list_resid.ndx -o $SUBGROUP_BASE/2016-08-05/dihedral/c368 \
           -cdc $SRCDIR/output/cdc_c2_2016-07-14.txt -onlyplot 396 404 518 522 547 606 630 \
            -chromo 368
    python $SRCDIR/07-24-plotsidechain_dihedrals.py -f $SCRATCH/2016-06-fmo500ns/dihedrals.trr \
            -resid $SRCDIR/output/dihedral_list_resid.ndx -o $SUBGROUP_BASE/2016-08-05/dihedral/c369 \
           -cdc $SRCDIR/output/cdc_c3_2016-07-14.txt -onlyplot 10 14 745 784 787 796 984 986 \
            -chromo 369
    python $SRCDIR/07-24-plotsidechain_dihedrals.py -f $SCRATCH/2016-06-fmo500ns/dihedrals.trr \
            -resid $SRCDIR/output/dihedral_list_resid.ndx -o $SUBGROUP_BASE/2016-08-05/dihedral/c370 \
           -cdc $SRCDIR/output/cdc_c4_2016-07-14.txt -onlyplot 809 \
            -chromo 370
    python $SRCDIR/07-24-plotsidechain_dihedrals.py -f $SCRATCH/2016-06-fmo500ns/dihedrals.trr \
            -resid $SRCDIR/output/dihedral_list_resid.ndx -o $SUBGROUP_BASE/2016-08-05/dihedral/c371 \
           -cdc $SRCDIR/output/cdc_c5_2016-07-14.txt -onlyplot 1 10 14 63 809 813 828 834 838 990 \
            -chromo 371
    python $SRCDIR/07-24-plotsidechain_dihedrals.py -f $SCRATCH/2016-06-fmo500ns/dihedrals.trr \
            -resid $SRCDIR/output/dihedral_list_resid.ndx -o $SUBGROUP_BASE/2016-08-05/dihedral/c372 \
           -cdc $SRCDIR/output/cdc_c6_2016-07-14.txt -onlyplot 94 166 174 188 \
            -chromo 372
    python $SRCDIR/07-24-plotsidechain_dihedrals.py -f $SCRATCH/2016-06-fmo500ns/dihedrals.trr \
            -resid $SRCDIR/output/dihedral_list_resid.ndx -o $SUBGROUP_BASE/2016-08-05/dihedral/c373 \
           -cdc $SRCDIR/output/cdc_c7_2016-07-14.txt -onlyplot 247 254 643 663 \
            -chromo 373
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
