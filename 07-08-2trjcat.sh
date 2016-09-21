#!/bin/bash

set -o errexit
set -o nounset

MYDIR="/mnt/pup1/jhaberstroh/data/2016-06-fmo500ns/2016-06-fmo500ns"
JOB=830203
FMO_CONF=$HOME/Jobs/2015-07-FMOconf

cd $MYDIR/output

# for i in {1..3}; do 
#     files=""
#     for j in `seq $((i-1))01 $((i))00`; do
#         files="$files $MYDIR/md500-${JOB}-R${j}.xtc"
#     done
#     trjcat -f $files -o ${MYDIR}/output/md500-${JOB}-MERGE${i}00.xtc
# 
# 
#     TEMP_MDP="$(mktemp -p $MYDIR/output).mdp"
# (
# cat <<'MDPCONTENTS'
# # MDP generated by shell script
# integrator = md
# MDPCONTENTS
# ) > $TEMP_MDP
# 
# # Make the __CENTER.gro configuration from 4BCL.gro
#     trjconv -f $FMO_CONF/4BCL.gro -s $FMO_CONF/nvt/nvt.tpr \
#             -center -pbc mol -ur compact -o $MYDIR/output/__CENTER.gro << SYSTEM
# 19
# 0
# SYSTEM
# 
#     grompp -f $TEMP_MDP -c $MYDIR/output/__CENTER.gro -o $MYDIR/output/__CENTER.tpr \
#             -p $FMO_CONF/4BCL_pp.top
#     rm $TEMP_MDP
# 
#      trjconv -f md500-${JOB}-MERGE${i}00.xtc -o md500-${JOB}-PBC${i}00.xtc \
#              -s $MYDIR/output/__CENTER.tpr -dt 100 -center -ur compact -pbc mol << SYSTEM
# 1
# 0
# SYSTEM
# 
#      trjconv -f md500-${JOB}-PBC${i}00.xtc   -o md500-${JOB}-FIT${i}00.xtc \
#              -s $MYDIR/output/__CENTER.tpr -dt 100 -fit rot+trans << SYSTEM
# 19
# 0
# SYSTEM
# 
#     trjconv -f md500-${JOB}-CENTER${i}00.xtc -o md500-${JOB}-SEP-${i}-.gro \
#             -s ~/data/2015-09-FMO_conf/em/em.tpr -sep -dt 100 << SYSTEM 
# 0
# SYSTEM
# done

trjcat -f $MYDIR/md500-${JOB}-R{1..525}.xtc -o $MYDIR/output/job525_md500.xtc -dt 100

# Make the TEMP_MDP file!
TEMP_MDP="$(mktemp -p $MYDIR/output).mdp"
(
cat <<'MDPCONTENTS'
# MDP generated by shell script
integrator = md
MDPCONTENTS
) > $TEMP_MDP

# Make the __CENTER.gro configuration from 4BCL.gro
trjconv -f $FMO_CONF/4BCL.gro -s $FMO_CONF/nvt/nvt.tpr \
        -center -pbc mol -ur compact -o $MYDIR/output/__CENTER.gro << SYSTEM
19
0
SYSTEM

# Make the __CENTER.tpr configuration from 4BCL.gro
grompp -f $TEMP_MDP -c $MYDIR/output/__CENTER.gro -o $MYDIR/output/__CENTER.tpr \
       -p $FMO_CONF/4BCL_pp.top

# Center the protein for nice periodic boundaries
trjconv -f job525_md500.xtc     -o job525_md500-PBC.xtc \
        -s $MYDIR/output/__CENTER.tpr -dt 100 -center -ur compact -pbc mol << SYSTEM
1
0
SYSTEM

# Fit & rotate with the BCL molecules at the center
trjconv -f job525_md500-PBC.xtc -o job525_md500-FIT.xtc \
        -s $MYDIR/output/__CENTER.tpr -dt 100 -fit rot+trans << SYSTEM
19
0
SYSTEM

# Split trajectory into atoms
trjconv -f job525_md500-FIT.xtc -o job525_md500-SEP-.gro \
        -s $MYDIR/output/__CENTER.tpr -dt 100 -sep << SYSTEM 
0
SYSTEM

rm $TEMP_MDP
