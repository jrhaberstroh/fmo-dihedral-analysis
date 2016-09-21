#!/bin/bash
set -o nounset
set -o errexit

# Specify the run parameters
for chromo in 367 368 369 370 371 372 373; do
    CHROMO=$chromo
    FNAME=md_long
    
    # Specify the major folders
    SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
    TOPOLOGY=$HOME/data/2015-09-FMO_conf/
    TRAJECTORIES=$HOME/data/2015-09-FmoEd/structures/
    OUTDIR=/home/jhaberstroh/data/2016-06-fmo500ns/2016-08-cdc_rerun
    
    # Specify gromacs files to use
    MDRUN="mpirun -np 6 /home/jhaberstroh/local/gromacs_umb_mpi-4.6.7/bin/mdrun_umb_mpi"
    GROMPP=grompp
    
    # Specify the input files
    START=$TOPOLOGY/4BCL.gro
    MDP=$SRCDIR/output/CDC_${CHROMO}.mdp
    INDEX=$TOPOLOGY/index.ndx
    TOP=$TOPOLOGY/4BCL_pp.top
    
	cat <<- MDP > $MDP
	title           = CDC rerun for c$CHROMO
	; Run parameters
	integrator      = md            ; leap-frog integrator
	
	; Bond parameters
	continuation    = no            ; first dynamics run
	constraint_algorithm = lincs    ; holonomic constraints
	constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
	lincs_iter      = 1             ; accuracy of LINCS
	lincs_order     = 4             ; also related to accuracy
	
	; Cutoffs and electrostatics
	pbc             = xyz           ; 3-D PBC
	nstlist         = 10
	ns_type         = grid          ; search neighboring grid cells
	rlist           = 1.2           ; short-range neighborlist cutoff (in nm)
	rcoulomb        = 1.2           ; short-range electrostatic cutoff (in nm)
	rvdw            = 1.2           ; short-range van der Waals cutoff (in nm)
	DispCorr        = EnerPres      ; account for cut-off vdW scheme
	coulombtype     = PME           ; Particle Mesh Ewald for long-range electrostatics
	pme_order       = 4             ; cubic interpolation
	fourierspacing  = 0.24          ; grid spacing for FFT
	
	; Pulling for CDC
	pull                    = constant-force
	pull-nstxout            = 1
	pull-nstfout            = 1
	pull-geometry           = direction-periodic
	pull-dim                = Y Y Y
	pull-ngroups            = 1
	pull-group1 = BCL_BCL_$CHROMO
	pull-k1 = 0.0
	MDP
    
    $GROMPP -c $START \
            -n $INDEX -p $TOP \
            -f $MDP \
            -o $OUTDIR/$FNAME \
            -po $OUTDIR/$FNAME \
            -maxwarn 1
    $MDRUN -v -deffnm $OUTDIR/$FNAME -rerun ${TRAJECTORIES}/${FNAME}.xtc \
            >> $OUTDIR/cdc${CHROMO}_50nsX10ps.cdc
    
    
    rm $OUTDIR/*.edr
    rm $OUTDIR/*.log
    rm $OUTDIR/*.mdp
    rm $OUTDIR/*.tpr
    rm $OUTDIR/*.xvg
    rm $OUTDIR/\#*
    
done
