#!/bin/bash
set -o nounset
set -o errexit
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

$SRCDIR/07-14-a4plotrmsf.py
## $SRCDIR/07-14-a4plotgyrate.py
$SRCDIR/07-14-a4plotrotamers.py
##t Temporary commented for debugging
$SRCDIR/07-14-b4gapRMSF.sh
##t $SRCDIR/07-14-b5gap_exclusion.py

for f in {1..7}; do
    GAP=$SUBGROUP_BASE/2016-07-15/c${f}_gapRMSF.png
    BACK=$SUBGROUP_BASE/2016-07-15/backbone_rmsf.png
    OUT=$SUBGROUP_BASE/2016-07-15/superimpose_c${f}.png
    composite -compose atop $GAP $BACK $OUT
done

OUTDIR=$SUBGROUP_BASE/2016-07-15

PC=$SUBGROUP_BASE/2016-07-15/pctable_overlay.png
PCI=$SUBGROUP_BASE/2016-07-15/pctable_overlay_inv.png
BACK=$SUBGROUP_BASE/2016-07-15/backbone_rmsf_wide.png
composite $PC  $BACK -compose atop -geometry '99.25%x93%+38+15' $OUTDIR/superimpose_bbpc.png
composite $PCI $BACK -compose atop -geometry '99.25%x93%+38+15' $OUTDIR/superimpose_bbpci.png
    
ROTAMER=$SUBGROUP_BASE/2016-07-15/rotamersPPO_4BCL_wide.png
composite $PC  $ROTAMER -compose atop -geometry '99.25%x93%+38+15' $OUTDIR/superimpose_rotpc.png
composite $PCI $ROTAMER -compose atop -geometry '99.25%x93%+38+15' $OUTDIR/superimpose_rotpci.png


