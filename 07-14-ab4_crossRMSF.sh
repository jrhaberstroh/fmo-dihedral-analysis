#!/bin/bash
set -o nounset
set -o errexit
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

for f in {1..7}; do
    GAP=$SUBGROUP_BASE/2016-07-15/c${f}_gapRMSF.png
    BACK=$SUBGROUP_BASE/2016-07-15/backbone_rmsf.png
    OUT=$SUBGROUP_BASE/2016-07-15/superimpose_c${f}.png
    composite -compose atop -blend 20 $GAP $BACK $OUT
done
