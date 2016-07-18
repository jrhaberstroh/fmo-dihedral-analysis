#!/bin/bash
ORIGIN_SCRATCH="/scratch2/scratchdirs/jhabers"
LOCAL_DATA="$HOME/data/"
NAME="2016-06-fmo500ns"

rsync -r jhabers@edison.nersc.gov:${ORIGIN_SCRATCH}/${NAME}  \
        ${LOCAL_DATA}/${NAME}
