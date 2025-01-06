#!/bin/bash

set -x
set -euo pipefail

if [ $# -ne 2 ]
  then
    echo "usage: $0 <chr> <input>"
    exit -1
fi


CHR=$1
INPUT_BAM=$2
BIN="docker run -u $(id -u):$(id -g) --rm -v /data:/data glimpse2-dev"

DATA=/data

PREFIX=/data/Runs

PANEL_PREFIX=$PREFIX/reference_panel
SPLIT_REF=$PANEL_PREFIX/split/1000GP
CHUNKS=$PREFIX/chunks/${CHR}.txt
OUTPUT=$PREFIX/output

mkdir -p $OUTPUT


while IFS="" read -r LINE || [ -n "$LINE" ]; 
 do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    CHR=$(echo ${LINE} | cut -d" " -f2)
    REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
    REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)

    $BIN GLIMPSE2_phase \
     --bam-file ${INPUT_BAM} \
     --reference ${SPLIT_REF}_${CHR}_${REGS}_${REGE}.bin \
     --output ${OUTPUT}/${CHR}_${REGS}_${REGE}.bcf 
 done < ${CHUNKS}

