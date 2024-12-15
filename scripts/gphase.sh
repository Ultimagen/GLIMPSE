#!/bin/bash

if [ $# -ne 2 ]
  then
    echo "usage: $0 <chr> <input>"
    exit -1
fi


CHR=$1
INPUT_BAM=$2
BIN="sudo docker run --rm -v .:/data -v ~/data/ug_bam/input/n1:/input glimpse2-dev /bin/"
DATA=/data
INPUT=/input
SPLIT_REF=reference_panel/split/1000GP
CHUNKS=chunks/${CHR}.txt


while IFS="" read -r LINE || [ -n "$LINE" ]; 
 do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    CHR=$(echo ${LINE} | cut -d" " -f2)
    REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
    REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
    OUT=impute/imputed

    echo $BIN/GLIMPSE2_phase \
     --bam-file $INPUT/${INPUT_BAM} \
     --reference $DATA/${SPLIT_REF}_${CHR}_${REGS}_${REGE}.bin \
     --output $DATA/${OUT}_${CHR}_${REGS}_${REGE}.bcf 
 done < ${CHUNKS}

