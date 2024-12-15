#!/bin/bash

if [ $# -ne 1 ]
  then
    echo "usage: $0 <input>"
    exit -1
fi


INPUT_BAM=$1

cat chromosomes | parallel -j 1 ./gphase.sh {1} $1 | parallel
