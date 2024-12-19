#!/bin/bash

set -x
set -euo pipefail

if [ $# -ne 1 ]
  then
    echo "usage: $0 <chr>"
    exit -1
fi


CHR=$1

BCFTOOLS="docker run -u $(id -u):$(id -g) --rm -v /data:/data staphb/bcftools"
BIN="docker run -u $(id -u):$(id -g) --rm -v /data:/data glimpse2-dev"

PANEL_NAME=CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz
PANEL_SITE=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased

PREFIX=/data/Runs

PANEL_PREFIX=$PREFIX/reference_panel
PANEL=$PANEL_PREFIX/$PANEL_NAME
PREP_PANEL=$PANEL_PREFIX/1000GP.${CHR}.bcf
SITES_PANEL=$PANEL_PREFIX/1000GP.${CHR}.sites.vcf.gz
SPLIT_REF=$PANEL_PREFIX/split/1000GP

MAPS_SRC=~/GLIMPSE/maps/genetic_maps.b38
MAP=$PREFIX/maps/genetic_maps.b38/${CHR}.b38.gmap.gz
CHUNKS=$PREFIX/chunks/${CHR}.txt

NPROC=$(nproc)

mkdir -p chunks $PANEL_PREFIX/split impute ligate maps/genetic_maps.b38

if [ ! -f $MAP ]; then
  echo "Copy map ... $MAPS_SRC -> $MAP"
  cp $MAPS_SRC/$(basename -- "$MAP") $MAP
fi

# if [ ! -f $PANEL ]; then
 echo "Download panel ... $PANEL"
 wget --directory-prefix=$PANEL_PREFIX -c $PANEL_SITE/${PANEL_NAME}{,.tbi}
# fi

if [ ! -f $PREP_PANEL ]; then
 echo "Prepare reference panel ... $PANEL -> $PREP_PANEL"

  $BCFTOOLS bash -c \
      "bcftools norm -m -any $PANEL -Ou --threads $NPROC | \
      bcftools view -m 2 -M 2 -v snps --threads $NPROC -Ob -o $PREP_PANEL"

 $BCFTOOLS bcftools index -f $PREP_PANEL --threads $NPROC
fi

if [ ! -f $SITES_PANEL ]; then
 echo "Extract sites from the reference panel ... $PREP_PANEL -> $SITES_PANEL"

 $BCFTOOLS bcftools view -G -Oz --threads $NPROC -o $SITES_PANEL $PREP_PANEL
 $BCFTOOLS bcftools index -f $SITES_PANEL --threads $NPROC 
fi

if [ ! -f $CHUNKS ]; then
 echo "Split the genome into chunks ... $SITES_PANEL -> $CHUNKS"
 $BIN GLIMPSE2_chunk \
    --input $SITES_PANEL \
    --region $CHR \
    --output $CHUNKS \
    --map $MAP \
    --sequential
fi

if true; then
 echo "Create binary reference model ... $CHUNKS -> $SPLIT_REF"
 while IFS="" read -r LINE || [ -n "$LINE" ];
 do
  printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
  IRG=$(echo $LINE | cut -d" " -f3)
  ORG=$(echo $LINE | cut -d" " -f4)

  $BIN GLIMPSE2_split_reference \
  	--reference $PREP_PANEL \
  	--map $MAP \
  	--input-region ${IRG} \
  	--output-region ${ORG} \
  	--output $SPLIT_REF
 done <$CHUNKS
fi
