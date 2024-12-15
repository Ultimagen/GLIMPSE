#!/bin/bash

if [ $# -ne 1 ]
  then
    echo "usage: $0 <chr>"
    exit -1
fi


CHR=$1
BIN="sudo docker run --rm -v .:/data glimpse2-dev /bin/"
DATA=/data
PANEL_NAME=CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz
PANEL_SITE=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased
PANEL=reference_panel/$PANEL_NAME
PREP_PANEL=reference_panel/1000GP.${CHR}.bcf
SITES_PANEL=reference_panel/1000GP.${CHR}.sites.vcf.gz
SPLIT_REF=reference_panel/split/1000GP
MAP=maps/genetic_maps.b38/${CHR}.b38.gmap.gz
CHUNKS=chunks/${CHR}.txt
MAPS_SRC=~/src/GLIMPSE/maps/genetic_maps.b38
MAPS_DST=maps/genetic_maps.b38

NPROC=4

mkdir -p chunks reference_panel/split impute ligate maps/genetic_maps.b38

if [ ! -f $MAP ]; then
  echo "Copy map ... $MAPS_SRC -> $MAP"
  cp $MAPS_SRC/$(basename -- "$MAP") $MAP
fi

if [ ! -f $PANEL ]; then
 echo "Download panel ... $PANEL"
 wget --directory-prefix=reference_panel -c $PANEL_SITE/${PANEL_NAME}{,.tbi}
fi

if [ ! -f $PREP_PANEL ]; then
 echo "Prepare reference panel ... $PANEL -> $PREP_PANEL"
 bcftools norm -m -any $PANEL -Ou --threads $NPROC | \
	bcftools view -m 2 -M 2 -v snps --threads $NPROC -Ob -o $PREP_PANEL
 bcftools index -f $PREP_PANEL --threads $NPROC
fi

if [ ! -f $SITES_PANEL ]; then
 echo "Extract sites from the reference panel ... $PREP_PANEL -> $SITES_PANEL"
 bcftools view -G -Oz -o $SITES_PANEL $PREP_PANEL
 bcftools index -f $SITES_PANEL
fi

if [ ! -f $CHUNKS ]; then
 echo "Split the genome into chunks ... $SITES_PANEL -> $CHUNKS"
 $BIN/GLIMPSE2_chunk \
	--input $DATA/$SITES_PANEL \
	--region $CHR \
	--output $DATA/$CHUNKS \
	--map /data/$MAP \
	--sequential
fi

if true; then
 echo "Create binary reference model ... $CHUNKS -> $DATA/$SPLIT_REF"
 while IFS="" read -r LINE || [ -n "$LINE" ];
 do
  printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
  IRG=$(echo $LINE | cut -d" " -f3)
  ORG=$(echo $LINE | cut -d" " -f4)

  $BIN/GLIMPSE2_split_reference \
  	--reference $DATA/$PREP_PANEL \
  	--map $DATA/$MAP \
  	--input-region ${IRG} \
  	--output-region ${ORG} \
  	--output $DATA/$SPLIT_REF
 done <$CHUNKS
fi
