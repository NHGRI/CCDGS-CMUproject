#!/bin/bash

# Software Path
plink=/mnt/hanchardlab/bin/plink_beta/plink

# Directories
WORKDIR=/mnt/hanchardlab/lesediw/CAfGEN_GWAS/RAW_DATA

# Files 
RAWFILE=/mnt/hanchardlab/CAFGEN/AGBL_MERGE_2017.11_2018.07_A3_Oct2_2018/PLINK_021018_1113/AGBL_MERGE_2017.11_2018.07_A3_Oct2_2018

$plink \
--file $RAWFILE \
--make-bed \
--out RAW_DATA

RAWDATA=$WORKDIR/RAW_DATA

KEEP=/mnt/hanchardlab/lesediw/CAfGEN_GWAS/SupportFiles/CafgenKeepSamples.txt  # list of samples you want to keep
UpdateIDs=/mnt/hanchardlab/lesediw/CAfGEN_GWAS/SupportFiles/UpdateIDs.txt # if you want to change sample IDs se <https://www.cog-genomics.org/plink/1.9/data#update_indiv>
UpdatePhenoSex=/mnt/hanchardlab/lesediw/CAfGEN_GWAS/SupportFiles/UpdatePhenoSex.txt # updating phenotype and sex of samples  <https://www.cog-genomics.org/plink/1.9/input#pheno>

cd $WORKDIR

$plink \
--bfile $RAWDATA \
--keep $KEEP \
--make-bed \
--out RAWCafgenSamplesOnly

$plink \
--bfile RAWCafgenSamplesOnly \
--update-ids $UpdateIDs \
--make-bed \
--out RAWCafgen-Updated-IDs

$plink \
--bfile RAWCafgen-Updated-IDs \
--update-sex $UpdatePhenoSex \
--pheno $UpdatePhenoSex --pheno-name PHENO \
--make-bed \
--out RAWCafgen-Updated-IDs-Pheno-Sex

