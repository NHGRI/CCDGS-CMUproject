#!/bin/bash

# Software Path
plink=/mnt/hanchardlab/bin/plink_beta/plink

# Directories
WORKDIR=/mnt/hanchardlab/lesediw/CAfGEN_GWAS/QualityControl2

# Files 

RAWDATA=/mnt/hanchardlab/lesediw/CAfGEN_GWAS/RAW_DATA/RAWCafgenGWAS
EXCLUDE=/mnt/hanchardlab/lesediw/CAfGEN_GWAS/SupportFiles/Samples2Exclude.txt

PREFIX=CafgenGWAS

cd $WORKDIR

# Basic Stats

#$plink \
#--bfile $RAWDATA \
#--freq case-control \
#--missing \
#--het \
#--out $PREFIX-STATS



$plink \
--bfile $RAWDATA \
--chr 1-23 \
--mind 0.1 \
--maf 0.001 \
--geno 0.05 \
--hwe 0.000001 \
--make-bed \
--out $PREFIX-SnpClean4


##______________________________________________________________________________
##Sample QC

# Sample heterozygosity

$plink \
--bfile $PREFIX-SnpClean \
--split-x hg19 no-fail \
--make-bed \
--out $PREFIX-Split-X

$plink \
--bfile $PREFIX-Split-X \
--het \
--out $PREFIX

# LD pruning of SNPs

$plink \
--bfile $PREFIX-Split-X \
--indep-pairwise 100 10 0.2 \
--out $PREFIX-LD-100-10-02

$plink \
--bfile $PREFIX-Split-X \
--extract $PREFIX-LD-100-10-02.prune.in \
--make-bed \
--out $PREFIX-LD-Pruned

# Pairwise comparisions of samples, relatenes calculated

$plink \
--bfile $PREFIX-LD-Pruned \
--genome \
--threads 10 \
--out $PREFIX-genome

# Sex check

$plink \
--bfile $PREFIX-LD-Pruned \
--check-sex \
--out $PREFIX-SexCheck

$plink \
--bfile $PREFIX-Split-X \
--missing \
--out $PREFIX-SNPcleanMissing

#Sample Clean

$plink \
--bfile $PREFIX-Split-X \
--remove $EXCLUDE \
--make-bed \
--out $PREFIX-SampleClean

#convert to vcf

$plink \
--bfile $PREFIX-SampleClean \
--recode vcf-iid \
--out $PREFIX-CleanData


