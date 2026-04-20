#!/bin/bash

# SOFTWARE
#GATK=/mnt/hanchardlab/bin/gatk-4.0.4.0/gatk

GATK4=/home/bin/gatk-4.0.11.0/gatk
GATK3=/home/bin/GATK/GenomeAnalysisTK.jar

# FILES

REF=/home/lesediw/GATK_Resources_bundle_v28_b37/human_g1k_v37.fasta
DBSNP=/home/lesediw/GATK_Resources_bundle_v28_b37/dbsnp_138.b37.vcf
OMNI=/home/lesediw/GATK_Resources_bundle_v28_b37/1000G_omni2.5.b37.vcf
THOUSANDG=/home/lesediw/GATK_Resources_bundle_v28_b37/1000G_phase1.snps.high_confidence.b37.vcf
HAPMAP=/home/lesediw/GATK_Resources_bundle_v28_b37/hapmap_3.3.b37.vcf
MILLS=/home/lesediw/GATK_Resources_bundle_v28_b37/Mills_and_1000G_gold_standard.indels.b37.vcf

RAWCOMBINEDVCF=/home/lesediw/gVCFs/SampleRenamedGvcfs/ALL/ALL_Round2WES.raw.snps.indels.g.vcf.gz

JOINTCALLRAW=/home/lesediw/JointCall/JointCall_ALL_Round2WES.raw.snps.indels.g.vcf.gz


FilteredRecalibratedFile=/home/lesediw/JointCall/RAW_JointCalled.FilteredRecalibrated.vcf.gz

# DIRECTORIES  

WORKDIR=/home/lesediw/JointCall
temp=/home/lesediw/JointCall/temp

#___________________________________________________________________________

cd $WORKDIR

bcftools view -Oz --apply-filters .,PASS $FilteredRecalibratedFile > JointCalled_2_PASSonlyVariants.vcf.gz


