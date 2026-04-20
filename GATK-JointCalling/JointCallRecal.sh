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


# DIRECTORIES  

WORKDIR=/home/lesediw/JointCall
temp=/home/lesediw/JointCall/temp


PREFIX=RAW_JointCalled

NT=12
MEM=900G

cd $WORKDIR


# BUILDING RECALIBRATION MODELS
#-------------------------------

# SNPs

#$GATK4 --java-options "-Xmx900G" VariantRecalibrator \
#-R $REF \
#-V $PREFIX.Filtered.vcf.gz \
#--resource hapmap,known=false,training=true,truth=true,prior=15.0:$HAPMAP \
#--resource omni,known=false,training=true,truth=false,prior=12.0:$OMNI \
#--resource 1000G,known=false,training=true,truth=false,prior=10.0:$THOUSANDG \
#--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$DBSNP \
#-an QD \
#-an MQ \
#-an FS \
#-an SOR \
#-an MQRankSum \
#-an ReadPosRankSum \
#-an InbreedingCoeff \
#-mode SNP \
#-O $PREFIX.recalibrate_SNP.recal \
#--tranches-file $PREFIX.recalibrate_SNP.tranches \
#--rscript-file $PREFIX.recalibrate_SNP_plots.R
#
## INDELs
#
#$GATK4 --java-options "-Xmx900G" VariantRecalibrator \
#-R $REF \
#-V $PREFIX.Filtered.vcf.gz \
#--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$DBSNP \
#--resource mills,known=true,training=true,truth=true,prior=12.0:$MILLS \
#-an QD \
#-an FS \
#-an SOR \
#-an MQRankSum \
#-an ReadPosRankSum \
#-an InbreedingCoeff \
#-mode INDEL \
#--max-gaussians 4 \
#-O $PREFIX.recalibrate_INDEL.recal \
#--tranches-file $PREFIX.recalibrate_INDEL.tranches \
#--rscript-file $PREFIX.recalibrate_INDEL_plots.R

# APPLYING RECALLIBRATION TO CALL SET
#-------------------------------------

# SNPs

$GATK4 --java-options "-Xmx900G" ApplyVQSR \
-R $REF \
-V $PREFIX.Filtered.vcf.gz \
-O $PREFIX.recalibrated_snps_raw_indels.vcf.gz \
--truth-sensitivity-filter-level 99.0 \
--tranches-file $PREFIX.recalibrate_SNP.tranches \
--recal-file $PREFIX.recalibrate_SNP.recal \
-mode SNP


# INDELS

$GATK4 --java-options "-Xmx900G" ApplyVQSR \
-R $REF \
-V $PREFIX.recalibrated_snps_raw_indels.vcf.gz \
-O $PREFIX.FilteredRecalibrated.vcf.gz \
--truth-sensitivity-filter-level 99.0 \
--tranches-file $PREFIX.recalibrate_INDEL.tranches \
--recal-file $PREFIX.recalibrate_INDEL.recal \
-mode INDEL



#___________________________________________________________________________
