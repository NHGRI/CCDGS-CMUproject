#!/bin/bash

# SOFTWARE
#GATK=/mnt/hanchardlab/bin/gatk-4.0.4.0/gatk

GATK4=/home/bin/gatk-4.0.11.0/gatk

GATK3=/home/bin/GATK/GenomeAnalysisTK.jar

# FILES

REF=/home/lesediw/GATK_Resources_bundle_v28_b37/human_g1k_v37.fasta
DBSNP=/home/lesediw/GATK_Resources_bundle_v28_b37/dbsnp_138.b37.vcf

RAWCOMBINEDVCF=/home/lesediw/gVCFs/SampleRenamedGvcfs/ALL/ALL_Round2WES.raw.snps.indels.g.vcf.gz

JOINTCALLRAW=/home/lesediw/JointCall/JointCall_ALL_Round2WES.raw.snps.indels.g.vcf.gz


# DIRECTORIES  

WORKDIR=/home/lesediw/JointCall
temp=/home/lesediw/JointCall/temp


PREFIX=RAW_JointCalled

cd $WORKDIR

##____________________________________________________________________________
#echo -e "Extraction variant types $(date)"
#
#$GATK4 --java-options "-Xmx900G" SelectVariants \
#-R $REF \
#-V $JOINTCALLRAW \
#--select-type-to-include SNP \
#-O $PREFIX.SNPs.vcf.gz \
#--tmp-dir=$temp
#
#
#$GATK4 --java-options "-Xmx900G" SelectVariants \
#-R $REF \
#-V $JOINTCALLRAW \
#--select-type-to-exclude SNP \
#-O $PREFIX.nonSNPs.vcf.gz \
#--tmp-dir=$temp
#
#echo -e "Extraction variant types $(date)"

#______________________________________________________________________________
# CAfGEN filter

#$GATK4 --java-options "-Xmx900G" VariantFiltration \
#-R $REF \
#-V $PREFIX.SNPs.vcf.gz \
#--filter-name "LowQual" --filter-expression "QUAL < 20 || QD < 2.0" \
#--filter-name "LowDP" --filter-expression "DP < 10" \
#--filter-name "StandBias" --filter-expression "FS > 60.0" \
#--filter-name "LowMQ" --filter-expression "MQ < 40.0 || MQRankSum < -12.5" \
#--filter-name "ReadPos" --filter-expression "ReadPosRankSum < -8.0" \
#--genotype-filter-name "LowGQ" --genotype-filter-expression "GQ < 20" \
#--genotype-filter-name "LowDP" --genotype-filter-expression "DP < 10 " \
#-O $PREFIX.SNPs.Filtered.vcf.gz
#
#
#$GATK4 --java-options "-Xmx900G" VariantFiltration \
#-R $REF \
#-V $PREFIX.nonSNPs.vcf.gz \
#--filter-name "LowQual" --filter-expression "QUAL < 20 || QD < 2.0" \
#--filter-name "LowDP" --filter-expression "DP < 10" \
#--filter-name "StandBias" --filter-expression "FS > 60.0" \
#--filter-name "ReadPos" --filter-expression "ReadPosRankSum < -8.0" \
#--genotype-filter-name "LowGQ" --genotype-filter-expression "GQ < 20" \
#--genotype-filter-name "LowDP" --genotype-filter-expression "DP < 10 " \
#-O $PREFIX.nonSNPs.Filtered.vcf.gz


#--genotype-filter-name "LowAD" --genotype-filter-expression "AD < 3" \
#--genotype-filter-name "LowAD" --genotype-filter-expression "AD < 3" \
##_______________________________________________________________________________
# Merge




java -Djava.io.tmpdir=./temp -Xmx900G -jar $GATK3 \
-T CombineVariants \
-R $REF \
-V $PREFIX.SNPs.Filtered.vcf.gz \
-V $PREFIX.nonSNPs.Filtered.vcf.gz \
-o $PREFIX.Filtered.vcf.gz \
-genotypeMergeOptions UNSORTED

#___________________________________________________________________________
