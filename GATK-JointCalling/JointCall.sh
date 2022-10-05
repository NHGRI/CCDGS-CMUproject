#!/bin/bash

# SOFTWARE
#GATK=/mnt/hanchardlab/bin/gatk-4.0.4.0/gatk

GATK=/home/bin/gatk-4.0.11.0/gatk

# FILES

REF=/home/lesediw/GATK_Resources_bundle_v28_b37/human_g1k_v37.fasta
DBSNP=/home/lesediw/GATK_Resources_bundle_v28_b37/dbsnp_138.b37.vcf

RAWCOMBINEDVCF=/home/lesediw/gVCFs/SampleRenamedGvcfs/ALL/ALL_Round2WES.raw.snps.indels.g.vcf.gz


# DIRECTORIES  

WORKDIR=/home/lesediw/JointCall
temp=/home/lesediw/JointCall/temp

cd $WORKDIR

echo -e "START Joint call $(date)"

$GATK --java-options "-Xmx900G" GenotypeGVCFs \
-R $REF \
-D $DBSNP \
-V $RAWCOMBINEDVCF \
-O JointCall_ALL_Round2WES.raw.snps.indels.g.vcf.gz \
--tmp-dir=$temp


echo -e "END Joint call $(date)"
