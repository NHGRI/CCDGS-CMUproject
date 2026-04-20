#!/bin/bash

GATK=/mnt/hanchardlab/bin/gatk-4.0.4.0/gatk

#GATK=/mnt/hanchardlab/bin/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar

REF=/mnt/hanchardlab/lesediw/GATK_Resources_bundle_v28_b37/human_g1k_v37.fasta

DBSNP=/mnt/hanchardlab/lesediw/GATK_Resources_bundle_v28_b37/dbsnp_138.b37.vcf

  

INPUTDIR=/mnt/hanchardlab/CAFGEN/wes_batch1_rnd2_mapped/analysis_ready_bams
OUTPUTDIR=/mnt/hanchardlab/lesediw/JointCallingCafgen/

MEM=20G 


cd $INPUTDIR

for bamfile in BWR0302*.bam  ### bam file name started with BWR0302
do
  SAMPLEID=`echo $bamfile | cut -d '.' -f 1`
  echo -e "START gVCF for $SAMPLEID $(date)"
  $GATK --java-options "-Xmx20G" HaplotypeCaller --TMP_DIR=./tmp \
  -R $REF \
  -I $bamfile \
  -O $OUTPUTDIR/$SAMPLEID-2.raw.snps.indels.g.vcf \
  -ERC GVCF \
  --dbsnp $DBSNP \
  -A DepthPerAlleleBySample \
  -A GenotypeSummaries \
  -A ChromosomeCounts \
  --QUIET
  echo -e "END gVCF for $SAMPLEID $(date)"
done