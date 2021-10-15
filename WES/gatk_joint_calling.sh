#!/usr/bin/env bash
gatk GenotypeGVCFs -R /mnt/lustre/groups/CBBI1243/Data/hg19/ucsc.hg19.fasta -O edmond.vcf.gz -V edmond.gvcf.gz
