#!/usr/bin/env bash

gatk CombineGVCFs -R /mnt/lustre/groups/CBBI1243/Data/hg19/ucsc.hg19.fasta -O edmond.gvcf.gz \
     $(for i in *gvcf.gz; do echo -V $i; done | tr '\n' ' ' | sed 's/ $/\n/g')
