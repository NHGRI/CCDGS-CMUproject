#!/usr/bin/env bash

gatk SplitVcfs 	-I /mnt/lustre/groups/CBBI1243/EDMOND/Data/edmond.vcf.gz \
		--INDEL_OUTPUT indel.vcf.gz \
		--SNP_OUTPUT snp.vcf.gz \
   		--STRICT false
