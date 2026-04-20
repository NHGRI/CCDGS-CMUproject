#!/usr/bin/env bash

gatk SelectVariants 	-O snp.vcf \
			-V Test_annotate.vcf \
			--select-type-to-include SNP
