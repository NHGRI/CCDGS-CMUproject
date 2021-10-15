#!/usr/bin/env bash

perl /mnt/lustre/groups/1000genomes/annotation/annovar/table_annovar.pl \
     edmond_filtered.vcf.gz \
     /mnt/lustre/groups/1000genomes/annotation/annovar/humandb/ \
     -buildver hg19 \
     -out edmond \
     -remove \
     -protocol refGene,knownGene,ensGene,gnomad_exome,avsnp147,dbnsfp30a,exac03,esp6500siv2_all,dbscsnv11,intervar_20180118,ALL.sites.2014_10,clinvar_20180603,cytoBand \
     -operation g,g,g,f,f,f,f,f,f,f,f,f,r \
     -nastring '.' \
     -polish \
     -vcfinput
