#!/bin/bash

# this is the script used for ID converting based on the chr:pos information
# works on biowulf
module load bedops

# input file format is *bed
input=./rsID_converting_example.bed

#reference file
dbsnp_file=/data/Hanserv/Reference/dbSNP/hg19.dbSNP151.bed

# cmd
bedmap --echo --echo-map-id --delim '\t' $input_file $dbsnp_file > output.rsid.txt
