#! /bin/bash

# directories and paths
WORKDIR=/data/Hanserv/NMAMP
cd $WORKDIR

#### mode and modules
sinteractive --mem=32g --cpus-per-task=8 --time 24:00:00

#### coverage for the target sequencing data
#### 1. file prepare
# 1.1 prepare the bed file
cat 184163_870_amplicon_targets_BCM.csv | awk -F "," '{print $4"\t"$5"\t"$6"\t"$3"\t"$1"\t"$2}' | tail -n +2 \
>> 184163_870_amplicon_targets_BCM.bed

# 1.2 Liftover the bed file from hg19 to hg38, Ref: https://hpc.nih.gov/apps/crossmap.html
module load crossmap || fail "could not load crossmap module"
if [[ ! -f hg19ToHg38.over.chain.gz ]]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz  
fi
crossmap bed hg19ToHg38.over.chain.gz 184163_870_amplicon_targets_BCM.bed \
| awk -F "\t" '{print $8"\t"$9"\t"$10"\t"$11}' \
>> 184163_870_amplicon_targets_BCM_hg38.bed # region bed file on hg38

BEDFILE=./184163_870_amplicon_targets_BCM_hg38.bed
BAMFILES=./bam

# 1.3 gene bed files
for gene in $(cat $BEDFILE | awk -F "\t" '{print $4}' | uniq); do
chr=$(cat $BEDFILE | grep $gene | awk -F "\t" '{print $1}' | uniq)
start=$(cat $BEDFILE | grep $gene | awk -F "\t" '{print $2}' | sort -n | head -1)
end=$(cat $BEDFILE | grep $gene | awk -F "\t" '{print $3}' | sort -n | tail -1)
echo "$chr,$start,$end,$gene" | tr "," "\t"
done >> 184163_870_amplicon_targets_35_gene_hg38.bed

# 1.4 unique gene list
cat 184163_870_amplicon_targets_BCM.bed | awk -F "\t" '{print $4}' | uniq \
>> 184163_870_amplicon_targets_unique_gene_list.txt
# 1.5 target region length
cat 184163_870_amplicon_targets_BCM_hg38.bed | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$3-$2"\t"$4}' \
>> 184163_870_amplicon_targets_length.txt
# 1.6 gene bed files with gene length
cat 184163_870_amplicon_targets_35_gene_hg38.bed | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$3-$2"\t"$4}' \
>> 184163_870_amplicon_targets_35_gene_length.txt


#### 2. get coverage
mkdir coverage

# 2.1 check the bam files depth/coverage on for each lane
# samtools depth *bamfile*  |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'
RAWBAMFILES=./bwa_output

for i in $(ls $RAWBAMFILES); do
cov=$(samtools depth $RAWBAMFILES/$i |  awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR; print sqrt(sumsq/NR - (sum/NR)**2)}')
echo $i $cov
done >> ./coverage/raw_bam_coverage_per_lane.txt

# 2.2 check the merged bam files depth/coverage for each sample on target region
BAMFILES=./bam
for i in $(ls $BAMFILES/*merged_sorted.bam); do
SAMPLEID=$(basename $i | cut -d _ -f 1)
echo "samtools depth -b $BEDFILE $i "
done >> ./coverage/merged_bam_coverage_per_sample.txt

#
samtools depth -b $BEDFILE -H $BAMFILES/*merged_sorted.bam \
| awk -F '\t' '(NR==1) {split($0,header);N=0.0;next;} {N++;for(i=3;i<=NF;i++) a[i]+=int($i);} END { for(x in a) print header[x], a[x]/N;}'

# coverage after markduplicates and base recalibration
for i in $(ls $BAMFILES/*merged_sorted_markdup_bqsr.bam); do
SAMPLEID=$(basename $i | cut -d _ -f 1)
cov=$(samtools depth -b $BEDFILE $i |  awk '{c++;s+=$3}END{print s/c}')
cov_per_base=$(samtools depth -ab $BEDFILE $i |  awk '{c++;s+=$3}END{print s/c}')
echo $SAMPLEID $cov $cov_per_base
done >> ./coverage/merged_sorted_markdup_bqsr_bam_region_coverage.txt

# 2.3 coverage by bedtools
# Ref:https://www.biostars.org/p/279140/
# To get the mean depth of coverage for each interval specified in your BED file, use:
bedtools coverage -a $BEDFILE -b ./bam/PA1002_merged_sorted_markdup_bqsr.bam -mean
# If you wanted to get the overall mean, then just pipe this into awk and get the average of the average, something like:
bedtools coverage -a $BEDFILE -b ./bam/PA1002_merged_sorted_markdup_bqsr.bam -mean | awk '{total+=$7} END {print total/NR}'
# output the per base read depth for each region in the BED file
bedtools coverage -a $BEDFILE -b ./bam/PA1002_merged_sorted_markdup_bqsr.bam -d 

# coverage by bedtools on region
for i in $(ls $BAMFILES/*merged_sorted_markdup_bqsr.bam); do
SAMPLEID=$(basename $i | cut -d _ -f 1)
echo "bedtools coverage -a $BEDFILE -b $i -mean > ./coverage/${SAMPLEID}_MeanCoverage_region.bedgraph"
echo "bedtools coverage -a $BEDFILE -b $i -hist > ./coverage/${SAMPLEID}_PerBaseDepth_region.bedgraph"
done >> ./scripts/07_bedtools_coverage_region.swarm 

swarm --file ./scripts/07_bedtools_coverage_region.swarm --module bedtools -t 4 -g 8 --gres=lscratch:400 --time=24:00:00

bedtools coverage -a ./184163_870_amplicon_targets_BCM_hg38.bed -b ./bam/PA1004_merged_sorted_markdup_bqsr.bam -mean > ./coverage/PA1004_MeanCoverage_region.bedgraph
bedtools coverage -a ./184163_870_amplicon_targets_BCM_hg38.bed -b ./bam/PA1004_merged_sorted_markdup_bqsr.bam -hist > ./coverage/PA1004_PerBaseDepth_region.bedgraph
bedtools coverage -a ./184163_870_amplicon_targets_BCM_hg38.bed -b ./bam/PA1081_merged_sorted_markdup_bqsr.bam -mean > ./coverage/PA1081_MeanCoverage_region.bedgraph
bedtools coverage -a ./184163_870_amplicon_targets_BCM_hg38.bed -b ./bam/PA1081_merged_sorted_markdup_bqsr.bam -hist > ./coverage/PA1081_PerBaseDepth_region.bedgraph

# PerBaseDepth output file column
# 1) depth
# 2) # bases at depth
# 3) size of A
# 4) % of A at depth

#######################################################################
## MeanCoverage
# Median coverage on the region for the cohort
for j in $(ls ./coverage/PA*MeanCoverage_region.bedgraph); do
SAMPLEID=$(basename $j | cut -d _ -f 1)
cat $j | awk '{print $0"\t""'$SAMPLEID'""\t"$4"-'$SAMPLEID'"}'
done >> ./coverage/all_samples_MeanCoverage_region.bedgraph

# sort by gene and then by coverage
cat ./coverage/all_samples_MeanCoverage_region.bedgraph | sort -k4,4 -k6,6 -k5,5n \
>> ./coverage/all_samples_MeanCoverage_region_sorted.bedgraph

# get the median and mean coverage for each gene each sample
ALLCOVFILE=./coverage/all_samples_MeanCoverage_region_sorted.bedgraph
echo "gene,id,gene_id,median,mean" | tr "," "\t" >> ./coverage/all_samples_MeanCoverage_region_median.txt
for gene_id in $(cat $ALLCOVFILE | awk -F "\t" '{print $7}' | uniq); do
cov=$(cat $ALLCOVFILE | grep $gene_id | awk -F "\t" '{print $5}')
#median=$(echo $cov | sort -n | awk ' NR>1 {a[++i]=$2} END {if(i%2==1)print a[(i+1)/2]; else print (a[i/2]+a[i/2+1])/2}' RS=" ")
median=$(echo $cov | awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' RS=" ")
value=$(echo "$cov" | awk -f ./scripts/calculate.awk )
#echo "${cov[@]/%/$'\n'}"
#median=$(echo $cov | sort -n | awk ' NR>1 {a[++i]=$2} END {if(i%2==1)print a[(i+1)/2]; else print (a[i/2]+a[i/2+1])/2}')
mean=$(echo $cov | awk '{s+=$1} END {print s/NR}' RS=" ")
gene=$(echo $gene_id | cut -d- -f1)
id=$(echo $gene_id | cut -d- -f2)
echo $gene,$id,$gene_id,$median,$mean | tr "," "\t"
done >> ./coverage/all_samples_MeanCoverage_region_median.txt 

# get the average median coverage for each gene
MEDIANFILE=./coverage/all_samples_MeanCoverage_region_median.txt
echo "gene,mean" | tr "," "\t" >> ./coverage/all_samples_MeanCoverage_region_median_mean.txt
for gene in $(cat $MEDIANFILE | awk -F "\t" '{print $1}' | uniq); do
cov=$(cat $MEDIANFILE | grep $gene | awk -F "\t" '{print $4}')
mean=$(echo $cov | awk '{s+=$1} END {print s/NR}' RS=" ")
echo $gene,$mean| tr "," "\t"
done >> ./coverage/all_samples_MeanCoverage_region_median_mean.txt 

#######################################################################
## PerBaseDepth
# Median coverage on the region for the cohort
PERBASEFILE=./coverage/PA*PerBaseDepth_region.bedgraph
for k in $(ls $PERBASEFILE); do
SAMPLEID=$(basename $k | cut -d _ -f 1)
cat $k | grep -v all | awk '{print $0"\t""'$SAMPLEID'""\t"$4"-'$SAMPLEID'"}'
done >> ./coverage/all_samples_PerBaseDepth_region.bedgraph

# sort by gene and then by coverage
cat ./coverage/all_samples_PerBaseDepth_region.bedgraph | sort -k4,4 -k9,9 -k5,5n \
>> ./coverage/all_samples_PerBaseDepth_region_sorted.bedgraph

# get the median and mean coverage for each gene each sample
ALLCOVFILE=./coverage/all_samples_PerBaseDepth_region_sorted.bedgraph
echo "gene,id,gene_id,median,mean" | tr "," "\t" >> ./coverage/all_samples_PerBaseDepth_region_median.txt
for gene_id in $(cat $ALLCOVFILE | awk -F "\t" '{print $10}' | uniq); do
cov=$(cat $ALLCOVFILE | grep $gene_id | awk -F "\t" '{print $5}')
#median=$(echo $cov | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
median=$(echo $cov | awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' RS=" ")
mean=$(echo $cov | awk '{s+=$1} END {print s/NR}' RS=" ")
gene=$(echo $gene_id | cut -d- -f1)
id=$(echo $gene_id | cut -d- -f2)
echo $gene,$id,$gene_id,$median,$mean | tr "," "\t"
done >> ./coverage/all_samples_PerBaseDepth_region_median.txt

# get the average median coverage for each gene
MEDIANFILE=./coverage/all_samples_PerBaseDepth_region_median.txt
echo "gene,mean" | tr "," "\t" >> ./coverage/all_samples_PerBaseDepth_region_median_mean.txt
for gene in $(cat $MEDIANFILE | awk -F "\t" '{print $1}' | uniq); do
cov=$(cat $MEDIANFILE | grep $gene | awk -F "\t" '{print $4}')
mean=$(echo $cov | awk '{s+=$1} END {print s/NR}' RS=" ")
echo $gene,$mean| tr "," "\t"
done >> ./coverage/all_samples_PerBaseDepth_region_median_mean.txt

#######################################################################


PERBASEFILE=./coverage/PA*PerBaseDepth_region.bedgraph
for k in $(ls $PERBASEFILE); do
SAMPLEID=$(basename $k | cut -d _ -f 1)
gcov=$(cat $k | grep "all" | awk -F "\t" '{print $2}')
frac=$(cat $k | grep "all" | awk -F "\t" '{print $5}')
cov_cumul=$(echo $frac | awk '{s=0; for(i=1;i<=NF;i++) $i=s+=$i; $i=s}1')
# gcov_cumul=$(echo $(( 1 - $cov_cumul)) | bc)
# echo 1 $cov_cumul | awk '{print $1 - $2}'
# echo "$a $b" | awk '{print $1 - $2}'
echo $gcov,$frac,$cov_cumul,$gcov_cumul | tr "," "\t" 
done >> ./coverage/${SAMPLEID}_all_gov_frac_covcumul.txt


#######################################################################
## Proportion of bases with <10X coverage in >95% of samples
# get "all" rows from all samples
PERBASEFILE=./coverage/PA*PerBaseDepth_region.bedgraph
for k in $(ls $PERBASEFILE); do
SAMPLEID=$(basename $k | cut -d _ -f 1)
cat $k | grep "all" | awk '{print $0"\t""'$SAMPLEID'"}'
done >> ./coverage/all_samples_PerBaseDepth_region_ALL.bedgraph

# sort
cat ./coverage/all_samples_PerBaseDepth_region_ALL.bedgraph | sort -k2,2n \
 >> ./coverage/all_samples_PerBaseDepth_region_ALL_sorted.bedgraph


##################
## proportion of bases with <10x coverage in 95% of samples per gene
# lower then 10x coverage reads and fraction
for c in $(cat ./coverage/all_samples_PerBaseDepth_region_sorted.bedgraph | awk -F "\t" '{print $6}'); do
if [ "$c" -lt 10 ]; then
echo $c
fi
done

# get the rows that cov < 10x
cat ./coverage/all_samples_PerBaseDepth_region_sorted.bedgraph | awk '$5<10' \
>> ./coverage/all_samples_PerBaseDepth_region_sorted_lt10x.txt

LOWCOVFILE=./coverage/all_samples_PerBaseDepth_region_sorted_lt10x.txt
for gene in $(cat $LOWCOVFILE | awk -F "\t" '{print $4}' | uniq); do
depth=$(cat $LOWCOVFILE | grep $gene | awk -F "\t" '{print $5}')
d_sum=$(echo $depth | awk '{s+=$1} END {printf "%.0f\n", s}' RS=" ")
length=$(cat $LOWCOVFILE | grep $gene | awk -F "\t" '{print $7}')
l_sum=$(echo $length | awk '{s+=$1} END {printf "%.0f\n", s}' RS=" ")
echo $gene,$d_sum,$l_sum | tr "," "\t" 
done >> ./coverage/all_samples_PerBaseDepth_region_sorted_lt10x_sum.txt
