#! /bin/bash

# directories and paths
WORKDIR=/data/Hanserv/NMAMP
cd $WORKDIR

#### mode and modules
sinteractive --mem=32g --cpus-per-task=8 --time 24:00:00

#### 1) alignment
mkdir bwa_output

# for all of the fastq files alignment, make a swarm file for higher speed
BWAINDEX=/fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa
FASTQDIR=./NMAMP-277946675/FASTQ_Generation_2021-07-07_23_07_44Z-436755319

for i in $(ls $FASTQDIR); do
FILE1=$(ls $FASTQDIR/$i | head -n 1)
FILE2=$(ls $FASTQDIR/$i | tail -n 1)
NAME=$(echo $FILE1 | cut -d _ -f 1-3)
SAMPLEID=$(echo $FILE1 | cut -d _ -f 1)
header=$(zcat $FASTQDIR/$i/$FILE1 | head -n 1)
id=$(echo $header | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
sm=$(echo $header | cut -f 1-2 -d":" | sed 's/@//' | sed 's/:/_/g')
pu=$(echo $header | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/./g')
idx=$(echo $header | grep -Eo "[ATGCN]+$")
echo "bwa mem -M -O 30 -E 4 -T 20 -v 3 -t 36 -R $(echo \"@RG\\tID:$id\\tSM:$SAMPLEID\\tPU:$pu"."$idx\\tLB:$sm".lib1"\\tPL:ILLUMINA\") \
$BWAINDEX $FASTQDIR/$i/$FILE1 $FASTQDIR/$i/$FILE2 \
| samblaster -M | samtools fixmate - - | samtools sort -t 36 -O bam -o ./bwa_output/${NAME}_sorted.bam"
done >> ./scripts/01_bwa.swarm

swarm --file ./scripts/01_bwa.swarm --module bwa,samblaster,samtools --gb-per-process 32 -t 8 --verbose 1 --time 8:00:00 

#### 2) merge the 4 replicates bam files by samples
mkdir bam

# 2.1 get the unique sample list
for i in $(ls $RAWBAMFILES); do 
SAMPLEID=$(echo $i | cut -d _ -f 1)
echo $SAMPLEID
done >> smapleList.txt

cat smapleList.txt | uniq > sample_list.txt
rm smapleList.txt

# 2.2 get the merge swarm file
for s in $(cat sample_list.txt); do 
echo "samtools merge ./bam/${s}_merged.bam `find $RAWBAMFILES -name ${s}_*_sorted.bam | xargs`"
done >> ./scripts/02_merge.swarm

swarm --file ./scripts/02_merge.swarm --module samtools --gb-per-process 32 -t 8 --verbose 1 --time 2:00:00 

#### 3) sort and index the merged bam files
BAMFILES=./bam

for b in $(ls $BAMFILES); do
SAMPLEID=$(echo $b | cut -d . -f 1)
samtools sort $BAMFILES/$b -t 36 -O bam -o $BAMFILES/${SAMPLEID}_sorted.bam
samtools index $BAMFILES/${SAMPLEID}_sorted.bam
done
