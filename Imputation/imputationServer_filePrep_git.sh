#! /bin/bash

# mode and modules
sinteractive --mem=32g 
module load plink
module load bcftools

# 1. directories and paths
IMPDIR=/data/Hanserv/NMGTH3A/Imputation
name=RAW_DATA 

cd $IMPDIR

# 2. data for imputation
# 2.1 copy the input file / data for imputation
cp ../PLINK_output_v4/$name.bed $IMPDIR/data
cp ../PLINK_output_v4/$name.bim $IMPDIR/data
cp ../PLINK_output_v4/$name.fam $IMPDIR/data

# 2.2 prepare the supporting files
mkdir supporting_files
cd supporting_files

wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

unzip HRC-1000G-check-bim-v4.2.7.zip
gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

# 2.3 Create a frequency file
cd $IMPDIR/data 
plink --freq --bfile RAW_DATA --out RAW_DATA

# 2.4 Script to check plink .bim files against HRC/1000G for strand, id names, positions, alleles, ref/alt assignment
perl ../supporting_files/HRC-1000G-check-bim.pl \
-b RAW_DATA.bim -f RAW_DATA.frq -r \
../supporting_files/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

# 2.5 clean up and split by chr
sh Run-plink.sh

# 2.6 Create vcf using vcfCooker
# git clone vcfCooker to /data/Hanserv/tools: git clone https://github.com/statgen/gotcloud.git
# 2.6.1 prepare vcfCooker
----------------------------------------
cd /data/$USER
git clone https://gcc02.safelinks.protection.outlook.com/?url=https%3A%2F%2Fgithub.com%2Fstatgen%2Fgotcloud.git&amp;data=05%7C01%7Cyixing.han%40nih.gov%7C81182309543b498b2f8e08da2c4733f3%7C14b77578977342d58507251ca2dc2b06%7C0%7C0%7C637870981560300251%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C3000%7C%7C%7C&amp;sdata=wBEVFmXbPYFJmOOnQWTwzueBSsC11vWHHaQHkboT7%2FI%3D&amp;reserved=0
cd gotcloud/src
module load cmake
make
cd bin
./vcfCooker --help
----------------------------------------
cd /data/Hanserv/tools/
git clone https://github.com/statgen/gotcloud.git
cd gotcloud/src
module load cmake
make
cd bin
pwd

# tool path
/data/Hanserv/tools/gotcloud/src/bin/vcfCooker

cd $IMPDIR/data 

# 2.6.2 copy genome.fa file from /fdb/igenome so it can be indexed
cd /data/Hanserv/Reference
mkdir hg19
cd hg19
cp /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa .
REFHG19=/data/Hanserv/Reference/hg19/genome.fa
 
# 2.6.3
# test on one chr
/data/Hanserv/tools/gotcloud/src/bin/vcfCooker --in-bfile RAW_DATA-updated-chr1 --ref $REFHG19 --out RAW_DATA-updated-chr1_test --write-vcf

# for all of the chr
for i in $(seq 1 23); do 
echo "/data/Hanserv/tools/gotcloud/src/bin/vcfCooker \
--in-bfile $IMPDIR/data/${name}-updated-chr${i} \
--ref $REFHG19 \
--out $IMPDIR/data/${name}-updated-chr${i}_clean \
--write-vcf" 
done >> ./scripts/05_vcfCooker.swarm

swarm --file ./scripts/05_vcfCooker.swarm -t 2 -g 4 --time=4:00:00

# 2.7 remove "chr" prefix as the Imputation Server set with chr for hg38 and without chr for hg19
# test for one file
(grep ^"#" RAW_DATA-updated-chr22_clean.vcf; grep -v ^"#" RAW_DATA-updated-chr22_clean.vcf | sed 's:^chr::ig' | sort -k1,1n -k2,2n) | bgzip -c > RAW_DATA-updated-chr22_clean_test.vcf.gz

for i in $(seq 1 23); do
echo "(grep ^\"\#\" $IMPDIR/data/${name}-updated-chr${i}_clean.vcf; \
grep -v ^\"\#\" $IMPDIR/data/${name}-updated-chr${i}_clean.vcf | sed 's:^chr::ig' | sort -k1,1n -k2,2n) | \
bgzip -c > $IMPDIR/data/${name}-updated-chr${i}_clean.vcf.gz"
done >> ./scripts/06_sortzipVCF.swarm

swarm --file ./scripts/06_sortzipVCF.swarm -t 2 -g 4 --time=4:00:00

# 2.8 checkVCF
# test on one chr
/usr/bin/python /data/Hanserv/tools/checkVCF/checkVCF.py -r /data/Hanserv/tools/checkVCF/hs37d5.fa \
-o /data/Hanserv/NMGTH3A/Imputation/data/RAW_DATA-updated-chr1_test  /data/Hanserv/NMGTH3A/Imputation/data/RAW_DATA-updated-chr1_test.vcf.gz

# for all files
for i in $(seq 1 23); do 
echo "/usr/bin/python /data/Hanserv/tools/checkVCF/checkVCF.py \
-r /data/Hanserv/tools/checkVCF/hs37d5.fa \
-o $IMPDIR/data/${name}-updated-chr${i}_clean $IMPDIR/data/${name}-updated-chr${i}_clean.vcf.gz" 
done >> ./scripts/07_checkVCF.swarm

swarm --file ./scripts/07_checkVCF.swarm -t 2 -g 4 --time=4:00:00

# 2.9 zip vcf, ref: https://www.biostars.org/p/59492/
for i in $(seq 1 23); do 
echo "bgzip -c $IMPDIR/data/${name}-updated-chr${i}_clean.vcf > $IMPDIR/data/${name}-updated-chr${i}_clean.vcf.gz
tabix -p vcf $IMPDIR/data/${name}-updated-chr${i}_clean.vcf.gz" 
done >> ./scripts/09_zipVCF.swarm

swarm --file ./scripts/09_zipVCF.swarm -t 2 -g 4 --time=4:00:00

# 3. submit to the imputation server and download after receive the email. Unzip the files and move to QC step.
# https://imputation.biodatacatalyst.nhlbi.nih.gov/#!
# https://imputationserver.sph.umich.edu/
