#!/bin/bash

# Make DIR and copy over
mkdir -p $1/
cp EUR_1KG_phase3_samples.tsv $1/
cd $1/

# Download 1KG data
wget -m -np -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/

# Download plink
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
unzip *.zip

if [ $2 = "EUR" ];
then
	for i in {1..22}
	do
		./plink --recode 12 transpose --vcf-half-call missing --vcf ALL.chr$i\_GRCh38.genotypes.20170504.vcf.gz -keep EUR_1KG_phase3_samples.tsv --out EUR.1KGphase3.GRCh38.chr$i
	done
else
	
	for i in {1..22}
        do
                ./plink --recode 12 transpose --vcf-half-call missing --vcf ALL.chr$i\_GRCh38.genotypes.20170504.vcf.gz --out ALL.1KGphase3.GRCh38.chr$i
        done
fi

gzip *.tped
