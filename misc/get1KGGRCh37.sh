#!/bin/bash

# Make DIR and copy over
mkdir -p $1/
cp EUR_1KG_phase3_samples.tsv $1/
cd $1/

# Download 1KG data
wget -r -nd -l1 --no-parent ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/

# Download plink
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
unzip *.zip

function pwait() {
    while [ $(jobs -p | wc -l) -ge $3 ]; do
        sleep 1
    done
}


if [ $2 = "EUR" ];
then
	for i in {1..22}
	do

		# Plink convert
		./plink --recode 12 transpose --vcf-half-call missing --vcf ALL.chr$i.phase3\_shapeit2\_mvncall\_integrated\_v5a.20130502.genotypes.vcf.gz -keep EUR_1KG_phase3_samples.tsv --out EUR.1KG.GRCh37.chr$i &
		pwait $3
	done
else
	
	for i in {1..22}
        do

		# Plink convert
                ./plink --recode 12 transpose --vcf-half-call missing --vcf ALL.chr$i.phase3\_shapeit2\_mvncall\_integrated\_v5a.20130502.genotypes.vcf.gz --out ALL.1KG.GRCh37.chr$i &
		pwait $3
        done
fi

wait

gzip *.tped
