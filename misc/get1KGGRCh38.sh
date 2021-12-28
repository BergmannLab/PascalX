#!/bin/bash

# Make DIR and copy over
mkdir -p $1/
cp EUR_1KG_phase3_samples.tsv $1/
cp updateRSID.py $1/
cd $1/

# Download 1KG data
wget -r -nd -l1 --no-parent ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/

# Download snpdb
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz

# Download plink
if [ $4 = "tped" ];
then
	wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
	unzip *.zip
fi

function pwait() {
    while [ $(jobs -p | wc -l) -ge $1 ]; do
        sleep 1
    done
}

for i in {1..22}
do
	# Add rsids
        python3 updateRSID.py $i CCDG\_14151\_B01\_GRM\_WGS\_2020-08-05\_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz All_20180418.vcf.gz &
	pwait $3
done

wait

if [ $4 = "tped" ];
then
	if [ $2 = "EUR" ];
	then
		for i in {1..22}
		do
			./plink --recode 12 transpose --vcf-half-call missing --vcf CCDG\_14151\_B01\_GRM\_WGS\_2020-08-05\_chr$i.filtered.shapeit2-duohmm-phased.snpid.vcf.gz -keep EUR_1KG_phase3_samples.tsv --out EUR.1KG.GRCh38.chr$i &
			pwait $3
		done
	else
	
		for i in {1..22}
        	do
                	./plink --recode 12 transpose --vcf-half-call missing --vcf CCDG\_14151\_B01\_GRM\_WGS\_2020-08-05\_chr$i.filtered.shapeit2-duohmm-phased.snpid.vcf.gz --out ALL.1KG.GRCh38.chr$i &
			pwait $3
        	done
	fi

	wait
	gzip *.tped

else
        for i in {1..22}
        do
		mv CCDG\_14151\_B01\_GRM\_WGS\_2020-08-05\_chr$i.filtered.shapeit2-duohmm-phased.snpid.vcf.gz $2.1KG.GRCh38.chr$i.vcf.gz
	done
fi

