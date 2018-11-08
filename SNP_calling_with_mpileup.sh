#script for SNP calling
#MF 11/1/2018
#bcftools v. 1.9


cd /media/megan/"New Volume"/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments

mkdir mpileupANDvcftools_output

bcftools mpileup -f ./Hzea_genome/GCA_002150865.1_Hzea_1.0_genomic.fna -b ./Hzea_BamFiles.txt -I -d 70 --threads 6 -O u -o ./mpileupANDvcftools_output/FieldHzea.bcf

bcftools view ./samtoolsANDvcftools_output/FieldHzea_calling.bcf > ./samtoolsANDvcftools_output/FieldHzea.vcf #conversion to vcf format


