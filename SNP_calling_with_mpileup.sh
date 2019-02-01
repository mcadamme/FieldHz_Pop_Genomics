#script for SNP calling
#MF 11/1/2018
#bcftools v. 1.9


cd /media/megan/"New Volume"/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments

mkdir mpileupANDvcftools_output

bcftools mpileup -f ./Hzea_genome/GCA_002150865.1_Hzea_1.0_genomic.fna -b ./Hzea_BamFiles.txt -I -d 70 --threads 6 -O u -o ./mpileupANDvcftools_output/FieldHzea.bcf

bcftools call -vmO v -o ./mpileupANDvcftools_output/FieldHzea_variantsonly.vcf ./mpileupANDvcftools_output/FieldHzea.bcf #SNP calling and conversion to vcf format

#here I filtered my called SNPs in preparation for analysis in R

vcftools --vcf FieldHzea_variantsonly.vcf --recode --out thinned_FieldHzea_variantsonly.vcf --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75
 


