#SNP calling script for MD and LA samples.
#MF 7/22/2020
#bcftools v. 1.9


cd /media/megan/"New Volume1"/Hzea_WGRS_Bowtie2_output/two_state_2017_samps

#LA has year in name, MD has Bt and non-Bt

ls *.bam > WGRS_Hzea_2017_BamFiles.txt

bcftools mpileup -f /media/megan/"New Volume1"/Hzea_genome/GCA_002150865.1_Hzea_1.0_genomic.fna -b ./WGRS_Hzea_2017_BamFiles.txt -I -d 70 --threads 3 -O u -o /home/megan/Hzea_2017_samps.bcf

bcftools call -vmO v -o /home/megan/Hzea_2017_samps_variantsOnly.vcf /home/megan/Hzea_2017_samps.bcf #SNP calling and conversion to vcf format

#filtering my called SNPs in preparation for analysis

vcftools --vcf /home/megan/Hzea_2017_samps_variantsOnly.vcf --recode --out ./Hzea_2017_samps_variantsonly.vcf --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --thin 25000 --remove-indels #produces 14697 SNPs

R CMD BATCH /home/megan/scripts/Field_HZ_Pop_Genomics/pairwise_fst_LA_MD_2017.R


