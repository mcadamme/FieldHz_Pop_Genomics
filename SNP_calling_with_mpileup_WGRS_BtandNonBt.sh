#SNP calling script for WGRS Bt and NonBt H.zea samples collected in MD.
#MF 1/15/2019
#bcftools v. 1.9


cd /media/megan/"New Volume1"/Hzea_WGRS_Bowtie2_output/BtandNonBt_alignmentFiles

mkdir WGRS_mpileupANDvcftools_output

bcftools mpileup -f /media/megan/"New Volume"/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/Hzea_genome/GCA_002150865.1_Hzea_1.0_genomic.fna -b ./WGRS_Hzea_BtandNonBt_BamFiles.txt -I -d 70 --threads 3 -O u -o ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea.bcf

bcftools call -vmO v -o ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly.vcf ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea.bcf #SNP calling and conversion to vcf format

#filtering my called SNPs in preparation for analysis

vcftools --vcf ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly.vcf --recode --out ./WGRS_mpileupANDvcftools_output/thinned_BtandNonBt_Hzea_variantsonly.vcf --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

#Getting genome-wide Allele frequencies for each pop
vcftools --vcf ./WGRS_mpileupANDvcftools_output/thinned_BtandNonBt_Hzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/Bt_WGRS_Hzea.txt  --out  ./WGRS_mpileupANDvcftools_output/thinned_FieldHzea_variantsonly_BtOnly.freq --freq

vcftools --vcf ./WGRS_mpileupANDvcftools_output/thinned_BtandNonBt_Hzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/NonBt_WGRS_Hzea.txt  --out  ./WGRS_mpileupANDvcftools_output/thinned_FieldHzea_variantsonly_NonBtOnly.freq --freq

#Getting allele freqs just for the outlier chromosomes

mkdir ./WGRS_mpileupANDvcftools_output/KZ118395.1_only

vcftools --vcf ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly.vcf --recode --out ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/KZ118395.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118395.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

vcftools --vcf ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/KZ118395.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/Bt_WGRS_Hzea.txt  --out  ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/thinned_FieldHzea_variantsonly_BtOnly.freq --freq

vcftools --vcf ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/KZ118395.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/NonBt_WGRS_Hzea.txt  --out  ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/thinned_FieldHzea_variantsonly_NonBtOnly.freq --freq

vcftools --vcf ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/KZ118395.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/Bt_WGRS_Hzea.txt  --out  ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/thinned_FieldHzea_variantsonly_BtOnly --het

vcftools --vcf ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/KZ118395.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/NonBt_WGRS_Hzea.txt  --out  ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/thinned_FieldHzea_variantsonly_NonBtOnly --het

/home/megan/src/vcflib/bin/wcFst --target 0,1,2,3,4,5,6,7  --background 8,9,10,11,12,13,14,15 --file ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/KZ118395.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./WGRS_mpileupANDvcftools_output/KZ118395.1_only/KZ118395.1_wcFst_BtandNonBt

#sliding window pfst analysis

/home/megan/src/vcflib/bin/pFst --target 0,1,2,3,4,5,6,7  --background 8,9,10,11,12,13,14,15 --file ./WGRS_mpileupANDvcftools_output/thinned_BtandNonBt_Hzea_variantsonly.vcf.recode.vcf --type PL > ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly_pFST_all

/home/megan/src/vcflib/bin/smoother --file ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly_pFST_all -o pFst > ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly_all.smoothed


#sliding window wcfst analysis

/home/megan/src/vcflib/bin/wcFst --target 0,1,2,3,4,5,6,7  --background 8,9,10,11,12,13,14,15 --file ./WGRS_mpileupANDvcftools_output/thinned_BtandNonBt_Hzea_variantsonly.vcf.recode.vcf --type PL > ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly_wcFST_all

/home/megan/src/vcflib/bin/smoother --file ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly_wcFST_all -o wcFst  -w 10000 -s 1000 > ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly_wcST_10kb_all.smoothed 

/home/megan/src/vcflib/bin/smoother --file ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly_wcFST_all -o wcFst  -w 40000 -s 10000 > ./WGRS_mpileupANDvcftools_output/BtandNonBt_Hzea_variantsonly_wcST_40kb_all.smoothed 



