#Generating VCFs for Cry-associated scaffolds
#MF 04/15/2021


cd ~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS


#first looking only at KZ118241.1 (Chr1) with some filtering
#mkdir KZ118241.1_only

vcftools --vcf ./vcf_files/FieldHzea_variantsonly.vcf --recode --out ./KZ118241.1_only/KZ118241.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118241.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.9 --remove-indels --thin 500

#KZ118067.1 (Chr6) with some filtering
#mkdir KZ118067.1_only

vcftools --vcf ./vcf_files/FieldHzea_variantsonly.vcf --recode --out ./KZ118067.1_only/KZ118067.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118067.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.9 --remove-indels --thin 500

#KZ116099.1 (Chr7) with some filtering
#mkdir KZ116099.1_only

vcftools --vcf ./vcf_files/FieldHzea_variantsonly.vcf --recode --out ./KZ116099.1_only/KZ116099.1_thinned_FieldHzea_variantsonly.vcf --chr KZ116099.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.9 --remove-indels --thin 500

#KZ117975.1 (Chr7) with some filtering
#mkdir KZ117975.1_only

vcftools --vcf ./vcf_files/FieldHzea_variantsonly.vcf --recode --out ./KZ117975.1_only/KZ117975.1_thinned_FieldHzea_variantsonly.vcf --chr KZ117975.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.9 --remove-indels --thin 500

#KZ117108.1 (Chr7) with some filtering - this one is CRY2!!!
#mkdir KZ117108.1_only

vcftools --vcf ./vcf_files/FieldHzea_variantsonly.vcf --recode --out ./KZ117108.1_only/KZ117108.1_thinned_FieldHzea_variantsonly.vcf --chr KZ117108.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.9 --remove-indels --thin 500


#KZ118133.1 (Chr11) with some filtering - left this one thinned every 250 bp because there was only 1 SNP left in the region after thinning to every 500
#mkdir KZ118133.1_only

vcftools --vcf ./vcf_files/FieldHzea_variantsonly.vcf --recode --out ./KZ118133.1_only/KZ118133.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118133.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.9 --remove-indels --thin 250


#KZ118395.1 (Chr5) with some filtering
#mkdir KZ118395.1_only

vcftools --vcf ./vcf_files/FieldHzea_variantsonly.vcf --recode --out ./KZ118395.1_only/KZ118395.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118395.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.9 --remove-indels --thin 500

#KZ117131.1 (Chr5) with some filtering
#mkdir KZ117131.1_only

vcftools --vcf ./vcf_files/FieldHzea_variantsonly.vcf --recode --out ./KZ117131.1_only/KZ117131.1_thinned_FieldHzea_variantsonly.vcf --chr KZ117131.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.9 --remove-indels --thin 500


