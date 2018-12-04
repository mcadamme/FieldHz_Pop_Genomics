#to get comparisons by population

cd Downloads

vcftools --vcf FieldHzea_variantsonly.vcf --out thinned_FieldHzea_variantsonly.vcf --recode-info-all --minDP 3 --maf 0.05 --max-missing 0.5



#2002 and 2007
bcftools view Filtered_FieldHzea.recode.vcf -o Filtered_FieldHzea2002and2007.vcf -S /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2007.txt

vcftools --vcf Filtered_FieldHzea2002and2007.vcf --out thinned_FieldHzea2002and2007 --recode --minDP 3 --maf 0.05 --max-missing 0.5 --remove-indels

#2002 and 2012
bcftools view Filtered_FieldHzea.recode.vcf -o Filtered_FieldHzea2002and2012.vcf -S /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2012.txt --force-samples

vcftools --vcf Filtered_FieldHzea2002and2012.vcf --out thinned_FieldHzea2002and2012 --recode --minDP 3 --maf 0.05 --max-missing 0.5 --remove-indels

#2002 and 2016
bcftools view Filtered_FieldHzea.recode.vcf -o Filtered_FieldHzea2002and2016.vcf -S /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2016.txt --force-samples

vcftools --vcf Filtered_FieldHzea2002and2016.vcf --out thinned_FieldHzea2002and2016 --recode --minDP 3 --maf 0.05 --max-missing 0.5 --remove-indels


