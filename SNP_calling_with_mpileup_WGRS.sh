#script for SNP calling
#MF 01/20/2020
#bcftools v. 1.9


cd /media/megan/"New Volume"/Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020

mkdir WGRS_mpileupANDvcftools_output_v1

bcftools mpileup -f /media/megan/"New Volume"/Hzea_genome/GCA_002150865.1_Hzea_1.0_genomic.fna -b ./WGRS_Hzea_BamFiles.txt -I -d 70 --threads 6 -O u -o /media/megan/"New Volume"/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/FieldHzea.bcf


cd /media/megan/"New Volume"/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1

bcftools call -vmO v -o ./FieldHzea_variantsonly.vcf ./FieldHzea.bcf #SNP calling and conversion to vcf format

#filtering my called SNPs in preparation for analysis
vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./thinned_FieldHzea_variantsonly.vcf --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels


#running vcftools on bcftools output to filter and get pop stats.
#I first examined contigs with important candidate genes for Bt resistance to identify markers present, if any. 


#first looking only at KZ118424.1 (with tetraspannin) with some filtering
mkdir KZ118424.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ118424.1_only/FieldHzea2002 --chr KZ118424.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ118424.1_only/FieldHzea2012 --chr KZ118424.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ118424.1_only/FieldHzea2017 --chr KZ118424.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ118424.1_only/KZ118424.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118424.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ118424.1_only/KZ118424.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ118424.1_only/KZ118424.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ118424.1_only/KZ118424.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ118424.1_only/KZ118424.1_wcFst_2002_2017.smoothed

#looking at KZ118817.1 (with map4k4) with some filtering
mkdir KZ118817.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ118817.1_only/FieldHzea2002 --chr KZ118817.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ118817.1_only/FieldHzea2012 --chr KZ118817.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ118817.1_only/FieldHzea2017 --chr KZ118817.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ118817.1_only/KZ118817.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118817.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ118817.1_only/KZ118817.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ118817.1_only/KZ118817.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ118817.1_only/KZ118817.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ118817.1_only/KZ118817.1_wcFst_2002_2017.smoothed

                                                                                                                                                            

#looking at KZ118590.1 with ABCB1
mkdir KZ118590.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ118590.1_only/FieldHzea2002 --chr KZ118590.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ118590.1_only/FieldHzea2012 --chr KZ118590.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ118590.1_only/FieldHzea2017 --chr KZ118590.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ118590.1_only/KZ118590.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118590.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ118590.1_only/KZ118590.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ118590.1_only/KZ118590.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ118590.1_only/KZ118590.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ118590.1_only/KZ118590.1_wcFst_2002_2017.smoothed



#looking at KZ118207.1 with ABCA2
mkdir KZ118207.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ118207.1_only/FieldHzea2002 --chr KZ118207.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ118207.1_only/FieldHzea2012 --chr KZ118207.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ118207.1_only/FieldHzea2017 --chr KZ118207.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ118207.1_only/KZ118207.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118207.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ118207.1_only/KZ118207.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ118207.1_only/KZ118207.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ118207.1_only/KZ118207.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ118207.1_only/KZ118207.1_wcFst_2002_2017.smoothed



#looking at KZ118195.1 with CAD2
mkdir KZ118195.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ118195.1_only/FieldHzea2002 --chr KZ118195.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ118195.1_only/FieldHzea2012 --chr KZ118195.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ118195.1_only/FieldHzea2017 --chr KZ118195.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ118195.1_only/KZ118195.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118195.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ118195.1_only/KZ118195.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ118195.1_only/KZ118195.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ118195.1_only/KZ118195.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ118195.1_only/KZ118195.1_wcFST_2002_2017.smoothed



#looking only at KZ117463.1_only (with Cad86c) with some filtering
mkdir KZ117463.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ117463.1_only/FieldHzea2002 --chr KZ117463.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ117463.1_only/FieldHzea2012 --chr KZ117463.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ117463.1_only/FieldHzea2017 --chr KZ117463.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ117463.1_only/KZ117463.1_thinned_FieldHzea_variantsonly.vcf --chr KZ117463.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ117463.1_only/KZ117463.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ117463.1_only/KZ117463.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ117463.1_only/KZ117463.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ117463.1_only/KZ117463.1_wcFst_2002_2017.smoothed



#looking only at KZ118297.1_only (with ABCC2) with some filtering
mkdir KZ118297.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ118297.1_only/FieldHzea2002 --chr KZ118297.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ118297.1_only/FieldHzea2012 --chr KZ118297.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ118297.1_only/FieldHzea2017 --chr KZ118297.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ118297.1_only/KZ118297.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118297.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ118297.1_only/KZ118297.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ118297.1_only/KZ118297.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ118297.1_only/KZ118297.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ118297.1_only/KZ118297.1_wcFst_2002_2017.smoothed


#looking only at KZ117563.1_only (with HevCaLP orthologue) with some filtering
mkdir KZ117563.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ117563.1_only/FieldHzea2002 --chr KZ117563.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ117563.1_only/FieldHzea2012 --chr KZ117563.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ117563.1_only/FieldHzea2017 --chr KZ117563.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ117563.1_only/KZ117563.1_thinned_FieldHzea_variantsonly.vcf --chr KZ117563.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ117563.1_only/KZ117563.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ117563.1_only/KZ117563.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ117563.1_only/KZ117563.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ117563.1_only/KZ117563.1_wcFst_2002_2017.smoothed


#looking only at KZ117832.1_only (with ALP2) with some filtering

mkdir KZ117832.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ117832.1_only/FieldHzea2002 --chr KZ117832.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ117832.1_only/FieldHzea2012 --chr KZ117832.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ117832.1_only/FieldHzea2017 --chr KZ117832.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ117832.1_only/KZ117832.1_thinned_FieldHzea_variantsonly.vcf --chr KZ117832.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ117832.1_only/KZ117832.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ117832.1_only/KZ117832.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ117832.1_only/KZ117832.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ117832.1_only/KZ117832.1_wcFst_2002_2017.smoothed


#looking only at KZ118301.1_only (with APN) with some filtering
mkdir KZ118301.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ118301.1_only/FieldHzea2002 --chr KZ118301.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ118301.1_only/FieldHzea2012 --chr KZ118301.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ118301.1_only/FieldHzea2017 --chr KZ118301.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ118301.1_only/KZ118301.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118301.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ118301.1_only/KZ118301.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ118301.1_only/KZ118301.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ118301.1_only/KZ118301.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ118301.1_only/KZ118301.1_wcFst_2002_2017.smoothed



#looking at KZ118395.1 with carboxypeptidaseQ
mkdir KZ118395.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ118395.1_only/FieldHzea2002 --chr KZ118395.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ118395.1_only/FieldHzea2012 --chr KZ118395.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ118395.1_only/FieldHzea2017 --chr KZ118395.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ118395.1_only/KZ118395.1_thinned_FieldHzea_variantsonly.vcf --chr KZ118395.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ118395.1_only/KZ118395.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ118395.1_only/KZ118395.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ118395.1_only/KZ118395.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ118395.1_only/KZ118395.1_wcFst_2002_2017.smoothed


#looking at KZ117131.1 with cyp333b
mkdir KZ117131.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --out ./KZ117131.1_only/FieldHzea2002 --chr KZ117131.1 --freq                                                                                                         
                                                                                                                                                      
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2012.txt --out ./KZ117131.1_only/FieldHzea2012 --chr KZ117131.1 --freq

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --out ./KZ117131.1_only/FieldHzea2017 --chr KZ117131.1 --freq

vcftools --vcf ./FieldHzea_variantsonly.vcf --recode --out ./KZ117131.1_only/KZ117131.1_thinned_FieldHzea_variantsonly.vcf --chr KZ117131.1 --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --remove-indels

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./KZ117131.1_only/KZ117131.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./KZ117131.1_only/KZ117131.1_wcFst_2002_2017

/home/megan/src/vcflib/bin/smoother --file ./KZ117131.1_only/KZ117131.1_wcFst_2002_2017 -o wcFst -w 10000 > ./KZ117131.1_only/KZ117131.1_wcFst_2002_2017.smoothed



#whole genome sliding window wcfst analysis 2002 & 2017 - diff window sizes

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./2002and2017_wcFST_all

/home/megan/src/vcflib/bin/smoother --file ./2002and2017_wcFST_all -o wcFst -w 5000 > ./2002and2017_5kb_wcFST_all.smoothed

/home/megan/src/vcflib/bin/smoother --file ./2002and2017_wcFST_all -o wcFst -w 10000 > ./2002and2017_10kb_wcFST_all.smoothed

/home/megan/src/vcflib/bin/smoother --file ./2002and2017_wcFST_all -o wcFst -w 20000 -s 5000> ./2002and2017_20kb_wcFST_all.smoothed

/home/megan/src/vcflib/bin/smoother --file ./2002and2017_wcFST_all -o wcFst -w 40000 -s 5000 > ./2002and2017_40kb_wcFST_all.smoothed

#whole genome sliding window wcfst analysis 2002 & 2012 - diff window sizes

/home/megan/src/vcflib/bin/wcFst --target 13,14,15,16,17,18,19,20,21,22,23  --background 0,1,2,3,4,5,6,7,8,9,10,11,12 --file ./thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./2002and2012_wcFST_all

/home/megan/src/vcflib/bin/smoother --file ./2002and2012_wcFST_all -o wcFst -w 5000 > ./2002and2012_5kb_wcFST_all.smoothed

/home/megan/src/vcflib/bin/smoother --file ./2002and2012_wcFST_all -o wcFst -w 10000 > ./2002and2012_10kb_wcFST_all.smoothed

/home/megan/src/vcflib/bin/smoother --file ./2002and2012_wcFST_all -o wcFst -w 20000 -s 5000 > ./2002and2012_20kb_wcFST_all.smoothed

/home/megan/src/vcflib/bin/smoother --file ./2002and2012_wcFST_all -o wcFst -w 40000 -s 5000 > ./2002and2012_40kb_wcFST_all.smoothed

#whole genome sliding window wcfst analysis 2012 & 2017 - diff window sizes

/home/megan/src/vcflib/bin/wcFst --target 24,25,26,27,28,29,30,31,32,33,34  --background 13,14,15,16,17,18,19,20,21,22,23 --file ./thinned_FieldHzea_variantsonly.vcf.recode.vcf --type PL > ./2012and2017_wcFST_all

/home/megan/src/vcflib/bin/smoother --file ./2012and2017_wcFST_all -o wcFst -w 5000 > ./2012and2017_5kb_wcFST_all.smoothed

/home/megan/src/vcflib/bin/smoother --file ./2012and2017_wcFST_all -o wcFst -w 10000 > ./2012and2017_10kb_wcFST_all.smoothed

/home/megan/src/vcflib/bin/smoother --file ./2012and2017_wcFST_all -o wcFst -w 20000 -s 5000 > ./2012and2017_20kb_wcFST_all.smoothed

/home/megan/src/vcflib/bin/smoother --file ./2012and2017_wcFST_all -o wcFst -w 40000 -s 5000 > ./2012and2017_20kb_wcFST_all.smoothed


