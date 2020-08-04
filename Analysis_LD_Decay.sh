#script to plot LD Decay by year
#08042020


cd /media/megan/"New Volume1"/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1

#getting vcf files for 2002 and 2017
vcftools --vcf ./thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2002.txt --recode --out ./thinned_FieldHzea_variantsonly_2002 --max-missing 0.75 --remove-indels --thin 500

vcftools --vcf ./thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/WGS_2017.txt --recode --out ./thinned_FieldHzea_variantsonly_2017 --max-missing 0.75 --remove-indels --thin 500



cd /media/megan/easystore

#Getting r2 values
/home/megan/src/plink_linux_x86_64_20200616/plink --vcf /media/megan/"New Volume1"/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/thinned_FieldHzea_variantsonly_2002.recode.vcf --r2 --double-id --allow-extra-chr --ld-window-r2 0 --ld-window 99999 --ld-window-kb 40 

/home/megan/src/plink_linux_x86_64_20200616/plink --vcf /media/megan/"New Volume1"/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v1/thinned_FieldHzea_variantsonly_2017.recode.vcf --r2 --double-id --allow-extra-chr --ld-window-r2 0 --ld-window 99999 --ld-window-kb 40 
