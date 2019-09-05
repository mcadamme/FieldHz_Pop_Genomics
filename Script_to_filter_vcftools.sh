#running vcftools to filter and get pop stats.

#first for BCFtools output

cd /media/megan/"New Volume"/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/mpileupANDvcftools_output


#first looking only at KZ118424.1 (with tetraspannin) with some filtering
mkdir KZ118424.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2002 --chr KZ118424.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2002 --chr KZ118424.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2002 --chr KZ118424.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2007 --chr KZ118424.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2007 --chr KZ118424.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2007 --chr KZ118424.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2012 --chr KZ118424.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2012 --chr KZ118424.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2012 --chr KZ118424.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2016 --chr KZ118424.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2016 --chr KZ118424.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118424.1_only/FieldHzea2016 --chr KZ118424.1 --geno-depth

#looking at KZ117720.1
mkdir KZ117720.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2002 --chr KZ117720.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2002 --chr KZ117720.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2002 --chr KZ117720.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2007 --chr KZ117720.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2007 --chr KZ117720.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2007 --chr KZ117720.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2012 --chr KZ117720.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2012 --chr KZ117720.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2012 --chr KZ117720.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2016 --chr KZ117720.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2016 --chr KZ117720.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117720.1_only/FieldHzea2016 --chr KZ117720.1 --geno-depth

#looking at KZ118207.1
mkdir KZ118207.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2002 --chr KZ118207.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2002 --chr KZ118207.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2002 --chr KZ118207.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2007 --chr KZ118207.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2007 --chr KZ118207.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2007 --chr KZ118207.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2012 --chr KZ118207.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2012 --chr KZ118207.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2012 --chr KZ118207.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2016 --chr KZ118207.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2016 --chr KZ118207.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118207.1_only/FieldHzea2016 --chr KZ118207.1 --geno-depth

#looking at KZ118195.1
mkdir KZ118195.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2002 --chr KZ118195.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2002 --chr KZ118195.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2002 --chr KZ118195.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2007 --chr KZ118195.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2007 --chr KZ118195.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2007 --chr KZ118195.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2012 --chr KZ118195.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2012 --chr KZ118195.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2012 --chr KZ118195.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2016 --chr KZ118195.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2016 --chr KZ118195.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118195.1_only/FieldHzea2016 --chr KZ118195.1 --geno-depth

#looking only at KZ117463.1_only (with Cad86c) with some filtering

mkdir KZ117463.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2002 --chr KZ117463.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2002 --chr KZ117463.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2002 --chr KZ117463.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2007 --chr KZ117463.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2007 --chr KZ117463.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2007 --chr KZ117463.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2012 --chr KZ117463.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2012 --chr KZ117463.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2012 --chr KZ117463.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2016 --chr KZ117463.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2016 --chr KZ117463.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117463.1_only/FieldHzea2016 --chr KZ117463.1 --geno-depth


#looking only at KZ118297.1_only (with ABCC2) with some filtering

mkdir KZ118297.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2002 --chr KZ118297.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2002 --chr KZ118297.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2002 --chr KZ118297.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2007 --chr KZ118297.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2007 --chr KZ118297.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2007 --chr KZ118297.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2012 --chr KZ118297.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2012 --chr KZ118297.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2012 --chr KZ118297.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2016 --chr KZ118297.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2016 --chr KZ118297.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118297.1_only/FieldHzea2016 --chr KZ118297.1 --geno-depth


#looking only at KZ117563.1_only (with HevCaLP orthologue) with some filtering

mkdir KZ117563.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2002 --chr KZ117563.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2002 --chr KZ117563.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2002 --chr KZ117563.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2007 --chr KZ117563.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2007 --chr KZ117563.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2007 --chr KZ117563.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2012 --chr KZ117563.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2012 --chr KZ117563.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2012 --chr KZ117563.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2016 --chr KZ117563.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2016 --chr KZ117563.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117563.1_only/FieldHzea2016 --chr KZ117563.1 --geno-depth


#looking only at KZ117832.1_only (with ALP2) with some filtering

mkdir KZ117832.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2002 --chr KZ117832.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2002 --chr KZ117832.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2002 --chr KZ117832.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2007 --chr KZ117832.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2007 --chr KZ117832.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2007 --chr KZ117832.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2012 --chr KZ117832.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2012 --chr KZ117832.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2012 --chr KZ117832.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2016 --chr KZ117832.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2016 --chr KZ117832.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ117832.1_only/FieldHzea2016 --chr KZ117832.1 --geno-depth

#looking only at KZ118301.1_only (with APN) with some filtering

mkdir KZ118301.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2002 --chr KZ118301.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2002 --chr KZ118301.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2002 --chr KZ118301.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2007 --chr KZ118301.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2007 --chr KZ118301.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2007 --chr KZ118301.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2012 --chr KZ118301.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2012 --chr KZ118301.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2012 --chr KZ118301.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2016 --chr KZ118301.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2016 --chr KZ118301.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118301.1_only/FieldHzea2016 --chr KZ118301.1 --geno-depth

#looking only at KZ118395.1_only (strong outlier) with some filtering

mkdir KZ118395.1_only

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2002 --chr KZ118395.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2002 --chr KZ118395.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2002 --chr KZ118395.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2007 --chr KZ118395.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2007 --chr KZ118395.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2007 --chr KZ118395.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2012 --chr KZ118395.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2012 --chr KZ118395.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2012 --chr KZ118395.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2016 --chr KZ118395.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2016 --chr KZ118395.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --out ./KZ118395.1_only/FieldHzea2016 --chr KZ118395.1 --geno-depth

###Looking at all scaffolds
#2002 samples
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --thin 200 --out FieldHzea2002 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --thin 200 --out FieldHzea2002 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --thin 200 --out FieldHzea2002 --site-pi
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --thin 200 --out thinned_FieldHzea2002 --recode


#2007 samples
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --thin 200 --out FieldHzea2007 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --thin 200 --out FieldHzea2007 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --thin 200 --out FieldHzea2007 --site-pi
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt --thin 200 --out thinned_FieldHzea2007 --recode

#2012 samples
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --thin 200 --out FieldHzea2012 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --thin 200 --out FieldHzea2012 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --thin 200 --out FieldHzea2012 --site-pi
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt --thin 200 --out thinned_FieldHzea2012 --recode

#2016 samples
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --thin 200 --out FieldHzea2016 --freq
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --thin 200 --out FieldHzea2016 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --thin 200 --out FieldHzea2016 --site-pi
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt --thin 200 --out thinned_FieldHzea2016 --recode


#run of vcftools version of pairwise fst per site - will check if this makes sense with R analysis.
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --out 2002and2007 --thin 200 --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --out 2002and2012 --thin 200 --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --out 2002and2016 --thin 200 --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt

##getting pairwise datasets for R analysis
##generating combined popfiles
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2007.txt
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2012.txt
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2016.txt
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007and2012.txt
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007and2016.txt
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012and2016.txt

cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2007.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012and2016.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/allpops_no_low.txt



vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --out thinned_FieldHzea2002and2007 --thin 200 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2007.txt --recode

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --out thinned_FieldHzea2002and2012 --thin 200 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2012.txt --recode

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --out thinned_FieldHzea2002and2016 --thin 200 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2016.txt --recode

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --out thinned_FieldHzea2007and2012 --thin 200 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007and2012.txt --recode

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --out thinned_FieldHzea2007and2016 --thin 200 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007and2016.txt --recode

vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --out thinned_FieldHzea2012and2016 --thin 200 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012and2016.txt --recode


#full set thinned and filtered.
vcftools --vcf thinned_FieldHzea_variantsonly.vcf.recode.vcf --thin 200 --out thinned_FieldHzea_allpops --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/allpops_no_low.txt --recode



###Now using Stacks dataset - see supplement

cd /media/megan/"New Volume"/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/stacks_output

#first looking only at KZ118424.1 (with tetraspannin) with some filtering
mkdir KZ118424.1_only

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2002 --chr KZ118424.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2002 --chr KZ118424.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2002 --chr KZ118424.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2007 --chr KZ118424.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2007 --chr KZ118424.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2007 --chr KZ118424.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2012 --chr KZ118424.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2012 --chr KZ118424.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2012 --chr KZ118424.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2016 --chr KZ118424.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2016 --chr KZ118424.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118424.1_only/FieldHzea2016 --chr KZ118424.1 --geno-depth

#looking at KZ117720.1
mkdir KZ117720.1_only

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2002 --chr KZ117720.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2002 --chr KZ117720.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2002 --chr KZ117720.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2007 --chr KZ117720.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2007 --chr KZ117720.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2007 --chr KZ117720.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2012 --chr KZ117720.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2012 --chr KZ117720.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2012 --chr KZ117720.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2016 --chr KZ117720.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2016 --chr KZ117720.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117720.1_only/FieldHzea2016 --chr KZ117720.1 --geno-depth

#looking at KZ118207.1
mkdir KZ118207.1_only

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2002 --chr KZ118207.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2002 --chr KZ118207.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2002 --chr KZ118207.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2007 --chr KZ118207.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2007 --chr KZ118207.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2007 --chr KZ118207.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2012 --chr KZ118207.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2012 --chr KZ118207.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2012 --chr KZ118207.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2016 --chr KZ118207.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2016 --chr KZ118207.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118207.1_only/FieldHzea2016 --chr KZ118207.1 --geno-depth

#looking at KZ118195.1
mkdir KZ118195.1_only

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2002 --chr KZ118195.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2002 --chr KZ118195.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2002 --chr KZ118195.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2007 --chr KZ118195.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2007 --chr KZ118195.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2007 --chr KZ118195.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2012 --chr KZ118195.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2012 --chr KZ118195.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2012 --chr KZ118195.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2016 --chr KZ118195.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2016 --chr KZ118195.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118195.1_only/FieldHzea2016 --chr KZ118195.1 --geno-depth

#looking only at KZ117463.1_only (with Cad86c) with some filtering

mkdir KZ117463.1_only

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2002 --chr KZ117463.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2002 --chr KZ117463.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2002 --chr KZ117463.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2007 --chr KZ117463.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2007 --chr KZ117463.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2007 --chr KZ117463.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2012 --chr KZ117463.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2012 --chr KZ117463.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2012 --chr KZ117463.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2016 --chr KZ117463.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2016 --chr KZ117463.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117463.1_only/FieldHzea2016 --chr KZ117463.1 --geno-depth


#looking only at KZ118297.1_only (with ABCC2) with some filtering

mkdir KZ118297.1_only

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2002 --chr KZ118297.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2002 --chr KZ118297.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2002 --chr KZ118297.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2007 --chr KZ118297.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2007 --chr KZ118297.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2007 --chr KZ118297.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2012 --chr KZ118297.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2012 --chr KZ118297.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2012 --chr KZ118297.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2016 --chr KZ118297.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2016 --chr KZ118297.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118297.1_only/FieldHzea2016 --chr KZ118297.1 --geno-depth


#looking only at KZ117563.1_only (with HevCaLP orthologue) with some filtering

mkdir KZ117563.1_only

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2002 --chr KZ117563.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2002 --chr KZ117563.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2002 --chr KZ117563.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2007 --chr KZ117563.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2007 --chr KZ117563.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2007 --chr KZ117563.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2012 --chr KZ117563.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2012 --chr KZ117563.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2012 --chr KZ117563.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2016 --chr KZ117563.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2016 --chr KZ117563.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117563.1_only/FieldHzea2016 --chr KZ117563.1 --geno-depth


#looking only at KZ117832.1_only (with ALP2) with some filtering

mkdir KZ117832.1_only

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2002 --chr KZ117832.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2002 --chr KZ117832.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2002 --chr KZ117832.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2007 --chr KZ117832.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2007 --chr KZ117832.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2007 --chr KZ117832.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2012 --chr KZ117832.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2012 --chr KZ117832.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2012 --chr KZ117832.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2016 --chr KZ117832.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2016 --chr KZ117832.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ117832.1_only/FieldHzea2016 --chr KZ117832.1 --geno-depth

#looking only at KZ118301.1_only (with APN) with some filtering

mkdir KZ118301.1_only

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2002 --chr KZ118301.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2002 --chr KZ118301.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2002 --chr KZ118301.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2007 --chr KZ118301.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2007 --chr KZ118301.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2007 --chr KZ118301.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2012 --chr KZ118301.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2012 --chr KZ118301.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2012 --chr KZ118301.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2016 --chr KZ118301.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2016 --chr KZ118301.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118301.1_only/FieldHzea2016 --chr KZ118301.1 --geno-depth

#looking only at KZ118395.1_only (strong outlier) with some filtering

mkdir KZ118395.1_only

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2002 --chr KZ118395.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2002 --chr KZ118395.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2002 --chr KZ118395.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2007 --chr KZ118395.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2007 --chr KZ118395.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2007 --chr KZ118395.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2012 --chr KZ118395.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2012 --chr KZ118395.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2012 --chr KZ118395.1 --geno-depth
                                                                                                                                                             
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2016 --chr KZ118395.1 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2016 --chr KZ118395.1 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --out ./KZ118395.1_only/FieldHzea2016 --chr KZ118395.1 --geno-depth

###Looking at all scaffolds
#2002 samples
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2002 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2002 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2002 --site-pi
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --thin 50 --out thinned_FieldHzea2002 --recode


#2007 samples
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2007 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2007 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2007 --site-pi
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt --thin 50 --out thinned_FieldHzea2007 --recode

#2012 samples
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2012 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2012 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2012 --site-pi
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt --thin 50 --out thinned_FieldHzea2012 --recode

#2016 samples
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2016 --freq
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2016 --site-mean-depth
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --thin 50 --out FieldHzea2016 --site-pi
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt --thin 50 --out thinned_FieldHzea2016 --recode


#run of vcftools version of pairwisefst per site - will check if this makes sense with R analysis.
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --out 2002and2007 --thin 50 --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --out 2002and2012 --thin 50 --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --out 2002and2016 --thin 50 --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt --weir-fst-pop /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt

##getting pairwise datasets for R analysis
##generating combined popfiles
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2007_stacks.txt
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2012_stacks.txt
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002_noreps_or_low_stacks.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2016_stacks.txt
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007and2012_stacks.txt
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007_noreps_or_low_stacks.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007and2016_stacks.txt
cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012_noreps_or_low_stacks.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2016_noreps_or_low_stacks.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012and2016_stacks.txt

cat /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2007_stacks.txt /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012and2016_stacks.txt > /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/allpops_no_low_stacks.txt



vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --out thinned_FieldHzea2002and2007 --thin 50 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2007_stacks.txt --recode

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --out thinned_FieldHzea2002and2012 --thin 50 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2012_stacks.txt --recode

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --out thinned_FieldHzea2002and2016 --thin 50 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2016_stacks.txt --recode

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --out thinned_FieldHzea2007and2012 --thin 50 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007and2012_stacks.txt --recode

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --out thinned_FieldHzea2007and2016 --thin 50 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007and2016_stacks.txt --recode

vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --out thinned_FieldHzea2012and2016 --thin 50 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012and2016_stacks.txt --recode


#full set thinned and filtered.
vcftools --vcf thinned_FieldHzea_variantsonly_stacks.vcf.recode.vcf --thin 50 --out thinned_FieldHzea_allpops --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/allpops_no_low_stacks.txt --recode






