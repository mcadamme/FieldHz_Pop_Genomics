#Filtering out top loading DAPC SNPs
#06162020 MF

cd /media/megan/"New Volume1"/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/mpileupANDvcftools_output

#top 1
vcftools --vcf thinned_FieldHzea_allpops.recode.vcf --exclude-positions ./DAPC2_top1per.txt --out top1_varContr_FieldHzea --recode

#top 5
vcftools --vcf thinned_FieldHzea_allpops.recode.vcf --exclude-positions ./DAPC2_top5per.txt --out top5_varContr_FieldHzea --recode

#top 2.5
vcftools --vcf thinned_FieldHzea_allpops.recode.vcf --exclude-positions ./DAPC2_top2p5per.txt --out top2p5_varContr_FieldHzea --recode

