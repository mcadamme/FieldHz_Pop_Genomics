#This is the script I used to run the stacks ref_map pipeline for comparison of genotyping with samtools mpileup, then further filtering with vcftools.
#2/1/19
#MF

cd /media/megan/"New Volume"/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments

#mkdir ./stacks_output

ref_map.pl -T 6 --samples ./ --popmap /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/all_pops_no_reps_for_Stacks.txt -o ./stacks_output -X "gstacks: --min-mapq 5" -X "populations: --vcf --fasta_samples --vcf_haplotypes"
 
vcftools --vcf ./stacks_output/populations.snps.vcf --recode --out ./stacks_output/thinned_FieldHzea_variantsonly_stacks.vcf --minDP 3  --min-alleles 2 --max-alleles 2 --maf 0.05 --max-missing 0.75 --keep /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/allpops_no_low_stacks.txt

#Rerunning populations to get haplotype information.
