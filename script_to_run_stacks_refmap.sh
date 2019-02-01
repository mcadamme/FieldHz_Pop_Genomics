#This is the script I used to run the stacks ref_map pipeline for comparison of genotyping with samtools mpileup, then further filtering with vcftools.
#2/1/19
#MF

cd /media/megan/"New Volume"/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments

mkdir ./stacks_output

ref_map.pl -T 6 --samples ./ --popmap /home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/all_pops_no_reps_for_Stacks.txt -o ./stacks_output -X populations:"-r 0.75 -min_maf 0.05 --vcf --fasta_samples --vcf_haplotypes"
 
vcftools --vcf ./stacks_output/XXX.vcf --recode --out ./stacks_output/thinned_FieldHzea_variantsonly_stacks.vcf --minDP 3  --min-alleles 2 --max-alleles 2
