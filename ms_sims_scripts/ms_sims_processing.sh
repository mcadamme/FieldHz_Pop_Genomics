#This is the bash script I used to process my ms simulation file to get Fst values.  It relies on several R scripts, which must be accessible to calculate and plot FST values for each simulation.

#paths to files must change depending on which dataset is being used
#splitting files into unique genotype matrices

cd /PATH/TO/OUTPUT_DIR
#example - cd /home/megan/ms_sims/FieldPops/N0_140thou_10kb

cat /PATH/TO/IN_FILE | awk 'BEGIN {RS="//";FS="\n";OFS="\n"}{$1=$2=$3="";print $0 > NR ".txt"}'
#example - cat ./N0_140thou.txt | awk 'BEGIN {RS="//";FS="\n";OFS="\n"}{$1=$2=$3="";print $0 > NR ".txt"}'


#adding delimiter for each snp
for file in *.txt;
do
cat $file | gawk '{$1=$1}1' FPAT='.{1}' OFS=, > "$file"_output.txt
done

#removing empty files b/c some were generated, probably because of the nature of these sims.
find . -name "*.txt" -type 'f' -size -10c -delete


#removing few files that get renamed but mess up the analysis
rm /home/megan/ms_sims/FieldPops/N0_21thou_10kb/N0_21thou.txt_output.txt
rm /home/megan/ms_sims/FieldPops/N0_63thou_10kb/N0_63thou.txt_output.txt
rm /home/megan/ms_sims/FieldPops/N0_105thou_10kb/N0_105thou.txt_output.txt
rm /home/megan/ms_sims/FieldPops/N0_140thou_10kb/N0_140thou.txt_output.txt

#Next, generating FST distributions - can run these in parallel.

R CMD BATCH /home/megan/scripts/Field_HZ_Pop_Genomics/ms_sims_scripts/ms_sim_analysis_FSTvalsPart1.R
R CMD BATCH /home/megan/scripts/Field_HZ_Pop_Genomics/ms_sims_scripts/ms_sim_analysis_FSTvalsPart2.R
R CMD BATCH /home/megan/scripts/Field_HZ_Pop_Genomics/ms_sims_scripts/ms_sim_analysis_FSTvalsPart3.R
R CMD BATCH /home/megan/scripts/Field_HZ_Pop_Genomics/ms_sims_scripts/ms_sim_analysis_FSTvalsPart4.R



