#script to align H.zea WGRS reads with bowtie2 v.2.2.6
#01/14/2020
#MLF

cd /media/megan/"New Volume"/

mkdir Hzea_WGRS_Bowtie2_output
mkdir ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020

#end-to-end alignments with highest sensitivity

##2002
bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_1_S1_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_1_S1_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_01.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_2_S2_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_2_S2_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_02.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_3_S3_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_3_S3_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_03.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_4_S4_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_4_S4_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_04.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_5_S5_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_5_S5_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_05.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_6_S6_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_6_S6_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_06.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_7_S7_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_7_S7_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_07.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_8_S8_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_8_S8_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_08.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_9_S9_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_9_S9_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_09.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_10_S10_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_10_S10_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_10.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_11_S11_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_11_S11_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_11.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_12_S12_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_12_S12_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_12.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_13_S13_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2002_13_S13_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2002_13.sam

##2012
bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_1_S17_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_1_S17_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_01.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_2_S18_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_2_S18_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_02.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_3_S19_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_3_S19_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_03.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_4_S20_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_4_S20_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_04.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_5_S21_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_5_S21_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_05.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_6_S22_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_6_S22_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_06.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_7_S23_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2012_7_S23_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_07.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2012_9_S14_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2012_9_S14_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_09.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2012_10_S15_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2012_10_S15_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_10.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2012_11_S16_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2012_11_S16_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_11.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2012_12_S17_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2012_12_S17_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2012_12.sam

#2017
bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_1_S25_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_1_S25_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_01.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_3_S26_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_3_S26_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_02.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_5_S27_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_5_S27_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_03.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_6_S28_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_6_S28_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_04.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_7_S29_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_7_S29_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_05.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_8_S30_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_8_S30_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_06.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_9_S31_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/2017_9_S31_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_07.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2017_9_S19_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2017_9_S19_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_08.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2017_10_S20_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2017_10_S20_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_09.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2017_11_S21_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2017_11_S21_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_10.sam

bowtie2 -x ./Hzea_genome/Hzea_genome -1 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2017_12_S22_forward_paired.fastq.gz -2 ./Raw_Sequence_Data/filter_trimmed_Hzea_WGRS_Field_allyears/AMD_2017_12_S22_reverse_paired.fastq.gz --very-sensitive --threads 8 -S ./Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020/HZ_2017_11.sam

