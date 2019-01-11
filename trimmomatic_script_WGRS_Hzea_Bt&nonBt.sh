#script to trim off the low-quality ends of my Dively reads with Trimmomatic v. 0.38.
#01/08/19
#MF

cd /media/megan/"New Volume"/Raw_Sequence_Data/H.zea_WGRS_DivelySamples/GouldFritz/

mkdir trimmed_BtandNonBt

#Bt samples
java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Bt_1_S1_R1_001.fastq.gz Bt_1_S1_R2_001.fastq.gz ./trimmed_BtandNonBt/Bt_1_S1_forward_paired.fq.gz ./trimmed_BtandNonBt/Bt_1_S1_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Bt_1_S1_reverse_paired.fq.gz ./trimmed_BtandNonBt/Bt_1_S1_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Bt_2_S2_R1_001.fastq.gz Bt_2_S2_R2_001.fastq.gz ./trimmed_BtandNonBt/Bt_2_S2_forward_paired.fq.gz ./trimmed_BtandNonBt/Bt_2_S2_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Bt_2_S2_reverse_paired.fq.gz ./trimmed_BtandNonBt/Bt_2_S2_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Bt_3_S3_R1_001.fastq.gz Bt_3_S3_R2_001.fastq.gz ./trimmed_BtandNonBt/Bt_3_S3_forward_paired.fq.gz ./trimmed_BtandNonBt/Bt_3_S3_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Bt_3_S3_reverse_paired.fq.gz ./trimmed_BtandNonBt/Bt_3_S3_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Bt_4_S4_R1_001.fastq.gz Bt_4_S4_R2_001.fastq.gz ./trimmed_BtandNonBt/Bt_4_S4_forward_paired.fq.gz ./trimmed_BtandNonBt/Bt_4_S4_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Bt_4_S4_reverse_paired.fq.gz ./trimmed_BtandNonBt/Bt_4_S4_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Bt_5_S5_R1_001.fastq.gz Bt_5_S5_R2_001.fastq.gz ./trimmed_BtandNonBt/Bt_5_S5_forward_paired.fq.gz ./trimmed_BtandNonBt/Bt_5_S5_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Bt_5_S5_reverse_paired.fq.gz ./trimmed_BtandNonBt/Bt_5_S5_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Bt_6_S6_R1_001.fastq.gz Bt_6_S6_R2_001.fastq.gz ./trimmed_BtandNonBt/Bt_6_S6_forward_paired.fq.gz ./trimmed_BtandNonBt/Bt_6_S6_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Bt_6_S6_reverse_paired.fq.gz ./trimmed_BtandNonBt/Bt_6_S6_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Bt_7_S7_R1_001.fastq.gz Bt_7_S7_R2_001.fastq.gz ./trimmed_BtandNonBt/Bt_7_S7_forward_paired.fq.gz ./trimmed_BtandNonBt/Bt_7_S7_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Bt_7_S7_reverse_paired.fq.gz ./trimmed_BtandNonBt/Bt_7_S7_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Bt_8_S8_R1_001.fastq.gz Bt_8_S8_R2_001.fastq.gz ./trimmed_BtandNonBt/Bt_8_S8_forward_paired.fq.gz ./trimmed_BtandNonBt/Bt_8_S8_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Bt_8_S8_reverse_paired.fq.gz ./trimmed_BtandNonBt/Bt_8_S8_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

#Non-Bt samples
java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Non-Bt_1_S9_R1_001.fastq.gz Non-Bt_1_S9_R2_001.fastq.gz ./trimmed_BtandNonBt/Non-Bt_1_S9_forward_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_1_S9_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Non-Bt_1_S9_reverse_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_1_S9_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Non-Bt_2_S10_R1_001.fastq.gz Non-Bt_2_S10_R2_001.fastq.gz ./trimmed_BtandNonBt/Non-Bt_2_S10_forward_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_2_S10_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Non-Bt_2_S10_reverse_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_2_S10_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Non-Bt_3_S11_R1_001.fastq.gz Non-Bt_3_S11_R2_001.fastq.gz ./trimmed_BtandNonBt/Non-Bt_3_S11_forward_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_3_S11_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Non-Bt_3_S11_reverse_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_3_S11_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Non-Bt_4_S12_R1_001.fastq.gz Non-Bt_4_S12_R2_001.fastq.gz ./trimmed_BtandNonBt/Non-Bt_4_S12_forward_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_4_S12_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Non-Bt_4_S12_reverse_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_4_S12_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Non-Bt_5_S13_R1_001.fastq.gz Non-Bt_5_S13_R2_001.fastq.gz ./trimmed_BtandNonBt/Non-Bt_5_S13_forward_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_5_S13_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Non-Bt_5_S13_reverse_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_5_S13_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Non-Bt_6_S14_R1_001.fastq.gz Non-Bt_6_S14_R2_001.fastq.gz ./trimmed_BtandNonBt/Non-Bt_6_S14_forward_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_6_S14_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Non-Bt_6_S14_reverse_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_6_S14_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Non-Bt_7_S15_R1_001.fastq.gz Non-Bt_7_S15_R2_001.fastq.gz ./trimmed_BtandNonBt/Non-Bt_7_S15_forward_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_7_S15_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Non-Bt_7_S15_reverse_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_7_S15_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 Non-Bt_8_S16_R1_001.fastq.gz Non-Bt_8_S16_R2_001.fastq.gz ./trimmed_BtandNonBt/Non-Bt_8_S16_forward_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_8_S16_forward_unpaired.fq.gz ./trimmed_BtandNonBt/Non-Bt_8_S16_reverse_paired.fq.gz ./trimmed_BtandNonBt/Non-Bt_8_S16_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50
