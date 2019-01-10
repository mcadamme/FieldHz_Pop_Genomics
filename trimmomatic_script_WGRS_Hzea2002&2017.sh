#script to trim off the low-quality ends of my reads with Trimmomatic v. 0.38.
#12/11/18
#MF

cd /media/megan/"New Volume"/Raw_Sequence_Data/H.zea_WGRS_DivelySamples/GouldFritz/

mkdir trimmed

#2002 samples
java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2002_1_S17_R1_001.fastq.gz 2002_1_S17_R2_001.fastq.gz ./trimmed/2002_1_S17_forward_paired.fq.gz ./trimmed/2002_1_S17_forward_unpaired.fq.gz ./trimmed/2002_1_S17_reverse_paired.fq.gz ./trimmed/2002_1_S17_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2002_2_S18_R1_001.fastq.gz 2002_2_S18_R2_001.fastq.gz ./trimmed/2002_2_S18_forward_paired.fq.gz ./trimmed/2002_2_S18_forward_unpaired.fq.gz ./trimmed/2002_2_S18_reverse_paired.fq.gz ./trimmed/2002_2_S18_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2002_3_S19_R1_001.fastq.gz 2002_3_S19_R2_001.fastq.gz ./trimmed/2002_3_S19_forward_paired.fq.gz ./trimmed/2002_3_S19_forward_unpaired.fq.gz ./trimmed/2002_3_S19_reverse_paired.fq.gz ./trimmed/2002_3_S19_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2002_4_S20_R1_001.fastq.gz 2002_4_S20_R2_001.fastq.gz ./trimmed/2002_4_S20_forward_paired.fq.gz ./trimmed/2002_4_S20_forward_unpaired.fq.gz ./trimmed/2002_4_S20_reverse_paired.fq.gz ./trimmed/2002_4_S20_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2002_5_S21_R1_001.fastq.gz 2002_5_S21_R2_001.fastq.gz ./trimmed/2002_5_S21_forward_paired.fq.gz ./trimmed/2002_5_S21_forward_unpaired.fq.gz ./trimmed/2002_5_S21_reverse_paired.fq.gz ./trimmed/2002_5_S21_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2002_6_S22_R1_001.fastq.gz 2002_6_S22_R2_001.fastq.gz ./trimmed/2002_6_S22_forward_paired.fq.gz ./trimmed/2002_6_S22_forward_unpaired.fq.gz ./trimmed/2002_6_S22_reverse_paired.fq.gz ./trimmed/2002_6_S22_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2002_7_S23_R1_001.fastq.gz 2002_7_S23_R2_001.fastq.gz ./trimmed/2002_7_S23_forward_paired.fq.gz ./trimmed/2002_7_S23_forward_unpaired.fq.gz ./trimmed/2002_7_S23_reverse_paired.fq.gz ./trimmed/2002_7_S23_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2002_10_S24_R1_001.fastq.gz 2002_10_S24_R2_001.fastq.gz ./trimmed/2002_10_S24_forward_paired.fq.gz ./trimmed/2002_10_S24_forward_unpaired.fq.gz ./trimmed/2002_10_S24_reverse_paired.fq.gz ./trimmed/2002_10_S24_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

#2017
java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2017_1_S25_R1_001.fastq.gz 2017_1_S25_R2_001.fastq.gz ./trimmed/2017_1_S25_forward_paired.fq.gz ./trimmed/2017_1_S25_forward_unpaired.fq.gz ./trimmed/2017_1_S25_reverse_paired.fq.gz ./trimmed/2017_1_S25_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2017_3_S26_R1_001.fastq.gz 2017_3_S26_R2_001.fastq.gz ./trimmed/2017_3_S26_forward_paired.fq.gz ./trimmed/2017_3_S26_forward_unpaired.fq.gz ./trimmed/2017_3_S26_reverse_paired.fq.gz ./trimmed/2017_3_S26_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2017_5_S27_R1_001.fastq.gz 2017_5_S27_R2_001.fastq.gz ./trimmed/2017_5_S27_forward_paired.fq.gz ./trimmed/2017_5_S27_forward_unpaired.fq.gz ./trimmed/2017_5_S27_reverse_paired.fq.gz ./trimmed/2017_5_S27_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2017_6_S28_R1_001.fastq.gz 2017_6_S28_R2_001.fastq.gz ./trimmed/2017_6_S28_forward_paired.fq.gz ./trimmed/2017_6_S28_forward_unpaired.fq.gz ./trimmed/2017_6_S28_reverse_paired.fq.gz ./trimmed/2017_6_S28_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2017_7_S29_R1_001.fastq.gz 2017_7_S29_R2_001.fastq.gz ./trimmed/2017_7_S29_forward_paired.fq.gz ./trimmed/2017_7_S29_forward_unpaired.fq.gz ./trimmed/2017_7_S29_reverse_paired.fq.gz ./trimmed/2017_7_S29_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2017_8_S30_R1_001.fastq.gz 2017_8_S30_R2_001.fastq.gz ./trimmed/2017_8_S30_forward_paired.fq.gz ./trimmed/2017_8_S30_forward_unpaired.fq.gz ./trimmed/2017_8_S30_reverse_paired.fq.gz ./trimmed/2017_8_S30_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2017_9_S31_R1_001.fastq.gz 2017_9_S31_R2_001.fastq.gz ./trimmed/2017_9_S31_forward_paired.fq.gz ./trimmed/2017_9_S31_forward_unpaired.fq.gz ./trimmed/2017_9_S31_reverse_paired.fq.gz ./trimmed/2017_9_S31_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 2017_10_S32_R1_001.fastq.gz 2017_10_S32_R2_001.fastq.gz ./trimmed/2017_10_S32_forward_paired.fq.gz ./trimmed/2017_10_S32_forward_unpaired.fq.gz ./trimmed/2017_10_S32_reverse_paired.fq.gz ./trimmed/2017_10_S32_reverse_unpaired.fq.gz SLIDINGWINDOW:5:20 MINLEN:50
