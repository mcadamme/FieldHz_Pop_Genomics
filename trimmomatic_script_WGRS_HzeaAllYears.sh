##Script to run trimmomatic for H. zea whole genome resequencing data for Megan. 08/02/19 KLB

cd /Volumes/data_a3/Hzea/
mkdir Trimmed

## 2002
java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_1_S1_R1_001.fastq.gz AMD_2002_1_S1_R2_001.fastq.gz ./Trimmed/AMD_2002_1_S1_forward_paired.fastq.gz ./Trimmed/AMD_2002_1_S1_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_1_S1_reverse_paired.fastq.gz ./Trimmed/AMD_2002_1_S1_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_2_S2_R1_001.fastq.gz AMD_2002_2_S2_R2_001.fastq.gz ./Trimmed/AMD_2002_2_S2_forward_paired.fastq.gz ./Trimmed/AMD_2002_2_S2_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_2_S2_reverse_paired.fastq.gz ./Trimmed/AMD_2002_2_S2_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_3_S3_R1_001.fastq.gz AMD_2002_3_S3_R2_001.fastq.gz ./Trimmed/AMD_2002_3_S3_forward_paired.fastq.gz ./Trimmed/AMD_2002_3_S3_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_3_S3_reverse_paired.fastq.gz ./Trimmed/AMD_2002_3_S3_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_4_S4_R1_001.fastq.gz AMD_2002_4_S4_R2_001.fastq.gz ./Trimmed/AMD_2002_4_S4_forward_paired.fastq.gz ./Trimmed/AMD_2002_4_S4_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_4_S4_reverse_paired.fastq.gz ./Trimmed/AMD_2002_4_S4_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_5_S5_R1_001.fastq.gz AMD_2002_5_S5_R2_001.fastq.gz ./Trimmed/AMD_2002_5_S5_forward_paired.fastq.gz ./Trimmed/AMD_2002_5_S5_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_5_S5_reverse_paired.fastq.gz ./Trimmed/AMD_2002_5_S5_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_6_S6_R1_001.fastq.gz AMD_2002_6_S6_R2_001.fastq.gz ./Trimmed/AMD_2002_6_S6_forward_paired.fastq.gz ./Trimmed/AMD_2002_6_S6_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_6_S6_reverse_paired.fastq.gz ./Trimmed/AMD_2002_6_S6_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_7_S7_R1_001.fastq.gz AMD_2002_7_S7_R2_001.fastq.gz ./Trimmed/AMD_2002_7_S7_forward_paired.fastq.gz ./Trimmed/AMD_2002_7_S7_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_7_S7_reverse_paired.fastq.gz ./Trimmed/AMD_2002_7_S7_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_8_S8_R1_001.fastq.gz AMD_2002_8_S8_R2_001.fastq.gz ./Trimmed/AMD_2002_8_S8_forward_paired.fastq.gz ./Trimmed/AMD_2002_8_S8_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_8_S8_reverse_paired.fastq.gz ./Trimmed/AMD_2002_8_S8_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_9_S9_R1_001.fastq.gz AMD_2002_9_S9_R2_001.fastq.gz ./Trimmed/AMD_2002_9_S9_forward_paired.fastq.gz ./Trimmed/AMD_2002_9_S9_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_9_S9_reverse_paired.fastq.gz ./Trimmed/AMD_2002_9_S9_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_10_S10_R1_001.fastq.gz AMD_2002_10_S10_R2_001.fastq.gz ./Trimmed/AMD_2002_10_S10_forward_paired.fastq.gz ./Trimmed/AMD_2002_10_S10_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_10_S10_reverse_paired.fastq.gz ./Trimmed/AMD_2002_10_S10_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_11_S11_R1_001.fastq.gz AMD_2002_11_S11_R2_001.fastq.gz ./Trimmed/AMD_2002_11_S11_forward_paired.fastq.gz ./Trimmed/AMD_2002_11_S11_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_11_S11_reverse_paired.fastq.gz ./Trimmed/AMD_2002_11_S11_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_12_S12_R1_001.fastq.gz AMD_2002_12_S12_R2_001.fastq.gz ./Trimmed/AMD_2002_12_S12_forward_paired.fastq.gz ./Trimmed/AMD_2002_12_S12_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_12_S12_reverse_paired.fastq.gz ./Trimmed/AMD_2002_12_S12_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2002_13_S13_R1_001.fastq.gz AMD_2002_13_S13_R2_001.fastq.gz ./Trimmed/AMD_2002_13_S13_forward_paired.fastq.gz ./Trimmed/AMD_2002_13_S13_forward_unpaired.fastq.gz ./Trimmed/AMD_2002_13_S13_reverse_paired.fastq.gz ./Trimmed/AMD_2002_13_S13_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

#2012
java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2012_1_S17_R1_001.fastq.gz 2012_1_S17_R2_001.fastq.gz ./Trimmed/2012_1_S17_forward_paired.fastq.gz ./Trimmed/2012_1_S17_forward_unpaired.fastq.gz ./Trimmed/2012_1_S17_reverse_paired.fastq.gz ./Trimmed/2012_1_S17_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2012_2_S18_R1_001.fastq.gz 2012_2_S18_R2_001.fastq.gz ./Trimmed/2012_2_S18_forward_paired.fastq.gz ./Trimmed/2012_2_S18_forward_unpaired.fastq.gz ./Trimmed/2012_2_S18_reverse_paired.fastq.gz ./Trimmed/2012_2_S18_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2012_3_S19_R1_001.fastq.gz 2012_3_S19_R2_001.fastq.gz ./Trimmed/2012_3_S19_forward_paired.fastq.gz ./Trimmed/2012_3_S19_forward_unpaired.fastq.gz ./Trimmed/2012_3_S19_reverse_paired.fastq.gz ./Trimmed/2012_3_S19_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2012_4_S20_R1_001.fastq.gz 2012_4_S20_R2_001.fastq.gz ./Trimmed/2012_4_S20_forward_paired.fastq.gz ./Trimmed/2012_4_S20_forward_unpaired.fastq.gz ./Trimmed/2012_4_S20_reverse_paired.fastq.gz ./Trimmed/2012_4_S20_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2012_5_S21_R1_001.fastq.gz 2012_5_S21_R2_001.fastq.gz ./Trimmed/2012_5_S21_forward_paired.fastq.gz ./Trimmed/2012_5_S21_forward_unpaired.fastq.gz ./Trimmed/2012_5_S21_reverse_paired.fastq.gz ./Trimmed/2012_5_S21_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2012_6_S22_R1_001.fastq.gz 2012_6_S22_R2_001.fastq.gz ./Trimmed/2012_6_S22_forward_paired.fastq.gz ./Trimmed/2012_6_S22_forward_unpaired.fastq.gz ./Trimmed/2012_6_S22_reverse_paired.fastq.gz ./Trimmed/2012_6_S22_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2012_7_S23_R1_001.fastq.gz 2012_7_S23_R2_001.fastq.gz ./Trimmed/2012_7_S23_forward_paired.fastq.gz ./Trimmed/2012_7_S23_forward_unpaired.fastq.gz ./Trimmed/2012_7_S23_reverse_paired.fastq.gz ./Trimmed/2012_7_S23_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2012_10_S24_R1_001.fastq.gz 2012_10_S24_R2_001.fastq.gz ./Trimmed/2012_10_S24_forward_paired.fastq.gz ./Trimmed/2012_10_S24_forward_unpaired.fastq.gz ./Trimmed/2012_10_S24_reverse_paired.fastq.gz ./Trimmed/2012_10_S24_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2012_9_S14_R1_001.fastq.gz AMD_2012_9_S14_R2_001.fastq.gz ./Trimmed/AMD_2012_9_S14_forward_paired.fastq.gz ./Trimmed/AMD_2012_9_S14_forward_unpaired.fastq.gz ./Trimmed/AMD_2012_9_S14_reverse_paired.fastq.gz ./Trimmed/AMD_2012_9_S14_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2012_10_S15_R1_001.fastq.gz AMD_2012_10_S15_R2_001.fastq.gz ./Trimmed/AMD_2012_10_S15_forward_paired.fastq.gz ./Trimmed/AMD_2012_10_S15_forward_unpaired.fastq.gz ./Trimmed/AMD_2012_10_S15_reverse_paired.fastq.gz ./Trimmed/AMD_2012_10_S15_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2012_11_S16_R1_001.fastq.gz AMD_2012_11_S16_R2_001.fastq.gz ./Trimmed/AMD_2012_11_S16_forward_paired.fastq.gz ./Trimmed/AMD_2012_11_S16_forward_unpaired.fastq.gz ./Trimmed/AMD_2012_11_S16_reverse_paired.fastq.gz ./Trimmed/AMD_2012_11_S16_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2012_12_S17_R1_001.fastq.gz AMD_2012_12_S17_R2_001.fastq.gz ./Trimmed/AMD_2012_12_S17_forward_paired.fastq.gz ./Trimmed/AMD_2012_12_S17_forward_unpaired.fastq.gz ./Trimmed/AMD_2012_12_S17_reverse_paired.fastq.gz ./Trimmed/AMD_2012_12_S17_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

##2017
java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2017_1_S25_R1_001.fastq.gz 2017_1_S25_R2_001.fastq.gz ./Trimmed/2017_1_S25_forward_paired.fastq.gz ./Trimmed/2017_1_S25_forward_unpaired.fastq.gz ./Trimmed/2017_1_S25_reverse_paired.fastq.gz ./Trimmed/2017_1_S25_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2017_3_S26_R1_001.fastq.gz 2017_3_S26_R2_001.fastq.gz ./Trimmed/2017_3_S26_forward_paired.fastq.gz ./Trimmed/2017_3_S26_forward_unpaired.fastq.gz ./Trimmed/2017_3_S26_reverse_paired.fastq.gz ./Trimmed/2017_3_S26_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2017_5_S27_R1_001.fastq.gz 2017_5_S27_R2_001.fastq.gz ./Trimmed/2017_5_S27_forward_paired.fastq.gz ./Trimmed/2017_5_S27_forward_unpaired.fastq.gz ./Trimmed/2017_5_S27_reverse_paired.fastq.gz ./Trimmed/2017_5_S27_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2017_6_S28_R1_001.fastq.gz 2017_6_S28_R2_001.fastq.gz ./Trimmed/2017_6_S28_forward_paired.fastq.gz ./Trimmed/2017_6_S28_forward_unpaired.fastq.gz ./Trimmed/2017_6_S28_reverse_paired.fastq.gz ./Trimmed/2017_6_S28_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2017_7_S29_R1_001.fastq.gz 2017_7_S29_R2_001.fastq.gz ./Trimmed/2017_7_S29_forward_paired.fastq.gz ./Trimmed/2017_7_S29_forward_unpaired.fastq.gz ./Trimmed/2017_7_S29_reverse_paired.fastq.gz ./Trimmed/2017_7_S29_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2017_8_S30_R1_001.fastq.gz 2017_8_S30_R2_001.fastq.gz ./Trimmed/2017_8_S30_forward_paired.fastq.gz ./Trimmed/2017_8_S30_forward_unpaired.fastq.gz ./Trimmed/2017_8_S30_reverse_paired.fastq.gz ./Trimmed/2017_8_S30_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 2017_9_S31_R1_001.fastq.gz 2017_9_S31_R2_001.fastq.gz ./Trimmed/2017_9_S31_forward_paired.fastq.gz ./Trimmed/2017_9_S31_forward_unpaired.fastq.gz ./Trimmed/2017_9_S31_reverse_paired.fastq.gz ./Trimmed/2017_9_S31_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2017_9_S19_R1_001.fastq.gz AMD_2017_9_S19_R2_001.fastq.gz ./Trimmed/AMD_2017_9_S19_forward_paired.fastq.gz ./Trimmed/AMD_2017_9_S19_forward_unpaired.fastq.gz ./Trimmed/AMD_2017_9_S19_reverse_paired.fastq.gz ./Trimmed/AMD_2017_9_S19_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2017_10_S20_R1_001.fastq.gz AMD_2017_10_S20_R2_001.fastq.gz ./Trimmed/AMD_2017_10_S20_forward_paired.fastq.gz ./Trimmed/AMD_2017_10_S20_forward_unpaired.fastq.gz ./Trimmed/AMD_2017_10_S20_reverse_paired.fastq.gz ./Trimmed/AMD_2017_10_S20_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2017_11_S21_R1_001.fastq.gz AMD_2017_11_S21_R2_001.fastq.gz ./Trimmed/AMD_2017_11_S21_forward_paired.fastq.gz ./Trimmed/AMD_2017_11_S21_forward_unpaired.fastq.gz ./Trimmed/AMD_2017_11_S21_reverse_paired.fastq.gz ./Trimmed/AMD_2017_11_S21_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50

java -jar /Volumes/data_a3/Hzea/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar PE -phred33 -threads 8 AMD_2017_12_S22_R1_001.fastq.gz AMD_2017_12_S22_R2_001.fastq.gz ./Trimmed/AMD_2017_12_S22_forward_paired.fastq.gz ./Trimmed/AMD_2017_12_S22_forward_unpaired.fastq.gz ./Trimmed/AMD_2017_12_S22_reverse_paired.fastq.gz ./Trimmed/AMD_2017_12_S22_reverse_unpaired.fastq.gz SLIDINGWINDOW:5:20 MINLEN:50
