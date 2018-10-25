#script to merge PE reads
#MF 10/25/2018

#start in this directory: cd /media/megan/"New Volume"/Raw_Sequence_Data/Hz_FieldColl_ddRAD_Libraries

#AMD lib8 files first

cd AMD_lib8_raw
mkdir FLASH_output_AMD_lib8

~/src/FLASH-1.2.11-Linux-x86_64/flash pool1_S1_L001_R1_001.fastq.gz pool1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib8_pool1 >&1 | tee AMD_lib8_pool1_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash pool3_S2_L001_R1_001.fastq.gz pool3_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib8_pool3 >&1 | tee AMD_lib8_pool3_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash pool4_S3_L001_R1_001.fastq.gz pool4_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib8_pool4 >&1 | tee AMD_lib8_pool4_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash pool5_S4_L001_R1_001.fastq.gz pool5_S4_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib8_pool5 >&1 | tee AMD_lib8_pool5_flash.log

mv AMD_lib8* ./FLASH_output_AMD_lib8

###New File AMD lib7

cd ../AMD_lib7_raw
mkdir FLASH_output_AMD_lib7

~/src/FLASH-1.2.11-Linux-x86_64/flash AMD-7-1_S1_L001_R1_001.fastq.gz AMD-7-1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib7_pool1 >&1 | tee AMD_lib7_pool1_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash AMD-7-2_S2_L001_R1_001.fastq.gz AMD-7-2_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib7_pool1 >&1 | tee AMD_lib7_pool2_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash AMD-7-3_S3_L001_R1_001.fastq.gz AMD-7-3_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib7_pool3 >&1 | tee AMD_lib7_pool3_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash AMD-7-4_S4_L001_R1_001.fastq.gz AMD-7-4_S4_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib7_pool4 >&1 | tee AMD_lib7_pool4_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash AMD-7-5_S5_L001_R1_001.fastq.gz AMD-7-5_S5_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib7_pool5 >&1 | tee AMD_lib7_pool5_flash.log

mv AMD_lib7* ./FLASH_output_AMD_lib7

###New File AMD lib 6

cd ../AMD_lib6_raw
mkdir FLASH_output_AMD_lib6

~/src/FLASH-1.2.11-Linux-x86_64/flash AMD-6-1_S1_L001_R1_001.fastq.gz AMD-6-1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib6_pool1 >&1 | tee AMD_lib6_pool1_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash AMD-6-2_S2_L001_R1_001.fastq.gz AMD-6-2_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib6_pool1 >&1 | tee AMD_lib6_pool2_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash AMD-6-3_S3_L001_R1_001.fastq.gz AMD-6-3_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib6_pool3 >&1 | tee AMD_lib6_pool3_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash AMD-6-4_S4_L001_R1_001.fastq.gz AMD-6-4_S4_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib6_pool4 >&1 | tee AMD_lib6_pool4_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash AMD-6-5_S5_L001_R1_001.fastq.gz AMD-6-5_S5_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib6_pool5 >&1 | tee AMD_lib6_pool5_flash.log

mv AMD_lib6* ./FLASH_output_AMD_lib6


