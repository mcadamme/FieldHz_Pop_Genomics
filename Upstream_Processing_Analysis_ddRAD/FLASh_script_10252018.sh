#script to merge PE reads
#MF 10/25/2018

#start in this directory: cd /media/megan/"New Volume"/Raw_Sequence_Data/Hz_FieldColl_ddRAD_Libraries

#AMD lib8 files first
###NOTE: Pool 2 from lib 8 failed, Pool 1 does not contain much data.

cd AMD_lib8_raw
mkdir FLASH_output_AMD_lib8

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_pool1_S1_L001_R1_001.fastq.gz trimmed_pool1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib8_pool1 >&1 | tee AMD_lib8_pool1_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_pool3_S2_L001_R1_001.fastq.gz trimmed_pool3_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib8_pool3 >&1 | tee AMD_lib8_pool3_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_pool4_S3_L001_R1_001.fastq.gz trimmed_pool4_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib8_pool4 >&1 | tee AMD_lib8_pool4_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_pool5_S4_L001_R1_001.fastq.gz trimmed_pool5_S4_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib8_pool5 >&1 | tee AMD_lib8_pool5_flash.log

mv AMD_lib8* ./FLASH_output_AMD_lib8

###New File AMD lib7
###Lib7 pool 4 was also low yielding

cd ../AMD_lib7_raw
mkdir FLASH_output_AMD_lib7

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_AMD-7-1_S1_L001_R1_001.fastq.gz trimmed_AMD-7-1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib7_pool1 >&1 | tee AMD_lib7_pool1_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_AMD-7-2_S2_L001_R1_001.fastq.gz trimmed_AMD-7-2_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib7_pool2 >&1 | tee AMD_lib7_pool2_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_AMD-7-3_S3_L001_R1_001.fastq.gz trimmed_AMD-7-3_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib7_pool3 >&1 | tee AMD_lib7_pool3_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_AMD-7-4_S4_L001_R1_001.fastq.gz trimmed_AMD-7-4_S4_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib7_pool4 >&1 | tee AMD_lib7_pool4_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_AMD-7-5_S5_L001_R1_001.fastq.gz trimmed_AMD-7-5_S5_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib7_pool5 >&1 | tee AMD_lib7_pool5_flash.log

mv AMD_lib7* ./FLASH_output_AMD_lib7

###New File AMD lib 6

cd ../AMD_lib6_raw
mkdir FLASH_output_AMD_lib6

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_AMD-6-1_S1_L001_R1_001.fastq.gz trimmed_AMD-6-1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib6_pool1 >&1 | tee AMD_lib6_pool1_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_AMD-6-2_S2_L001_R1_001.fastq.gz trimmed_AMD-6-2_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib6_pool2 >&1 | tee AMD_lib6_pool2_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_AMD-6-3_S3_L001_R1_001.fastq.gz trimmed_AMD-6-3_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib6_pool3 >&1 | tee AMD_lib6_pool3_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_AMD-6-4_S4_L001_R1_001.fastq.gz trimmed_AMD-6-4_S4_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib6_pool4 >&1 | tee AMD_lib6_pool4_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_AMD-6-5_S5_L001_R1_001.fastq.gz trimmed_AMD-6-5_S5_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib6_pool5 >&1 | tee AMD_lib6_pool5_flash.log

mv AMD_lib6* ./FLASH_output_AMD_lib6


###New file AMD_lib5
#note Lib5 index2 had low read counts

cd ../AMD_lib5_raw
mkdir FLASH_output_AMD_lib5

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_6_S4_L001_R1_001.fastq.gz trimmed_6_S4_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib5_idx6 >&1 | tee AMD_lib5_idx6_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_4_S3_L001_R1_001.fastq.gz trimmed_4_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib5_idx4 >&1 | tee AMD_lib5_idx4_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_2_S2_L001_R1_001.fastq.gz trimmed_2_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib5_idx2 >&1 | tee AMD_lib5_idx2_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_1_S1_L001_R1_001.fastq.gz trimmed_1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib5_idx1 >&1 | tee AMD_lib5_idx1_flash.log


mv AMD_lib5* ./FLASH_output_AMD_lib5


###New file AMD_lib4

cd ../AMD_lib4_raw
mkdir FLASH_output_AMD_lib4

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_Idx1_S1_L001_R1_001.fastq.gz trimmed_Idx1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib4_idx1 >&1 | tee AMD_lib4_idx1_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_Idx2_S2_L001_R1_001.fastq.gz trimmed_Idx2_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib4_idx2 >&1 | tee AMD_lib4_idx2_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_Idx4_S3_L001_R1_001.fastq.gz trimmed_Idx4_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib4_idx4 >&1 | tee AMD_lib4_idx4_flash.log


mv AMD_lib4* ./FLASH_output_AMD_lib4


###New file AMD_lib3

cd ../AMD_lib3_raw
mkdir FLASH_output_AMD_lib3

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_Idx1_S1_L001_R1_001.fastq.gz trimmed_Idx1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib3_idx1 >&1 | tee AMD_lib3_idx1_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_Idx2_S2_L001_R1_001.fastq.gz trimmed_Idx2_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib3_idx2 >&1 | tee AMD_lib3_idx2_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_Idx4_S3_L001_R1_001.fastq.gz trimmed_Idx4_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib3_idx4 >&1 | tee AMD_lib3_idx4_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_Idx6_S4_L001_R1_001.fastq.gz trimmed_Idx6_S4_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib3_idx6 >&1 | tee AMD_lib3_idx6_flash.log

mv AMD_lib3* ./FLASH_output_AMD_lib3

###New file AMD_lib2

cd ../AMD_lib2_raw
mkdir FLASH_output_AMD_lib2

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_index1_S1_L001_R1_001.fastq.gz trimmed_index1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib2_idx1 >&1 | tee AMD_lib2_idx1_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_index2_S2_L001_R1_001.fastq.gz trimmed_index2_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib2_idx2 >&1 | tee AMD_lib2_idx2_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_index4_S3_L001_R1_001.fastq.gz trimmed_index4_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib2_idx4 >&1 | tee AMD_lib2_idx4_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_index6_S4_L001_R1_001.fastq.gz trimmed_index6_S4_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib2_idx6 >&1 | tee AMD_lib2_idx6_flash.log

mv AMD_lib2* ./FLASH_output_AMD_lib2

###New file AMD_lib1

cd ../AMD_lib1_raw
mkdir FLASH_output_AMD_lib1

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_index1_S1_L001_R1_001.fastq.gz trimmed_index1_S1_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib1_idx1 >&1 | tee AMD_lib1_idx1_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_index2_S2_L001_R1_001.fastq.gz trimmed_index2_S2_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib1_idx2 >&1 | tee AMD_lib1_idx2_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_index4_S3_L001_R1_001.fastq.gz trimmed_index4_S3_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib1_idx4 >&1 | tee AMD_lib1_idx4_flash.log

~/src/FLASH-1.2.11-Linux-x86_64/flash trimmed_index12_S4_L001_R1_001.fastq.gz trimmed_index12_S4_L001_R2_001.fastq.gz --max-overlap=200 -o AMD_lib1_idx12 >&1 | tee AMD_lib1_idx12_flash.log

mv AMD_lib1* ./FLASH_output_AMD_lib1


###Copied all Flash output folders to new directory (using gui) for demultiplexing.

