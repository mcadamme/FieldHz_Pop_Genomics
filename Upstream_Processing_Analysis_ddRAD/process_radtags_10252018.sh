#script to demultiplex with process_radtags
#10/25/18
#used stacks v. 2.2

#Note, started in /media/megan/"New Volume"/Hz_PopGen_ddRAD_demult/AMD_libs1thru8_extended_frags

#AMDlib8 - dropping pool 1 b/c of low quality
mkdir ./FLASH_output_AMD_lib8/demult8_3 ./FLASH_output_AMD_lib8/demult8_4 ./FLASH_output_AMD_lib8/demult8_5

#pool 3
process_radtags -f ./FLASH_output_AMD_lib8/AMD_lib8_pool3.extendedFrags.fastq -b ./barcodes/barcode_lib78 -o ./FLASH_output_AMD_lib8/demult8_3 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#pool 4
process_radtags -f ./FLASH_output_AMD_lib8/AMD_lib8_pool4.extendedFrags.fastq -b ./barcodes/barcode_lib78 -o ./FLASH_output_AMD_lib8/demult8_4 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#pool 5
process_radtags -f ./FLASH_output_AMD_lib8/AMD_lib8_pool5.extendedFrags.fastq -b ./barcodes/barcode_lib78 -o ./FLASH_output_AMD_lib8/demult8_5 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#AMDlib7 - dropping pool 4 b/c of low quality
mkdir ./FLASH_output_AMD_lib7/demult7_1 ./FLASH_output_AMD_lib7/demult7_2 ./FLASH_output_AMD_lib7/demult7_3 ./FLASH_output_AMD_lib7/demult7_5

#pool 1
process_radtags -f ./FLASH_output_AMD_lib7/AMD_lib7_pool1.extendedFrags.fastq -b ./barcodes/barcode_lib78 -o ./FLASH_output_AMD_lib7/demult7_1 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#pool 2
process_radtags -f ./FLASH_output_AMD_lib7/AMD_lib7_pool2.extendedFrags.fastq -b ./barcodes/barcode_lib78 -o ./FLASH_output_AMD_lib7/demult7_2 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#pool 3
process_radtags -f ./FLASH_output_AMD_lib7/AMD_lib7_pool3.extendedFrags.fastq -b ./barcodes/barcode_lib78 -o ./FLASH_output_AMD_lib7/demult7_3 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#pool 5
process_radtags -f ./FLASH_output_AMD_lib7/AMD_lib7_pool5.extendedFrags.fastq -b ./barcodes/barcode_lib78 -o ./FLASH_output_AMD_lib7/demult7_5 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#AMDlib6 
mkdir ./FLASH_output_AMD_lib6/demult6_1 ./FLASH_output_AMD_lib6/demult6_2 ./FLASH_output_AMD_lib6/demult6_3 ./FLASH_output_AMD_lib6/demult6_4 ./FLASH_output_AMD_lib6/demult6_5

#pool 1
process_radtags -f ./FLASH_output_AMD_lib6/AMD_lib6_pool1.extendedFrags.fastq -b ./barcodes/barcode_lib6 -o ./FLASH_output_AMD_lib6/demult6_1 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#pool 2
process_radtags -f ./FLASH_output_AMD_lib6/AMD_lib6_pool2.extendedFrags.fastq -b ./barcodes/barcode_lib6 -o ./FLASH_output_AMD_lib6/demult6_2 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#pool 3
process_radtags -f ./FLASH_output_AMD_lib6/AMD_lib6_pool3.extendedFrags.fastq -b ./barcodes/barcode_lib6 -o ./FLASH_output_AMD_lib6/demult6_3 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#pool 4
process_radtags -f ./FLASH_output_AMD_lib6/AMD_lib6_pool4.extendedFrags.fastq -b ./barcodes/barcode_lib6 -o ./FLASH_output_AMD_lib6/demult6_4 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#pool 5
process_radtags -f ./FLASH_output_AMD_lib6/AMD_lib6_pool5.extendedFrags.fastq -b ./barcodes/barcode_lib6 -o ./FLASH_output_AMD_lib6/demult6_5 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#AMDlib5 - idx 2 had low read counts, so eliminating.
mkdir ./FLASH_output_AMD_lib5/demult5_1 ./FLASH_output_AMD_lib5/demult5_4 ./FLASH_output_AMD_lib5/demult5_6 

#idx 1
process_radtags -f ./FLASH_output_AMD_lib5/AMD_lib5_idx1.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib5/demult5_1 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 4
process_radtags -f ./FLASH_output_AMD_lib5/AMD_lib5_idx4.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib5/demult5_4 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 6
process_radtags -f ./FLASH_output_AMD_lib5/AMD_lib5_idx6.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib5/demult5_6 -e 'ecoRI' -qD -y 'gzfastq' --inline_null


#AMDlib4 
mkdir ./FLASH_output_AMD_lib4/demult4_1 ./FLASH_output_AMD_lib4/demult4_2 ./FLASH_output_AMD_lib4/demult4_4 

#idx 1
process_radtags -f ./FLASH_output_AMD_lib4/AMD_lib4_idx1.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib4/demult4_1 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 2
process_radtags -f ./FLASH_output_AMD_lib4/AMD_lib4_idx2.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib4/demult4_2 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 4
process_radtags -f ./FLASH_output_AMD_lib4/AMD_lib4_idx4.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib4/demult4_4 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#AMDlib3 
mkdir ./FLASH_output_AMD_lib3/demult3_1 ./FLASH_output_AMD_lib3/demult3_2 ./FLASH_output_AMD_lib3/demult3_4 ./FLASH_output_AMD_lib3/demult3_6 

#idx 1
process_radtags -f ./FLASH_output_AMD_lib3/AMD_lib3_idx1.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib3/demult3_1 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 2
process_radtags -f ./FLASH_output_AMD_lib3/AMD_lib3_idx2.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib3/demult3_2 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 4
process_radtags -f ./FLASH_output_AMD_lib3/AMD_lib3_idx4.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib3/demult3_4 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 6
process_radtags -f ./FLASH_output_AMD_lib3/AMD_lib3_idx6.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib3/demult3_6 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#AMDlib2 
mkdir ./FLASH_output_AMD_lib2/demult2_1 ./FLASH_output_AMD_lib2/demult2_2 ./FLASH_output_AMD_lib2/demult2_4 ./FLASH_output_AMD_lib2/demult2_6 

#idx 1
process_radtags -f ./FLASH_output_AMD_lib2/AMD_lib2_idx1.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib2/demult2_1 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 2
process_radtags -f ./FLASH_output_AMD_lib2/AMD_lib2_idx2.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib2/demult2_2 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 4
process_radtags -f ./FLASH_output_AMD_lib2/AMD_lib2_idx4.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib2/demult2_4 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 6
process_radtags -f ./FLASH_output_AMD_lib2/AMD_lib2_idx6.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib2/demult2_6 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#AMDlib1 
mkdir ./FLASH_output_AMD_lib1/demult1_1 ./FLASH_output_AMD_lib1/demult1_2 ./FLASH_output_AMD_lib1/demult1_4 ./FLASH_output_AMD_lib1/demult1_12 

#idx 1
process_radtags -f ./FLASH_output_AMD_lib1/AMD_lib1_idx1.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib1/demult1_1 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 2
process_radtags -f ./FLASH_output_AMD_lib1/AMD_lib1_idx2.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib1/demult1_2 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 4
process_radtags -f ./FLASH_output_AMD_lib1/AMD_lib1_idx4.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib1/demult1_4 -e 'ecoRI' -qD -y 'gzfastq' --inline_null

#idx 12
process_radtags -f ./FLASH_output_AMD_lib1/AMD_lib1_idx12.extendedFrags.fastq -b ./barcodes/barcode_libs1thru5 -o ./FLASH_output_AMD_lib1/demult1_12 -e 'ecoRI' -qD -y 'gzfastq' --inline_null
