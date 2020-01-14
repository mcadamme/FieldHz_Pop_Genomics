#sam to bam conversion WGRS project, filtering and removed duplicates
#Counting reads in .bam files on 01/14/2020
#MF 

#used samtools 0.1.19 

cd /media/megan/"New Volume"/Hzea_WGRS_Bowtie2_output

for sample in *.sam
do
	echo $sample
	describer=$(echo ${sample} | sed 's/.sam//')
	echo $describer
	
	#convert file from SAM to BAM format
	samtools view -bS $sample -o ${describer}.uns.bam

	#Sort BAM file
	samtools sort ${describer}.uns.bam ${describer}

	#Index BAM file
	samtools index ${describer}.bam

	#Remove intermediate files
	rm ${describer}.uns.bam

done


#filtering bam files based on mapping quality (applies to v.2 files)

for sample in *_paired.bam
do
	echo $sample
	describer=$(echo ${sample} | sed 's/_paired.bam//')
	echo $describer
	
	#filtering based on quality score
	samtools view -bq 5 $sample -o ${describer}.uns.bam

	#Sort BAM file
	samtools sort ${describer}.uns.bam ${describer}_filtered.paired

	#Index BAM file
	samtools index ${describer}_filtered.paired.bam

	#Remove intermediate files
	rm ${describer}.uns.bam

done

#marking and removing duplicates (applies to v.2 files)
for sample in *_filtered.paired.bam
do
	echo $sample
	describer=$(echo ${sample} | sed 's/_filtered.paired.bam//')
	echo $describer
	
	#convert file from SAM to BAM format
	/home/megan/src/gatk-4.0.12.0/gatk MarkDuplicates -I $sample -M ${describer}.metrics --REMOVE_DUPLICATES true -O ${describer}.marked.uns.bam

	#Sort filtered BAM file
	samtools sort ${describer}.marked.uns.bam ${describer}.marked

	#Index BAM file
	samtools index ${describer}.marked.bam

	#Remove intermediate files
	rm ${describer}.marked.uns.bam

done

#getting a list of bam files
ls *marked.bam > WGRS_Hzea_BamFiles.txt

#getting primary aligned reads for each file.
samtools view -c -F 260 *marked.bam
