#sam to bam conversion WGRS project, filtering and removed duplicates
#Counting reads in .bam files on 01/14/2020
#MF 

#used samtools 0.1.19 

cd /media/megan/"New Volume"/Hzea_WGRS_Bowtie2_output/Alignment_by_years_01142020

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


#filtering bam files based on mapping quality

for sample in *.bam
do
	echo $sample
	describer=$(echo ${sample} | sed 's/.bam//')
	echo $describer
	
	#filtering based on quality score
	samtools view -bq 5 $sample -o ${describer}.uns.bam

	#Sort BAM file
	samtools sort ${describer}.uns.bam ${describer}_filtered

	#Index BAM file
	samtools index ${describer}_filtered.bam

	#Remove intermediate files
	rm ${describer}.uns.bam

done

#marking and removing duplicates
for sample in *_filtered.bam
do
	echo $sample
	describer=$(echo ${sample} | sed 's/_filtered.bam//')
	echo $describer
	
	#convert file from SAM to BAM format
	/home/megan/src/gatk-4.0.12.0/gatk MarkDuplicates -I $sample -M ${describer}.metrics --REMOVE_DUPLICATES true -O ${describer}_marked.uns.bam

	#Sort filtered BAM file
	samtools sort ${describer}_marked.uns.bam ${describer}_marked

	#Index BAM file
	samtools index ${describer}_marked.bam

	#Remove intermediate files
	rm ${describer}_marked.uns.bam

done

#sanity check - no dups left
/home/megan/src/gatk-4.0.12.0/gatk MarkDuplicates -I HZ_2002_01_marked.bam -M HZ_2002_01_marked.metrics --REMOVE_DUPLICATES false -O HZ_2002_01_marked.uns.bam
rm HZ_2002_01_marked.uns.bam


#getting a list of bam files
ls *marked.bam > WGRS_Hzea_BamFiles.txt


