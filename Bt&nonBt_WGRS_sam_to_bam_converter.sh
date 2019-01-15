#sam to bam conversion Bt and Non-Bt WGRS project 1/14/2019
#did some filtering and removed duplicates on 1/14/2019
#MF 

#used samtools 0.1.19 

mkdir /media/megan/"New Volume"/Hzea_WGRS_Bowtie2_output/BtandNonBt_alignmentFiles

cp /media/megan/"New Volume"/Hzea_WGRS_Bowtie2_output/*Bt*_paired.sam /media/megan/"New Volume"/Hzea_WGRS_Bowtie2_output/BtandNonBt_alignmentFiles

cd /media/megan/"New Volume"/Hzea_WGRS_Bowtie2_output/BtandNonBt_alignmentFiles

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

#marking and removing duplicates
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

#getting a list of sam files
ls *marked.bam > WGRS_Hzea_BtandNonBt_BamFiles.txt
