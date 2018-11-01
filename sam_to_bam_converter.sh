#sam to bam conversion
#MF 10/31/2018

cd /media/megan/"New Volume"/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments

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

#getting a list of sam files
ls *.bam > Hzea_BamFiles.txt

