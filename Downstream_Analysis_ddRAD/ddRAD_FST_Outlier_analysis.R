#this is my script to analyze my H. zea ddRADseq polymorphism data. 
#M.Fritz 11/23/18; updated 3/23/2021

#loading libraries
library(LEA);library(OutFLANK);library(adegenet);library(vcfR)


##### Loading scafname and gff files for outlier ID #####

#loading identifier file
scaf_names <- read.table("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS/scaffold_and_contig_names.txt", header = T)

#Loading gff3 file 
gff3 <- read.table("~/Genomes/HzOGS2-15205-fixed_note-added.gff3", sep="\t", stringsAsFactors = F)
head(gff3)


#Outlier Analysis of BCFtools output
setwd("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS/ddRAD")

#converting to the proper file format
#overall pop across years
vcf2geno(input.file = "thinned_FieldHzea_allpops.recode.vcf", output.file = "FieldHzea_allpops.geno")

Field_Hzea_allpops.geno_in <- read.fwf("FieldHzea_allpops.geno", width=rep(1,259))
Field_Hzea_allpops.geno <- t(Field_Hzea_allpops.geno_in)
Field_Hzea_allpops_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/allpops_no_low.txt")

colnames(Field_Hzea_allpops_meta) <- c("samp_names")
Field_Hzea_allpops_meta$pop <- regmatches(Field_Hzea_allpops_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_allpops_meta$samp_names))

#OutFLANK analysis for all years
OF_allpops_SNPs <- MakeDiploidFSTMat(Field_Hzea_allpops.geno, locusNames=seq(1, 14398, by=1), popNames=Field_Hzea_allpops_meta$pop) # Provide the transposed genotype file, locus names, and population names.
OF_out <- OutFLANK(FstDataFrame=OF_allpops_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=4, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps across years")


#getting positions of the outliers
vcf_allpops <- read.vcfR("thinned_FieldHzea_allpops.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf_allpops))

outliers_allyears_df <- data.frame(OF_out$results, header = T)
outliers_allyears_df_withpos <- cbind(vcfann[,c(1:2)],outliers_allyears_df)
outliers_allyears <- subset(outliers_allyears_df_withpos, OutlierFlag == "TRUE")
outliers_allyears$ChromPos <- paste0(outliers_allyears$CHROM,"_",outliers_allyears$POS)

outliers_allyears_uniqScafs <- unique(outliers_allyears$CHROM)
outliers_allyears_uniqScafs

write.table(outliers_allyears, file = "FileS1_outliers_allyears.txt", col.names = T, row.names = F)

min(outliers_allyears$FST)
max(outliers_allyears$FST)

#all_years qval 0.05
png("Fig2_Outflank_Analysis_allyears_qval0.05.png", units = "px", height = 600, width = 800)
par(mar = c(5,6,4,2))
plot(OF_out$results$He, OF_out$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.3), cex.axis = 2, cex.lab = 2.2)
points(OF_out$results$He[OF_out$results$qvalues<0.05], OF_out$results$FST[OF_out$results$qvalues<0.05], pch=21, col="black")

dev.off()


##### Checking scaffolds with outliers against LGs in the super-scaffolded assembly #####

ord_Chroms_final <- read.csv("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS/ord_Chroms_final.csv", header = T)

outliers_allyears_withLG <- merge(outliers_allyears, ord_Chroms_final, by.x = "CHROM", by.y = "V6")
nrow(outliers_allyears_withLG) #42 of 53 outliers could be mapped to Linkage Groups or chromosomes in the super-scaffolded assembly

#gives unique LGs with outliers
outliers_uniqLGs <- unique(outliers_allyears_withLG$chr)
outliers_uniqLGs

#what percentage of these outliers were also in the top 5% of variance contributing SNPs from the DAPC?
DAPC_top5 <- read.table("/media/megan/New Volume1/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/mpileupANDvcftools_output/DAPC2_top5per.txt", header = F)
colnames(DAPC_top5) <- c("CHROM", "POS")
DAPC_top5$ChromPos <- paste(DAPC_top5$CHROM, "_", DAPC_top5$POS)

merged_chromPos <- merge(outliers_allyears, DAPC_top5, by = "ChromPos") #17 of 22 mapped to chromosomes

#checking these 17 outliers and their positions.
mergedOLs <- merge(merged_chromPos, ord_Chroms_final, by.x = "CHROM.x", by.y = "V6")
mergedoutliers_uniqLGs <- unique(mergedOLs$chr)
mergedoutliers_uniqLGs
