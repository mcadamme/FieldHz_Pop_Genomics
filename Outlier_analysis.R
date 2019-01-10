#this is my script to analyze my H. zea ddRADseq polymorphism data. 
#M.Fritz 11/23/18

library(LEA);library(OutFLANK);library(adegenet);library(vcfR)

setwd("/media/megan/New Volume/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/mpileupANDvcftools_output")

#converting to the proper file format
vcf2geno(input.file = "thinned_FieldHzea2002and2007.recode.vcf", output.file = "FieldHzea2002and2007.geno")
vcf2geno(input.file = "thinned_FieldHzea2007and2012.recode.vcf", output.file = "FieldHzea2007and2012.geno")
vcf2geno(input.file = "thinned_FieldHzea2012and2016.recode.vcf", output.file = "FieldHzea2012and2016.geno")
vcf2geno(input.file = "thinned_FieldHzea2002and2016.recode.vcf", output.file = "FieldHzea2002and2016.geno")


#loading dataset for 2002 and 2007
Field_Hzea_2002and2007.geno_in <- read.fwf("FieldHzea2002and2007.geno", width=rep(1,141))
Field_Hzea_2002and2007.geno <- t(Field_Hzea_2002and2007.geno_in)
Field_Hzea_2002and2007_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2007.txt")

colnames(Field_Hzea_2002and2007_meta) <- c("samp_names")
Field_Hzea_2002and2007_meta$pop <- regmatches(Field_Hzea_2002and2007_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_2002and2007_meta$samp_names))

#OutFLANK analysis for 2002 and 2007
OF_2002and2007_SNPs <- MakeDiploidFSTMat(Field_Hzea_2002and2007.geno, locusNames=seq(1, 27590, by=1), popNames=Field_Hzea_2002and2007_meta$pop) # Provide the transposed genotype file, locus names, and population names.
OF_out1 <- OutFLANK(FstDataFrame=OF_2002and2007_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out1, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2002 & 2007")

outliers1 <- which(OF_out1$results$OutlierFlag=="TRUE")
print(outliers1)

#getting positions of the outliers
vcf2002and2007 <- read.vcfR("thinned_FieldHzea2002and2007.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2002and2007))
outliers_2002and2007 <- vcfann[outliers1,]

#write.table(outliers_2002and2007, file = "outliers_2002and2007.txt", col.names = T, row.names = F)

nrow(outliers_2002and2007)
length(unique(outliers_2002and2007$CHROM))
outliers_2002and2007$MARK <- paste(outliers_2002and2007$CHROM, outliers_2002and2007$POS, sep = "_")


#loading dataset for 2007 and 2012
Field_Hzea_2007and2012.geno_in <- read.fwf("FieldHzea2007and2012.geno", width=rep(1,143))
Field_Hzea_2007and2012.geno <- t(Field_Hzea_2007and2012.geno_in)
Field_Hzea_2007and2012_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2007and2012.txt")

colnames(Field_Hzea_2007and2012_meta) <- c("samp_names")
Field_Hzea_2007and2012_meta$pop <- regmatches(Field_Hzea_2007and2012_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_2007and2012_meta$samp_names))

#OutFLANK analysis for 2007 and 2012
OF_2007and2012_SNPs <- MakeDiploidFSTMat(Field_Hzea_2007and2012.geno, locusNames=seq(1, 27590, by=1), popNames=Field_Hzea_2007and2012_meta$pop) 
OF_out2 <- OutFLANK(FstDataFrame=OF_2007and2012_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out2, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2007 & 2012")

outliers2 <- which(OF_out2$results$OutlierFlag=="TRUE")
print(outliers2)

vcf2007and2012 <- read.vcfR("thinned_FieldHzea2007and2012.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2007and2012))
outliers_2007and2012 <- vcfann[outliers2,]

#write.table(outliers_2007and2012, file = "outliers_2007and2012.txt", col.names = T, row.names = F)

nrow(outliers_2007and2012)
length(unique(outliers_2007and2012$CHROM))
outliers_2007and2012$MARK <- paste(outliers_2007and2012$CHROM, outliers_2007and2012$POS, sep = "_")



#2012 through 2016
Field_Hzea_2012and2016.geno_in <- read.fwf("FieldHzea2012and2016.geno", width=rep(1,118))
Field_Hzea_2012and2016.geno <- t(Field_Hzea_2012and2016.geno_in)
Field_Hzea_2012and2016_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2012and2016.txt")

colnames(Field_Hzea_2012and2016_meta) <- c("samp_names")
Field_Hzea_2012and2016_meta$pop <- regmatches(Field_Hzea_2012and2016_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_2012and2016_meta$samp_names))

#OutFLANK analysis for 2012 and 2016
OF_2012and2016_SNPs <- MakeDiploidFSTMat(Field_Hzea_2012and2016.geno, locusNames=seq(1, 27590, by=1), popNames=Field_Hzea_2012and2016_meta$pop) 
OF_out3 <- OutFLANK(FstDataFrame=OF_2012and2016_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out3, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2012 & 2016")

outliers3 <- which(OF_out3$results$OutlierFlag=="TRUE")
print(outliers3)

vcf2012and2016 <- read.vcfR("thinned_FieldHzea2012and2016.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2012and2016))
outliers_2012and2016 <- vcfann[outliers3,]

#write.table(outliers_2012and2016, file = "outliers_2012and2016.txt", col.names = T, row.names = F)

nrow(outliers_2012and2016)
length(unique(outliers_2012and2016$CHROM))
outliers_2012and2016$MARK <- paste(outliers_2012and2016$CHROM, outliers_2012and2016$POS, sep = "_")



#2002 through 2016
Field_Hzea_2002and2016.geno_in <- read.fwf("FieldHzea2002and2016.geno", width=rep(1,116))
Field_Hzea_2002and2016.geno <- t(Field_Hzea_2002and2016.geno_in)
Field_Hzea_2002and2016_meta <- read.table("/home/megan/scripts/Field_HZ_Pop_Genomics/pop_files/2002and2016.txt")

colnames(Field_Hzea_2002and2016_meta) <- c("samp_names")
Field_Hzea_2002and2016_meta$pop <- regmatches(Field_Hzea_2002and2016_meta$samp_names, regexpr("20[[:digit:]]+", Field_Hzea_2002and2016_meta$samp_names))

#OutFLANK analysis for 2002 and 2016
OF_2002and2016_SNPs <- MakeDiploidFSTMat(Field_Hzea_2002and2016.geno, locusNames=seq(1, 27590, by=1), popNames=Field_Hzea_2002and2016_meta$pop) 
OF_out4 <- OutFLANK(FstDataFrame=OF_2002and2016_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=2, qthreshold=0.05)
OutFLANKResultsPlotter(OF_out4, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom = TRUE, titletext="Scan for selective sweeps 2002 & 2016")

outliers4 <- which(OF_out4$results$OutlierFlag=="TRUE")
print(outliers4)

vcf2002and2016 <- read.vcfR("thinned_FieldHzea2002and2016.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf2002and2016))
outliers_2002and2016 <- vcfann[outliers4,]

#write.table(outliers_2002and2016, file = "outliers_2002and2016.txt", col.names = T, row.names = F)

nrow(outliers_2002and2016)
length(unique(outliers_2002and2016$CHROM))
outliers_2002and2016$MARK <- paste(outliers_2002and2016$CHROM, outliers_2002and2016$POS, sep = "_")

#Outlier Analysis Figures for Pub

#2002&2016 qthreshold 0.05
png("Outflank_Analysis_2002and2016_q0.05.png", units = "px", height = 600, width = 800)
par(mar = c(5,6,4,2))
plot(OF_out4$results$He, OF_out4$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.5), cex.axis = 2, cex.lab = 2.2)
points(OF_out4$results$He[OF_out4$results$qvalues<0.05], OF_out4$results$FST[OF_out4$results$qvalues<0.05], pch=21, col="black")

dev.off()

#qthreshold 0.05
png("Outflank_Analysis_multipanel_years_q0.05.png", units = "px", height = 1200, width = 800)

par(mfrow = c(3,1))
par(mar = c(5,6,4,2))

plot(OF_out1$results$He, OF_out1$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out1$results$He[OF_out1$results$qvalues<0.05], OF_out1$results$FST[OF_out1$results$qvalues<0.05], pch=21, col="black")
title("A", cex.main = 3, adj  = 0.05, line = -5)

plot(OF_out2$results$He, OF_out2$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out2$results$He[OF_out2$results$qvalues<0.05], OF_out2$results$FST[OF_out2$results$qvalues<0.05], pch=21, col="black")
title("B", cex.main = 3, adj  = 0.05, line = -5)

plot(OF_out3$results$He, OF_out3$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out3$results$He[OF_out3$results$qvalues<0.05], OF_out3$results$FST[OF_out3$results$qvalues<0.05], pch=21, col="black")
title("C", cex.main = 3, adj  = 0.05, line = -5)

dev.off()

#qthreshold 0.01
png("Outflank_Analysis_multipanel_years_q0.01.png", units = "px", height = 1200, width = 800)

par(mfrow = c(3,1))
par(mar = c(5,6,4,2))

plot(OF_out1$results$He, OF_out1$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out1$results$He[OF_out1$results$qvalues<0.01], OF_out1$results$FST[OF_out1$results$qvalues<0.01], pch=21, col="black")
title("A", cex.main = 3, adj  = 0.05, line = -5)

plot(OF_out2$results$He, OF_out2$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out2$results$He[OF_out2$results$qvalues<0.01], OF_out2$results$FST[OF_out2$results$qvalues<0.01], pch=21, col="black")
title("B", cex.main = 3, adj  = 0.05, line = -5)

plot(OF_out3$results$He, OF_out3$results$FST, pch=20, col="grey", ylab = "FST", xlab = "Expected Heterozygosity", 
     ylim = c(0, 0.45), cex.axis = 2, cex.lab = 2.2)
points(OF_out3$results$He[OF_out3$results$qvalues<0.01], OF_out3$results$FST[OF_out3$results$qvalues<0.01], pch=21, col="black")
title("C", cex.main = 3, adj  = 0.05, line = -5)

dev.off()


####getting genes nearby outliers
#loading identifier file
scaf_names <- read.table("/media/megan/New Volume/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/Hzea_genome/scaffold_and_contig_names.txt",
                         header = T)

#adding scaffold names to outlier dataframes
outliers_2002and2007_withIDs <- merge(outliers_2002and2007, scaf_names, by = "CHROM")
outliers_2007and2012_withIDs <- merge(outliers_2007and2012, scaf_names, by = "CHROM")
outliers_2012and2016_withIDs <- merge(outliers_2012and2016, scaf_names, by = "CHROM")
outliers_2002and2016_withIDs <- merge(outliers_2002and2016, scaf_names, by = "CHROM")



#Now pulling out the parts of the gff3 file that are important
gff3 <- read.table("/media/megan/New Volume/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/Hzea_genome/HzOGS2-15205-fixed_note-added.gff3",
                   sep="\t", stringsAsFactors = F)
head(gff3)

###getting nearby genes

#2002 and 2007
outliers_2002and2007_withAnnot <- merge(outliers_2002and2007_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
outliers_2002and2007_withAnnot$POS <- as.integer(as.character(outliers_2002and2007_withAnnot$POS))
outliers_2002and2007_withAnnot <- subset(outliers_2002and2007_withAnnot, V3 == "gene" & V5 > (POS -50000) & V4 < (POS + 50000))
#write.table(outliers_2002and2007_withAnnot, file = "Sig_outliers_withAnnot_2002and2007.txt", row.names = F)

#2007 and 2012
outliers_2007and2012_withAnnot <- merge(outliers_2007and2012_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
outliers_2007and2012_withAnnot$POS <- as.integer(as.character(outliers_2007and2012_withAnnot$POS))
outliers_2007and2012_withAnnot <- subset(outliers_2007and2012_withAnnot, V3 == "gene" & V5 > (POS -50000) & V4 < (POS + 50000))
#write.table(outliers_2007and2012_withAnnot, file = "Sig_outliers_withAnnot_2007and2012.txt", row.names = F)

#2012 and 2016
outliers_2012and2016_withAnnot <- merge(outliers_2012and2016_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
outliers_2012and2016_withAnnot$POS <- as.integer(as.character(outliers_2012and2016_withAnnot$POS))
outliers_2012and2016_withAnnot <- subset(outliers_2012and2016_withAnnot, V3 == "gene" & V5 > (POS -50000) & V4 < (POS + 50000))
#write.table(outliers_2012and2016_withAnnot, file = "Sig_outliers_withAnnot_2012and2016.txt", row.names = F)

#2002 and 2016
outliers_2002and2016_withAnnot <- merge(outliers_2002and2016_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
outliers_2002and2016_withAnnot$POS <- as.integer(as.character(outliers_2002and2016_withAnnot$POS))
outliers_2002and2016_withAnnot <- subset(outliers_2002and2016_withAnnot, V3 == "gene" & V5 > (POS -50000) & V4 < (POS + 50000))
#write.table(outliers_2002and2016_withAnnot, file = "Sig_outliers_withAnnot_2002and2016.txt", row.names = F)


#subsetting by top sig outliers in two by-year comparisons
outliers_firsttwo <- merge(outliers_2002and2007_withAnnot,outliers_2002and2016_withAnnot, by = "MARK")
length(unique(outliers_firsttwo$MARK))
print(unique(outliers_firsttwo$MARK))
length(unique(outliers_firsttwo$CHROM.x))
print(unique(outliers_firsttwo$CHROM.x))


outliers_nexttwo <- merge(outliers_2007and2012_withAnnot,outliers_2002and2016_withAnnot, by = "MARK")
length(unique(outliers_nexttwo$MARK))
print(unique(outliers_nexttwo$MARK))
length(unique(outliers_nexttwo$CHROM.x))
print(unique(outliers_nexttwo$CHROM.x))


outliers_lasttwo <- merge(outliers_2012and2016_withAnnot,outliers_2002and2016_withAnnot, by = "MARK")
length(unique(outliers_lasttwo$MARK))
print(unique(outliers_lasttwo$MARK))
length(unique(outliers_lasttwo$CHROM.x))
print(unique(outliers_lasttwo$CHROM.x))


##looking at shared outliers across mult by-year comparisons
outliers_all1 <- merge(outliers_firsttwo,outliers_nexttwo, by = "MARK")
length(unique(outliers_all1$MARK))
print(unique(outliers_all1$MARK))

outliers_all2 <- merge(outliers_lasttwo,outliers_nexttwo, by = "MARK")
length(unique(outliers_all2$MARK))
print(unique(outliers_all2$MARK))
#this outlier has weird changes in allele frequency over time, inconsistent with a selective sweep.



######Here I'm examining outliers for the wgs experiment using the wcFst file
wgs_data1 <- read.table("/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/FieldHzea_wcFST_all", header = F)
wgs_Hi_wcFst <- subset(wgs_data1, V5 > 0.1)
head(wgs_Hi_wcFst)

#getting scaf names
wgs_withIDs <- merge(wgs_Hi_wcFst, scaf_names, by.x = "V1", by.y = "CHROM")

wgs_wc_withAnnot <- merge(wgs_withIDs, gff3, by.x = "Scaf_ID", by.y = "V1")
str(wgs_wc_withAnnot)

wgs_wc_withAnnot <- subset(wgs_wc_withAnnot, V3.y == "gene" & V5.y > (V2.x -50000) & V4.y < (V3.x + 50000))
#write.table(wgs_wc_withAnnot, file = "/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/Hi_wcFst_annotated.txt", row.names = F)

#looking at the intersection of the two experiments
#sig_outliers_bothExp1 <- merge(outliers_2002and2007_withAnnot, wgs_wc_withAnnot, by = "V9")
#sig_outliers_bothExp2 <- merge(outliers_2007and2012_withAnnot, wgs_wc_withAnnot, by = "V9")
sig_outliers_bothExp3 <- merge(outliers_2012and2016_withAnnot, wgs_wc_withAnnot, by = "V9")
#sig_outliers_bothExp4 <- merge(outliers_2002and2016_withAnnot, wgs_wc_withAnnot, by = "V9")

#print(unique(sig_outliers_bothExp1$V9))
#print(unique(sig_outliers_bothExp2$V9))
print(unique(sig_outliers_bothExp3$V9))
#print(unique(sig_outliers_bothExp4$V9))

#write.table(unique(sig_outliers_bothExp1$V9), file = "/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/sig_outliers_wcFst0.1_bothExp1.txt", row.names = F)
#write.table(unique(sig_outliers_bothExp2$V9), file = "/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/sig_outliers_wcFst0.1_bothExp2.txt", row.names = F)
write.table(unique(sig_outliers_bothExp3$V9), file = "/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/sig_outliers_wcFst0.1_bothExp3.txt", row.names = F)
#write.table(unique(sig_outliers_bothExp4$V9), file = "/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/sig_outliers_wcFst0.1_bothExp4.txt", row.names = F)

######Now examining outliers for the wgs experiment using the pFst file
wgs_data2 <- read.table("/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/FieldHzea_pFST_all", header = F)
wgs_Hi_pFst <- subset(wgs_data2, V3 < 0.05)

#getting scaf names
wgs_withIDs_pfst <- merge(wgs_Hi_pFst, scaf_names, by.x = "V1", by.y = "CHROM")

wgs_pfst_withAnnot <- merge(wgs_withIDs_pfst, gff3, by.x = "Scaf_ID", by.y = "V1")

str(wgs_pfst_withAnnot)

wgs_pfst_withAnnot <- subset(wgs_pfst_withAnnot, V3.y == "gene" & V5 > (V2.x -50000) & V4 < (V3.x + 50000))
#write.table(wgs_pfst_withAnnot, file = "/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/Hi_pfst_annotated.txt", row.names = F)


#looking at the intersection of the two experiments
#sig_outliers_bothExp1_pfst <- merge(outliers_2002and2007_withAnnot, wgs_pfst_withAnnot, by = "V9")
#sig_outliers_bothExp2_pfst <- merge(outliers_2007and2012_withAnnot, wgs_pfst_withAnnot, by = "V9")
sig_outliers_bothExp3_pfst <- merge(outliers_2012and2016_withAnnot, wgs_pfst_withAnnot, by = "V9")
#sig_outliers_bothExp4_pfst <- merge(outliers_2002and2016_withAnnot, wgs_pfst_withAnnot, by = "V9")

#print(unique(sig_outliers_bothExp1$V9))
#print(unique(sig_outliers_bothExp2$V9))
print(unique(sig_outliers_bothExp3_pfst$V9))
#print(unique(sig_outliers_bothExp4$V9))

#write.table(unique(sig_outliers_bothExp1_pfst$V9), file = "/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/sig_outliers_pFst0.1_bothExp1.txt", row.names = F)
#write.table(unique(sig_outliers_bothExp2_pfst$V9), file = "/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/sig_outliers_pFst0.1_bothExp2.txt", row.names = F)
write.table(unique(sig_outliers_bothExp3_pfst$V9), file = "/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/sig_outliers_pFst0.05_bothExp3.txt", row.names = F)
#write.table(unique(sig_outliers_bothExp4_pfst$V9), file = "/media/megan/New Volume/Hzea_WGRS_Bowtie2_output/WGRS_mpileupANDvcftools_output_v2/sig_outliers_pFst0.1_bothExp4.txt", row.names = F)
