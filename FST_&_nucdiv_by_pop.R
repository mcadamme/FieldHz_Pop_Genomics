#script to get pairwise pop-level FST values.

setwd("/media/megan/New Volume1/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments")

library(adegenet); library(pegas); library(ade4); library(vcfR); library(StAMPP); library(ape); library(parallel)

#first reading in data produced by BCFtools
data <- read.vcfR(file = "./mpileupANDvcftools_output/thinned_FieldHzea_allpops.recode.vcf")
head(data)
data@fix[1:10,1:5]

Hz.genlight <- vcfR2genlight(data, n.cores=2)

pop(Hz.genlight)<-regmatches(indNames(Hz.genlight), regexpr("20[[:digit:]]+", indNames(Hz.genlight)))
head(pop(Hz.genlight)) #checking to see output makes sense

Hz.genlight@ploidy <- as.integer(ploidy(Hz.genlight)) 

#Getting a prelim look at SNPs in BCFtools dataset.
glPlot(Hz.genlight, posi="topleft")


#calculating Fst values
Hz.fst<-stamppFst(Hz.genlight, nboots = 500, percent =95, nclusters=4) 

Hz.fst$Fsts
Hz.fst$Pvalues

#Nucleotide diversity plus confidence intervals for all pops
#loading datasets
pi_2002 <- read.table("./mpileupANDvcftools_output/FieldHzea2002.sites.pi", header = T)
pi_2007 <- read.table("./mpileupANDvcftools_output/FieldHzea2007.sites.pi", header = T)
pi_2012 <- read.table("./mpileupANDvcftools_output/FieldHzea2012.sites.pi", header = T)
pi_2016 <- read.table("./mpileupANDvcftools_output/FieldHzea2016.sites.pi", header = T)

het_2002 <-read.table("./mpileupANDvcftools_output/FieldHzea2002.het", header = T)
het_2007 <-read.table("./mpileupANDvcftools_output/FieldHzea2007.het", header = T)
het_2012 <-read.table("./mpileupANDvcftools_output/FieldHzea2012.het", header = T)
het_2016 <-read.table("./mpileupANDvcftools_output/FieldHzea2016.het", header = T)

#bootstrapping function
boot.fn <- function(x, N=5000) {
  Int.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  Int.CI <- quantile(Int.1, probs=c(0.025,0.975))
  Int.CI
}

mean(pi_2002$PI)
boot.fn(pi_2002$PI)

mean(pi_2007$PI)
boot.fn(pi_2007$PI)

mean(pi_2012$PI)
boot.fn(pi_2012$PI)

mean(pi_2016$PI)
boot.fn(pi_2016$PI)


mean(het_2002$F)
boot.fn(het_2002$F)

mean(het_2007$F)
boot.fn(het_2007$F)

mean(het_2012$F)
boot.fn(het_2012$F)

mean(het_2016$F)
boot.fn(het_2016$F)
