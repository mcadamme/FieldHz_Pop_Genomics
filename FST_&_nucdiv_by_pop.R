#script to get pairwise pop-level FST values.

setwd("/media/megan/New Volume/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/mpileupANDvcftools_output")

library(adegenet); library(pegas); library(ade4); library(vcfR); library(StAMPP)

data <- read.vcfR(file = "thinned_FieldHzea_allpops.recode.vcf")
head(data)
data@fix[1:10,1:5]

aa.genlight <- vcfR2genlight(data, n.cores=1)

pop(aa.genlight)<-regmatches(indNames(aa.genlight), regexpr("20[[:digit:]]+", indNames(aa.genlight)))
head(pop(aa.genlight)) #checking to see output makes sense

aa.genlight@ploidy <- as.integer(ploidy(aa.genlight)) 
aa.fst<-stamppFst(aa.genlight, nboots = 500, percent =95, nclusters=4) 

aa.fst$Fsts
aa.fst$Pvalues

####avg nucleotide diversity plus confidence intervals for all pops
#loading datasets
pi_2002 <- read.table("FieldHzea2002.sites.pi", header = T)
pi_2007 <- read.table("FieldHzea2007.sites.pi", header = T)
pi_2012 <- read.table("FieldHzea2012.sites.pi", header = T)
pi_2016 <- read.table("FieldHzea2016.sites.pi", header = T)

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
