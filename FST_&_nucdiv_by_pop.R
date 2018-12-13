#script to get pairwise pop-level FST values.

setwd("/media/megan/New Volume/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments/mpileupANDvcftools_output")

library(adegenet); library(pegas); library(ade4); library(vcfR); library(StAMPP)

data <- read.vcfR(file = "thinned_FieldHzea_allpops.recode.vcf")
head(data)
data@fix[1:10,1:5]

Hz.genlight <- vcfR2genlight(data, n.cores=2)

pop(Hz.genlight)<-regmatches(indNames(Hz.genlight), regexpr("20[[:digit:]]+", indNames(Hz.genlight)))
head(pop(Hz.genlight)) #checking to see output makes sense

Hz.genlight@ploidy <- as.integer(ploidy(Hz.genlight)) 
Hz.fst<-stamppFst(Hz.genlight, nboots = 500, percent =95, nclusters=4) 

Hz.fst$Fsts
Hz.fst$Pvalues

#Original command to get line plot. At prompt, put in 225 PCs and 4 clusters
#Hz_grp <- find.clusters(Hz.genlight, max.n.clust=20)

#Checking output of find clusters for 3 different "numbers of clusters" to validate pairwise Fst results with PCA
Hz_grp2 <- find.clusters(Hz.genlight, max.n.clust=20, n.pca = 225, n.clust = 2) #chose to retain 225 PCs & 2 clusters - BtRes and BtSus.
Hz_grp3 <- find.clusters(Hz.genlight, max.n.clust=20, n.pca = 225, n.clust = 3)
Hz_grp4 <- find.clusters(Hz.genlight, max.n.clust=20, n.pca = 225, n.clust = 4) #chose to retain 225 PCs & 4 clusters - 1 per original pop.


#these plots show why one or two clusters are the best fit. 
#R cannot assign individuals of a single original population to their own inferred population.
#Most individuals are similar enough that they get assigned to the first two inferred populations, even when 4 clusters are specified.

png(filename = "nclusts.png", units = "px", height = 600, width = 400)
par(mfrow = c(3,1))
nclust2 <- table.value(table(pop(Hz.genlight), Hz_grp2$grp), col.lab=paste("inf", 1:2), clabel.col = 1.8, clabel.row = 1.8, clegend = 1.5,
            row.lab=paste("ori", 1:4)) 
nclust3 <- table.value(table(pop(Hz.genlight), Hz_grp3$grp), col.lab=paste("inf", 1:3), clabel.col = 1.8, clabel.row = 1.8, clegend = 1.5,
                       row.lab=paste("ori", 1:4))
nclust4 <- table.value(table(pop(Hz.genlight), Hz_grp4$grp), col.lab=paste("inf", 1:4), clabel.col = 1.8, clabel.row = 1.8, clegend = 1.5,
                       row.lab=paste("ori", 1:4))

dev.off()

#This is the PCA showing the plots of the populations discriminated using PC1 and PC2

dapc1 <- dapc(Hz.genlight, Hz_grp2$grp, n.pca = 225, n.da = 10)
dapc1
scatter(dapc1)

dapc2 <- dapc(Hz.genlight, Hz_grp3$grp, n.pca = 225, n.da = 10)
dapc2
scatter(dapc2)

dapc3 <- dapc(Hz.genlight, Hz_grp4$grp, n.pca = 225, n.da = 10)
dapc3
scatter(dapc3)

myCol <- c("darkblue","darkgreen","orange","purple")

scatter(dapc1, ratio.pca=0.3, bg="white", pch=20,  cell=0,
        cstar=0, col=myCol[1:2], solid=.4, cex=3, clab=0,
        mstree=T, scree.da=TRUE, posi.pca="topright",
        leg=TRUE, posi.leg = "topleft", txt.leg=paste("Cluster",1:2))

scatter(dapc2, ratio.pca=0.3, bg="white", pch=20,  cell=0,
        cstar=0, col=myCol[1:3], solid=.4, cex=3, clab=0,
        mstree=T, scree.da=TRUE, posi.pca="topright",
        leg=TRUE, posi.leg = "topleft", txt.leg=paste("Cluster",1:3))

scatter(dapc3, ratio.pca=0.3, bg="white", pch=20,  cell=0,
        cstar=0, col=myCol, solid=.4, cex=3, clab=0,
        mstree=T, scree.da=TRUE, posi.pca="topright",
        leg=TRUE, posi.leg = "bottomleft", txt.leg=paste("Cluster",1:4))

#looking at loadings
loadingplot(dapc2$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)
dapc2$posterior

big_contrib <- dapc2$var.contr > 0.002

#When I allow for more than 2 clusters, it seems as though there are a few individuals from 2002 that do cluster together
#and away from the other individuals in the overall population. 
#Perhaps those individuals are migrants from South America, whereas the others are all from North American populations that have
#overwintered.  The two individuals are HZ_2002_8_37, and HZ_2002_early_33.  

#Nucleotide diversity plus confidence intervals for all pops
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
