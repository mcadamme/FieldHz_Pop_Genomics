#this is the script I used for my DAPC
#I used this to try to identify SNPs that explain the greatest genetic divergence across years
# 02/11/2018 MF


setwd("/media/megan/New Volume1/Hz_PopGen_ddRAD_demult/Bowtie_genome_alignments")

library(adegenet); library(pegas); library(ade4); library(vcfR); library(ape); library(parallel); library(reshape2); library(ggplot2)

#first reading in data produced by BCFtools
data <- read.vcfR(file = "./mpileupANDvcftools_output/thinned_FieldHzea_allpops.recode.vcf")
head(data)
data@fix[1:10,1:5]

Hz.genlight <- vcfR2genlight(data, n.cores=2)

pop(Hz.genlight)<-regmatches(indNames(Hz.genlight), regexpr("20[[:digit:]]+", indNames(Hz.genlight)))
head(pop(Hz.genlight)) #checking to see output makes sense

Hz.genlight@ploidy <- as.integer(ploidy(Hz.genlight)) 

#Examining using PCA clustering of populations
pca1 <- glPca(x = Hz.genlight, n.cores = 6) #getting PCs - based upon figure prompt, the first 6 PCs seem to explain a good deal genet var.
6

#axes 1&2
pca1_colplot <- colorplot(pca1$scores[,1:2],pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:50],2,1,2, posi="topright", inset=.05, ratio=.3)

#axes 1&3
pca1_colplot2 <- colorplot(pca1$scores[,c(1,3)],pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:50],2,1,2, posi="topright", inset=.05, ratio=.3)

#axes 1&4
pca1_colplot3 <- colorplot(pca1$scores[,c(1,4)],pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:50],2,1,2, posi="topright", inset=.05, ratio=.3)


#getting the proportion of variation explained by 225 PCs
(sum(pca1$eig[1:225]))/(sum(pca1$eig))
#getting the proportion of variation explained by 50 PCs
(sum(pca1$eig[1:50]))/(sum(pca1$eig))
#getting the proportion of variation explained by first 6 PCs
(sum(pca1$eig[1:6]))/(sum(pca1$eig))


#Checking output of find clusters for different "numbers of clusters" to validate pairwise Fst results with PCA and identify
#SNPs that explain divergence between populations.

Hz_grp2 <- find.clusters(Hz.genlight, max.n.clust=10, n.pca = 225, n.cores = 6) #chose to retain 225 PCs & 2 clusters - BtRes and BtSus.
2

Hz_grp3 <- find.clusters(Hz.genlight, max.n.clust=10, n.pca = 225, n.clust = 3, n.cores = 6)
Hz_grp4 <- find.clusters(Hz.genlight, max.n.clust=10, n.pca = 225, n.clust = 4, n.cores = 6) #chose to retain 225 PCs & 4 clusters - 1 per original pop.
Hz_grp5 <- find.clusters(Hz.genlight, max.n.clust=10, n.pca = 225, n.clust = 5, n.cores = 6)

png(filename = "./mpileupANDvcftools_output/nclusts.png", units = "px", height = 600, width = 400)
par(mfrow = c(4,1))
nclust2 <- table.value(table(pop(Hz.genlight), Hz_grp2$grp), col.lab=paste("inf", 1:2), clabel.col = 1.8, clabel.row = 1.8, clegend = 1.5,
                       row.lab=paste("ori", 1:4)) 
nclust3 <- table.value(table(pop(Hz.genlight), Hz_grp3$grp), col.lab=paste("inf", 1:3), clabel.col = 1.8, clabel.row = 1.8, clegend = 1.5,
                       row.lab=paste("ori", 1:4))
nclust4 <- table.value(table(pop(Hz.genlight), Hz_grp4$grp), col.lab=paste("inf", 1:4), clabel.col = 1.8, clabel.row = 1.8, clegend = 1.5,
                       row.lab=paste("ori", 1:4))
nclust5 <- table.value(table(pop(Hz.genlight), Hz_grp5$grp), col.lab=paste("inf", 1:5), clabel.col = 1.8, clabel.row = 1.8, clegend = 1.5,
                       row.lab=paste("ori", 1:4))
dev.off()

#these plots show why one or two clusters are the best fit. 
#R cannot assign individuals of a single original population to their own inferred population.
#Most individuals are similar enough that they get assigned to the first two inferred populations, even when 4 clusters are specified.


#Now for a DAPC

#setting color scheme
myCol <- c("#1B9E77", "#A6761D", "#E6AB02", "#7570B3")


#cross-validating number PCs to retain. Using K = 2 groupings.
set.seed(999)
Hz.genlight_2 <- xvalDapc(tab(Hz.genlight, NA.method = "mean"), Hz_grp2$grp, parallel = "multicore")
Hz.genlight_2[-1] 

set.seed(999)
system.time(Hz.genlight_2 <- xvalDapc(tab(Hz.genlight, NA.method = "mean"), Hz_grp2$grp,
                                      n.pca = 5:25, n.rep = 100,
                                      parallel = "multicore", ncpus = 6L))
Hz.genlight_2[-1] #6 PCs looks best, but 7 & 8 are still high

#cross-validating number PCs to retain. Using K = 3 groupings.
set.seed(999)
Hz.genlight_3 <- xvalDapc(tab(Hz.genlight, NA.method = "mean"), Hz_grp3$grp, parallel = "multicore")
Hz.genlight_3[-1] 

set.seed(999)
system.time(Hz.genlight_3 <- xvalDapc(tab(Hz.genlight, NA.method = "mean"), Hz_grp3$grp,
                                      n.pca = 5:25, n.rep = 100,
                                      parallel = "multicore", ncpus = 6L))
Hz.genlight_3[-1] #5 & 7 PCs very slightly better successful assignment, but 6 &8 are close

#cross-validating number PCs to retain. Using K = 4 groupings.
set.seed(999)
Hz.genlight_4 <- xvalDapc(tab(Hz.genlight, NA.method = "mean"), Hz_grp4$grp, parallel = "multicore")
Hz.genlight_4[-1] 

set.seed(999)
system.time(Hz.genlight_4 <- xvalDapc(tab(Hz.genlight, NA.method = "mean"), Hz_grp4$grp,
                                      n.pca = 5:25, n.rep = 100,
                                      parallel = "multicore", ncpus = 6L))
Hz.genlight_4[-1] # 8 PCs best, much better than lower....



#Using by-year (prior) clustering.
set.seed(999)
Hz.genlight_x <- xvalDapc(tab(Hz.genlight, NA.method = "mean"), Hz.genlight$pop, parallel = "multicore")
Hz.genlight_x[-1] #shows that the numbers of successful assignments are highest between 100 and 120 PCs.

set.seed(999)
system.time(Hz.genlight_x <- xvalDapc(tab(Hz.genlight, NA.method = "mean"), Hz.genlight$pop,
                              n.pca = 98:118, n.rep = 100,
                              parallel = "multicore", ncpus = 6L))
Hz.genlight_x[-1] #115 is best number PCs


#DAPC to get compoplot - first looking at original population groupings
dapc1<- dapc(Hz.genlight, Hz.genlight$pop, n.pca = 115, n.da = 10, var.loadings=T, pca.info=T, n.cores = 6) #with prior group assignments, cross-val PC num

png(filename = "./mpileupANDvcftools_output/Fig1_DAPC_withPrior_Pop_Assignments.png", units = "px", height = 600, width = 800)
scatter(dapc1, scree.da=TRUE, bg="white", posi.pca="topright", legend=FALSE, col=myCol)
dev.off()

compoplot(dapc1, col=myCol,lab="", subset = (rownames(dapc1$ind.coord))[c(1:259)], txt.leg=paste("group", 1:4), cex.lab = 1.5, cex.axis = 1.2) #more blue on the far left than on far right.

dapc1.results <- as.data.frame(dapc1$posterior)
dapc1.results$pop <- Hz.genlight$pop
dapc1.results$indNames <- Hz.genlight$ind.names

dapc1.results <- melt(dapc1.results)

colnames(dapc1.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc1.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = myCol) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + labs(y="Posterior Membership Probability", fill  = "Assigned Pop")
p <- p + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.y = element_text(size = (16),margin = margin(t = 0, r = 15, b = 0, l = 0)),
               axis.text.y = element_text(size = (12)),
               legend.title = element_text(size = (14)),
               legend.text = element_text(size = (12)),
               strip.text.x = element_text(size = 12, color = "black", face = "bold"))
p

loadingplot(dapc1$var.load, lab.jitter=1, threshold = quantile(dapc1$var.load, prob=0.999))


#This is a DAPC using the clusters derived from the k-means analysis
#Want to know whether there are any SNPs that are shared between DAPC1 as outliers, and DAPC2.
#DAPC2 specifies clustering into two populations, and I wonder if it is predictive of Bt resistance. 
dapc2 <- dapc(Hz.genlight, Hz_grp2$grp, n.pca = 8, n.da = 10, n.cores = 6, var.loadings = T)
dapc2
scatter(dapc2,scree.da=TRUE, bg="white", posi.pca="topright", legend=TRUE,txt.leg=paste("group", 1:2), col=myCol)

dapc2.results <- as.data.frame(dapc2$posterior)
dapc2.results$pop <- Hz.genlight$pop
dapc2.results$indNames <- Hz.genlight$ind.names

dapc2.results <- melt(dapc2.results)

colnames(dapc2.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc2.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = myCol) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + labs(y="Posterior Membership Probability", fill  = "Assigned Pop")
p <- p + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.y = element_text(size = (16),margin = margin(t = 0, r = 15, b = 0, l = 0)),
               axis.text.y = element_text(size = (12)),
               legend.title = element_text(size = (14)),
               legend.text = element_text(size = (12)),
               strip.text.x = element_text(size = 12, color = "black", face = "bold"))
p

#getting SNPs that distinguish between these two pops
loadingplot(dapc2$var.load, lab.jitter=1, threshold = quantile(dapc2$var.load, prob=0.995))


#DAPC3 specifies clustering into three populations by kmeans analysis.

dapc3 <- dapc(Hz.genlight, Hz_grp3$grp, n.pca = 8, n.da = 10, n.cores = 6, var.loadings = T)
dapc3
scatter(dapc3,scree.da=TRUE, bg="white", posi.pca="topright", legend=TRUE,txt.leg=paste("group", 1:3), col=myCol)

dapc3.results <- as.data.frame(dapc3$posterior)
dapc3.results$pop <- Hz.genlight$pop
dapc3.results$indNames <- Hz.genlight$ind.names

dapc3.results <- melt(dapc3.results)

colnames(dapc3.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc3.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = myCol) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + labs(y="Posterior Membership Probability", fill  = "Assigned Pop")
p <- p + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.y = element_text(size = (16),margin = margin(t = 0, r = 15, b = 0, l = 0)),
               axis.text.y = element_text(size = (12)),
               legend.title = element_text(size = (14)),
               legend.text = element_text(size = (12)),
               strip.text.x = element_text(size = 12, color = "black", face = "bold"))
p

#DAPC4 specifies clustering into four populations by kmeans analysis.

dapc4 <- dapc(Hz.genlight, Hz_grp4$grp, n.pca = 8, n.da = 10, n.cores = 6, var.loadings = T)
dapc4
scatter(dapc4,scree.da=TRUE, bg="white", posi.pca="topright", legend=TRUE,txt.leg=paste("group", 1:4), col=myCol)

dapc4.results <- as.data.frame(dapc4$posterior)
dapc4.results$pop <- Hz.genlight$pop
dapc4.results$indNames <- Hz.genlight$ind.names

dapc4.results <- melt(dapc4.results)

colnames(dapc4.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc4.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = myCol) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + labs(y="Posterior Membership Probability", fill  = "Assigned Pop")
p <- p + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.y = element_text(size = (16),margin = margin(t = 0, r = 15, b = 0, l = 0)),
               axis.text.y = element_text(size = (12)),
               legend.title = element_text(size = (14)),
               legend.text = element_text(size = (12)),
               strip.text.x = element_text(size = 12, color = "black", face = "bold"))
p
