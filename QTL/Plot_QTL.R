#Used this script to plot QTL & produce csv with SNPs at different sig levels 
#sig level data will be input for Reanalysis_Manhattan.R
#03302021 MF

library(CMplot); library(plyr)

setwd("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS")

#loading zea superscaffolds
Full_Ord_Scafs <- read.table(file = "Hzea_superScaf_genome.txt", header = T)

##### Function to subset Sig QTL at diff thresholds #####

sig_levels <- list(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001)
df <- data.frame()

getSig_QTL <- function(X, outfile){
  for (i in sig_levels){
    print(i)
    X$p_lrt <- as.numeric(as.character(X$p_lrt))
    temp_df <- subset(X, p_lrt < i)
    temp_df$sigLev <- rep(i, times = nrow(temp_df))
    print(head(temp_df))
    df <- rbind(df, temp_df)
    
  }
write.csv(df, file = outfile)
}

##### loading QTL SNPs & removing low freq allele s#####

BAP11A1_DD <- read.csv("BAP11A1_DD_all_snps.csv", header = T)#BC0805
str(BAP11A1_DD)
sub_BAP11A1_DD <- subset(BAP11A1_DD, allele_frequency > 0.099)

BAP11A1_CL <- read.csv("BAP11A1_CL_all_snps.csv", header = T)#providence
str(BAP11A1_CL)
sub_BAP11A1_CL <- subset(BAP11A1_CL, allele_frequency > 0.099)

DEO9A1_DD <- read.csv("DEO9A1_DD_all_snps.csv", header = T)#ObsII
str(DEO9A1_DD)
sub_DEO9A1_DD <- subset(DEO9A1_DD, allele_frequency > 0.099)

DEO9A1_CL <- read.csv("DEO9A1_CL_all_snps.csv", header = T)#Obs
str(DEO9A1_CL)
sub_DEO9A1_CL <- subset(DEO9A1_CL, allele_frequency > 0.099)


###### Applying getSig_QTL function ######

#one toxin family
getSig_QTL(sub_BAP11A1_DD, "sig_BAP11A1_DD.csv")
getSig_QTL(sub_BAP11A1_CL, "sig_BAP11A1_CL.csv")

#two toxin family
getSig_QTL(sub_DEO9A1_DD, "sig_DEO9A1_DD.csv")
getSig_QTL(sub_DEO9A1_CL, "sig_DEO9A1_CL.csv")

#getting numbers of sig SNPs per sig threshold
table(read.csv("sig_BAP11A1_DD.csv", header = T)[,14])
table(read.csv("sig_BAP11A1_CL.csv", header = T)[,14])
table(read.csv("sig_DEO9A1_DD.csv", header = T)[,14])
table(read.csv("sig_DEO9A1_CL.csv", header = T)[,14])

temp_df <- subset(DEO9A1_CL, p_lrt < 0.0001)#sanity check for last table.  Only obs at this thresh has low allele freq.


#checking on the shared scaffolds between CL and DD
sig_BAP11A1_DD <- read.csv("sig_BAP11A1_DD.csv", header = T)
sig_BAP11A1_CL <- read.csv("sig_BAP11A1_CL.csv", header = T)
sig_DEO9A1_DD <- read.csv("sig_BAP11A1_DD.csv", header = T)
sig_DEO9A1_CL <- read.csv("sig_BAP11A1_CL.csv", header = T)

#which scaffolds are shared between DD and CL treatments?
intersect(unique(sig_BAP11A1_CL$scaffold), unique(sig_BAP11A1_DD$scaffold))

#which scaffolds are shared between DD and CL treatments?
intersect(unique(sig_DEO9A1_CL$scaffold), unique(sig_DEO9A1_DD$scaffold))


##### QTL plotting function #####
plot_QTL <- function(X, name, sigLev){
  
merged_Map <- merge(Full_Ord_Scafs, X, by.x = "Scaf", by.y = "scaffold")
merged_Map <- merged_Map[order(merged_Map[,5], merged_Map[,6], merged_Map[,3], 
                               merged_Map[,1], merged_Map[,9]),]
with_artPos <- data.frame()

  for (i in seq(c(1:max(merged_Map$Chr)))){
    print(i)
    sub_chrom <- subset(merged_Map, Chr == i)
    artPos <- seq(0, (nrow(sub_chrom)-1), by = 1)
    comb <- data.frame(cbind(sub_chrom, artPos))
    with_artPos <- rbind(with_artPos, comb)
  }


For_plot <- data.frame(cbind(with_artPos$snp_id, with_artPos$Chr, with_artPos$radtag_position, with_artPos$p_lrt))
names(For_plot) <- c("SnpName", "Chr", "Pos", "Pval")
str(For_plot)#sanity check

SNPs <- subset(For_plot, Pval < sigLev)
CMplot(For_plot, type="p",plot.type="m", LOG10=TRUE, threshold=NULL,file="jpg",memo=name,dpi=300,
       col=c("grey50","grey70"), pch = 21, cex = 0.8, ylim = c(0,12), highlight=SNPs$SnpName,
       highlight.col="dark blue", highlight.cex=0.8, highlight.pch=19, file.output=TRUE,verbose=TRUE,width=14,height=6)
}

##### Drawing plots #####

#plot_QTL runs once per dataset - last argument is sig level
plot_QTL(sub_BAP11A1_DD, "BAP11A1_DD", 1e-3)
plot_QTL(sub_BAP11A1_CL, "BAP11A1_CL", 1e-3)
plot_QTL(sub_DEO9A1_DD, "DEO9A1_DD", 1e-3)
plot_QTL(sub_DEO9A1_CL, "DEO9A1_CL", 1e-3)

