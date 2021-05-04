#getting overlap between ddRAD outliers and QTL
#04212021 MF

setwd("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS")

#loading zea superscaffolds
Full_Ord_Scafs <- read.table(file = "Hzea_superScaf_genome.txt", header = T)

#loading ddRAD outlier datasets
FSToutliers_allyears <- read.table(file = "FileS1_outliers_allyears.txt", header = T)
head(FSToutliers_allyears)

DAPC_top5 <- read.table(file = "DAPC2_top5per.txt", header = F)
head(DAPC_top5)

#loading and filtering QTL data
sig_BAP11A1_DD <- read.csv(file = "sig_BAP11A1_DD.csv", header = T)
sig_DEO9A1_DD <- read.csv(file = "sig_DEO9A1_DD.csv", header = T)

sub_sig_BAP11A1_DD <- subset(sig_BAP11A1_DD, gamma > 0.01 & sigLev == 0.01)
head(sub_sig_BAP11A1_DD)
nrow(sub_sig_BAP11A1_DD)
sub_sig_DEO9A1_DD <- subset(sig_DEO9A1_DD, gamma > 0.01 & sigLev == 0.01)
head(sub_sig_DEO9A1_DD)
nrow(sub_sig_DEO9A1_DD)


