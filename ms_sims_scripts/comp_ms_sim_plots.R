#line plot comparing FST distributions from different simulations.
#MF 03/05/2020; updated 03/30/2021

library(ggplot2); library(dplyr)

#datasets comparing different generation sizes - set to 10kb sims, but can be changed to 40kb sims
data_1 <- read.table("~/ms_sims/FieldPops/N0_21thou_10kb/FST_msSims1.out", header = F)
data_2 <- read.table("~/ms_sims/FieldPops/N0_63thou_10kb/FST_msSims2.out", header = F)
data_3 <- read.table("~/ms_sims/FieldPops/N0_105thou_10kb/FST_msSims3.out", header = F)
data_4 <- read.table("~/ms_sims/FieldPops/N0_140thou_10kb/FST_msSims4.out", header = F)


df <- data.frame(c(data_1[c(1:10000),], data_2[c(1:10000),], data_3[c(1:10000),], data_4[c(1:10000),]))
sim_name <- rep(c("4","3","2","1"), times = 1, each = 10000)

df2 <- cbind(df, sim_name)
str(df2)
names(df2) <- c("V1", "sim_name")

df3 <- na.omit(df2)

#getting max fst val for each dataset
maxFST <- tapply(df3$V1, df3$sim_name, max)
maxFST

varFST <- tapply(df3$V1, df3$sim_name, var)
varFST

mu <- tapply(df3$V1, df3$sim_name, mean)
mu

mu <-unlist(mu)
sims <- c("4","3","2","1")
mu2 <- data.frame(cbind(mu, sims))
mu2$mu <-rev(as.numeric(as.character(mu2$mu)))


#comparing all datasets
png(filename = "~/ms_sims/Dist_Fst_Neutral_Final_Lines_diffPopSize_dens.png", units = "px", height = 600, width = 700)
a <- ggplot(df3, aes(x = df3$V1, fill = df3$sim_name)) +
  scale_fill_manual(name = " Present Pop Size", labels = c("n=140,000", "n=105,000", "n=63,000", "n=21,000"), values = c("#868686FF","#0073C2FF","#ff1f39","#EFC000FF")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x=expression(F["ST"]*" value"), y=expression("Density")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) 

  a + geom_density(mapping = aes(color = df3$sim_name), data = df3, stat = "density", alpha = 0.3) +
  scale_color_manual(name = expression("Mean F"["ST"]), labels = c("n=140,000", "n=105,000", "n=63,000", "n=21,000"), values = c("#868686FF","#0073C2FF","#ff1f39","#EFC000FF")) +
  geom_vline(aes(xintercept = mu, color = sims), data = mu2, linetype = "dashed")

dev.off()

#getting p-vals - relative to 6 s.d. threshold @10kb window for 2002-2017 
nrow(subset(data_1, V1 > 0.08))/nrow(na.omit(data_1))
nrow(subset(data_2, V1 > 0.08))/nrow(na.omit(data_2))
nrow(subset(data_3, V1 > 0.08))/nrow(na.omit(data_3))
nrow(subset(data_4, V1 > 0.08))/nrow(na.omit(data_4))

#relative to 6 s.d. threshold @40kb window for 2002-2017 
nrow(subset(data_1, V1 > 0.047))/nrow(na.omit(data_1))
nrow(subset(data_2, V1 > 0.047))/nrow(na.omit(data_2))
nrow(subset(data_3, V1 > 0.047))/nrow(na.omit(data_3))
nrow(subset(data_4, V1 > 0.047))/nrow(na.omit(data_4))
