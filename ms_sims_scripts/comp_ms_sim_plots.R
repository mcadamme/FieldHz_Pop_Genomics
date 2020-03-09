#line plot comparing FST distributions from different simulations.

library(ggplot2); library(dplyr)

#datasets comparing different generation sizes
data_2 <- read.table("~/ms_sims/FieldPops/Ne_fivethou_noSize/FST_msSims1.out", header = F)
data_1 <- read.table("~/ms_sims/FieldPops/Ne_twopointfivethou_NoSize/FST_msSims2.out", header = F)
data_3 <- read.table("~/ms_sims/FieldPops/Ne_tenthou_NoSize/FST_msSims3.out", header = F)
data_4 <- read.table("~/ms_sims/FieldPops/Ne_twentythou_NoSize/FST_msSims4.out", header = F)


df <- rbind(data_4, data_3, data_2, data_1)
sim_name <- rep(c("4","3","2","1"), times = 1, each = 10000)

df2 <- cbind(df, sim_name)
str(df2)

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
mu2$mu <- as.numeric(as.character(mu2$mu))


#comparing all datasets
png(filename = "~/ms_sims/Dist_Fst_Neutral_Final_Lines_diffPopSize_dens.png", units = "px", height = 600, width = 700)
a <- ggplot(df3, aes(x = df3$V1, fill = df3$sim_name)) +
  scale_fill_manual(name = " Present Pop Size", labels = c("n=20,000", "n=10,000", "n=5,000", "n=2,500"), values = c("#868686FF","#0073C2FF","#ff1f39","#EFC000FF")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x=expression(F["ST"]*" value"), y=expression("Density")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) 

#a + geom_histogram(binwidth = 0.005, alpha = 0.4) +
  a + geom_density(mapping = aes(color = df3$sim_name), data = df3, stat = "density", alpha = 0.1) +
  scale_color_manual(name = expression("Mean F"["ST"]), labels = c("n=20,000", "n=10,000", "n=5,000", "n=2,500"), values = c("#868686FF","#0073C2FF","#ff1f39","#EFC000FF")) +
  geom_vline(aes(xintercept = mu, color = sims), data = mu2, linetype = "dashed")

dev.off()


