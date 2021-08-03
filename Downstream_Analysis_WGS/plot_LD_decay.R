#plotting ld decay
#04192020 MF

#import the data
LD_2002 <- data.frame(read.table("/media/megan/easystore/Hz2002_plink.ld", header=T))
LD_2012 <- data.frame(read.table("/media/megan/easystore/Hz2012_plink.ld", header=T))
LD_2017 <- data.frame(read.table("/media/megan/easystore/Hz2017_plink.ld", header=T))

###### subsetting data to include both "neutral" and "swept" scaffolds ######

#2002
LD_2002_in3Mb.1 <- subset(LD_2002, CHR_A == "KZ118395.1")#top two scafs in 3Mb sweep on Chr5
LD_2002_in3Mb.2 <- subset(LD_2002, CHR_A == "KZ117131.1")

LD_2002_Cry1Chr6 <- subset(LD_2002, CHR_A == "KZ118067.1")
LD_2002_Cry1Chr7.1 <- subset(LD_2002, CHR_A == "KZ116099.1")
LD_2002_Cry1Chr7.2 <- subset(LD_2002, CHR_A == "KZ117975.1")
#LD_2002_Cry1Chr11 <- subset(LD_2002, CHR_A == "KZ118133.1") #this is the new 4

LD_2002_Cry2Chr7 <- subset(LD_2002, CHR_A == "KZ117108.1")

LD_2002_notSel <- subset(LD_2002, CHR_A != "KZ118395.1" & CHR_A != "KZ118483.1" & 
                           CHR_A != "KZ117237.1" & CHR_A != "KZ117131.1" & CHR_A != "KZ116099.1" &
                           CHR_A != "KZ117975.1" & CHR_A != "KZ118067.1" & CHR_A != "KZ117108.1")
head(LD_2002_notSel)


#2012
LD_2012_in3Mb.1 <- subset(LD_2012, CHR_A == "KZ118395.1")#top two scafs in 3Mb sweep on Chr5
LD_2012_in3Mb.2 <- subset(LD_2012, CHR_A == "KZ117131.1")

LD_2012_Cry1Chr6 <- subset(LD_2012, CHR_A == "KZ118067.1")
LD_2012_Cry1Chr7.1 <- subset(LD_2012, CHR_A == "KZ116099.1")
LD_2012_Cry1Chr7.2 <- subset(LD_2012, CHR_A == "KZ117975.1")
#LD_2012_Cry1Chr11 <- subset(LD_2012, CHR_A == "KZ118133.1")

LD_2012_Cry2Chr7 <- subset(LD_2012, CHR_A == "KZ117108.1")

LD_2012_notSel <- subset(LD_2012, CHR_A != "KZ118395.1" & CHR_A != "KZ118483.1" & 
                           CHR_A != "KZ117237.1" & CHR_A != "KZ117131.1" & CHR_A != "KZ116099.1" &
                           CHR_A != "KZ117975.1" & CHR_A != "KZ118067.1" & CHR_A != "KZ117108.1")
head(LD_2012_notSel)


#2017
LD_2017_in3Mb.1 <- subset(LD_2017, CHR_A == "KZ118395.1")#top two scafs in 3Mb sweep on Chr5
LD_2017_in3Mb.2 <- subset(LD_2017, CHR_A == "KZ117131.1")

LD_2017_Cry1Chr6 <- subset(LD_2017, CHR_A == "KZ118067.1")
LD_2017_Cry1Chr7.1 <- subset(LD_2017, CHR_A == "KZ116099.1")
LD_2017_Cry1Chr7.2 <- subset(LD_2017, CHR_A == "KZ117975.1")
#LD_2017_Cry1Chr11 <- subset(LD_2017, CHR_A == "KZ118133.1")

LD_2017_Cry2Chr7 <- subset(LD_2017, CHR_A == "KZ117108.1")

LD_2017_notSel <- subset(LD_2017, CHR_A != "KZ118395.1" & CHR_A != "KZ118483.1" & 
                           CHR_A != "KZ117237.1" & CHR_A != "KZ117131.1" & CHR_A != "KZ116099.1" &
                           CHR_A != "KZ117975.1" & CHR_A != "KZ118067.1" & CHR_A != "KZ117108.1")
head(LD_2017_notSel)


##### Getting distance between SNPs in kbps #####
#2002
#selected
LD_2002_in3Mb.1$distancekb <- with(LD_2002_in3Mb.1, LD_2002_in3Mb.1$BP_B-LD_2002_in3Mb.1$BP_A)/1000 ## the distance between snp1 and snp2 in kb
LD_2002_in3Mb.2$distancekb <- with(LD_2002_in3Mb.2, LD_2002_in3Mb.2$BP_B-LD_2002_in3Mb.2$BP_A)/1000
LD_2002_Cry1Chr7.1$distancekb <- with(LD_2002_Cry1Chr7.1, LD_2002_Cry1Chr7.1$BP_B-LD_2002_Cry1Chr7.1$BP_A)/1000
LD_2002_Cry1Chr7.2$distancekb <- with(LD_2002_Cry1Chr7.2, LD_2002_Cry1Chr7.2$BP_B-LD_2002_Cry1Chr7.2$BP_A)/1000 
LD_2002_Cry1Chr6$distancekb <- with(LD_2002_Cry1Chr6, LD_2002_Cry1Chr6$BP_B-LD_2002_Cry1Chr6$BP_A)/1000
#LD_2002_Cry1Chr11$distancekb <- with(LD_2002_Cry1Chr11, LD_2002_Cry1Chr11$BP_B-LD_2002_Cry1Chr11$BP_A)/1000

LD_2002_Cry2Chr7$distancekb <- with(LD_2002_Cry2Chr7, LD_2002_Cry2Chr7$BP_B-LD_2002_Cry2Chr7$BP_A)/1000 

#neutral
LD_2002_notSel$distancekb <- with(LD_2002_notSel, LD_2002_notSel$BP_B-LD_2002_notSel$BP_A)/1000 


#2012
#selected
LD_2012_in3Mb.1$distancekb <- with(LD_2012_in3Mb.1, LD_2012_in3Mb.1$BP_B-LD_2012_in3Mb.1$BP_A)/1000 ## the distance between snp1 and snp2 in kb
LD_2012_in3Mb.2$distancekb <- with(LD_2012_in3Mb.2, LD_2012_in3Mb.2$BP_B-LD_2012_in3Mb.2$BP_A)/1000
LD_2012_Cry1Chr7.1$distancekb <- with(LD_2012_Cry1Chr7.1, LD_2012_Cry1Chr7.1$BP_B-LD_2012_Cry1Chr7.1$BP_A)/1000
LD_2012_Cry1Chr7.2$distancekb <- with(LD_2012_Cry1Chr7.2, LD_2012_Cry1Chr7.2$BP_B-LD_2012_Cry1Chr7.2$BP_A)/1000 
LD_2012_Cry1Chr6$distancekb <- with(LD_2012_Cry1Chr6, LD_2012_Cry1Chr6$BP_B-LD_2012_Cry1Chr6$BP_A)/1000
#LD_2012_Cry1Chr11$distancekb <- with(LD_2012_Cry1Chr11, LD_2012_Cry1Chr11$BP_B-LD_2012_Cry1Chr11$BP_A)/1000

LD_2012_Cry2Chr7$distancekb <- with(LD_2012_Cry2Chr7, LD_2012_Cry2Chr7$BP_B-LD_2012_Cry2Chr7$BP_A)/1000 

#neutral
LD_2012_notSel$distancekb <- with(LD_2012_notSel, LD_2012_notSel$BP_B-LD_2012_notSel$BP_A)/1000 


#2017
#selected
LD_2017_in3Mb.1$distancekb <- with(LD_2017_in3Mb.1, LD_2017_in3Mb.1$BP_B-LD_2017_in3Mb.1$BP_A)/1000 ## the distance between snp1 and snp2 in kb
LD_2017_in3Mb.2$distancekb <- with(LD_2017_in3Mb.2, LD_2017_in3Mb.2$BP_B-LD_2017_in3Mb.2$BP_A)/1000
LD_2017_Cry1Chr7.1$distancekb <- with(LD_2017_Cry1Chr7.1, LD_2017_Cry1Chr7.1$BP_B-LD_2017_Cry1Chr7.1$BP_A)/1000
LD_2017_Cry1Chr7.2$distancekb <- with(LD_2017_Cry1Chr7.2, LD_2017_Cry1Chr7.2$BP_B-LD_2017_Cry1Chr7.2$BP_A)/1000 
LD_2017_Cry1Chr6$distancekb <- with(LD_2017_Cry1Chr6, LD_2017_Cry1Chr6$BP_B-LD_2017_Cry1Chr6$BP_A)/1000
#LD_2017_Cry1Chr11$distancekb <- with(LD_2017_Cry1Chr11, LD_2017_Cry1Chr11$BP_B-LD_2017_Cry1Chr11$BP_A)/1000

LD_2017_Cry2Chr7$distancekb <- with(LD_2017_Cry2Chr7, LD_2017_Cry2Chr7$BP_B-LD_2017_Cry2Chr7$BP_A)/1000 

#neutral
LD_2017_notSel$distancekb <- with(LD_2017_notSel, LD_2017_notSel$BP_B-LD_2017_notSel$BP_A)/1000 


##### Making 10kb bins #####
#2002
LD_2002_in3Mb.1$grp <- cut(LD_2002_in3Mb.1$distancekb, 0:10) ## bin 10kb
LD_2002_in3Mb.2$grp <- cut(LD_2002_in3Mb.2$distancekb, 0:10)
LD_2002_Cry1Chr6$grp <- cut(LD_2002_Cry1Chr6$distancekb, 0:10)
LD_2002_Cry1Chr7.1$grp <- cut(LD_2002_Cry1Chr7.1$distancekb, 0:10)
LD_2002_Cry1Chr7.2$grp <- cut(LD_2002_Cry1Chr7.2$distancekb, 0:10)
#LD_2002_Cry1Chr11$grp <- cut(LD_2002_Cry1Chr11$distancekb, 0:10)

LD_2002_Cry2Chr7$grp <- cut(LD_2002_Cry2Chr7$distancekb, 0:10)

LD_2002_notSel$grp <- cut(LD_2002_notSel$distancekb, 0:10)

#2012
LD_2012_in3Mb.1$grp <- cut(LD_2012_in3Mb.1$distancekb, 0:10) ## bin 10kb
LD_2012_in3Mb.2$grp <- cut(LD_2012_in3Mb.2$distancekb, 0:10)
LD_2012_Cry1Chr6$grp <- cut(LD_2012_Cry1Chr6$distancekb, 0:10)
LD_2012_Cry1Chr7.1$grp <- cut(LD_2012_Cry1Chr7.1$distancekb, 0:10)
LD_2012_Cry1Chr7.2$grp <- cut(LD_2012_Cry1Chr7.2$distancekb, 0:10)
#LD_2012_Cry1Chr11$grp <- cut(LD_2012_Cry1Chr11$distancekb, 0:10)

LD_2012_Cry2Chr7$grp <- cut(LD_2012_Cry2Chr7$distancekb, 0:10)

LD_2012_notSel$grp <- cut(LD_2012_notSel$distancekb, 0:10)


#2017
LD_2017_in3Mb.1$grp <- cut(LD_2017_in3Mb.1$distancekb, 0:10) ## bin 10kb
LD_2017_in3Mb.2$grp <- cut(LD_2017_in3Mb.2$distancekb, 0:10)
LD_2017_Cry1Chr6$grp <- cut(LD_2017_Cry1Chr6$distancekb, 0:10)
LD_2017_Cry1Chr7.1$grp <- cut(LD_2017_Cry1Chr7.1$distancekb, 0:10)
LD_2017_Cry1Chr7.2$grp <- cut(LD_2017_Cry1Chr7.2$distancekb, 0:10)
#LD_2017_Cry1Chr11$grp <- cut(LD_2017_Cry1Chr11$distancekb, 0:10)

LD_2017_Cry2Chr7$grp <- cut(LD_2017_Cry2Chr7$distancekb, 0:10)

LD_2017_notSel$grp <- cut(LD_2017_notSel$distancekb, 0:10)



#A quick spot check before R-squared
head(LD_2002_in3Mb.1)
head(LD_2002_notSel)

head(LD_2017_in3Mb.1)
head(LD_2017_notSel)

##getting mean R-squared values per kb
#2002
r2means2002_in3Mb.1 <- with(LD_2002_in3Mb.1, tapply(LD_2002_in3Mb.1$R2, LD_2002_in3Mb.1$grp, FUN = mean))
head(r2means2002_in3Mb.1) #sanity check

r2means2002_in3Mb.2 <- with(LD_2002_in3Mb.2, tapply(LD_2002_in3Mb.2$R2, LD_2002_in3Mb.2$grp, FUN = mean))
r2means2002_Cry1Chr6 <- with(LD_2002_Cry1Chr6, tapply(LD_2002_Cry1Chr6$R2, LD_2002_Cry1Chr6$grp, FUN = mean))
r2means2002_Cry1Chr7.1 <- with(LD_2002_Cry1Chr7.1, tapply(LD_2002_Cry1Chr7.1$R2, LD_2002_Cry1Chr7.1$grp, FUN = mean))
r2means2002_Cry1Chr7.2 <- with(LD_2002_Cry1Chr7.2, tapply(LD_2002_Cry1Chr7.2$R2, LD_2002_Cry1Chr7.2$grp, FUN = mean))
#r2means2002_Cry1Chr11 <- with(LD_2002_Cry1Chr11, tapply(LD_2002_Cry1Chr11$R2, LD_2002_Cry1Chr11$grp, FUN = mean))

r2means2002_Cry2Chr7 <- with(LD_2002_Cry2Chr7, tapply(LD_2002_Cry2Chr7$R2, LD_2002_Cry2Chr7$grp, FUN = mean))

r2means2002_notSel <- with(LD_2002_notSel, tapply(LD_2002_notSel$R2, LD_2002_notSel$grp, FUN = mean))
head(r2means2002_notSel) #sanity check


#2012
r2means2012_in3Mb.1 <- with(LD_2012_in3Mb.1, tapply(LD_2012_in3Mb.1$R2, LD_2012_in3Mb.1$grp, FUN = mean))
head(r2means2012_in3Mb.1) #sanity check

r2means2012_in3Mb.2 <- with(LD_2012_in3Mb.2, tapply(LD_2012_in3Mb.2$R2, LD_2012_in3Mb.2$grp, FUN = mean))
r2means2012_Cry1Chr6 <- with(LD_2012_Cry1Chr6, tapply(LD_2012_Cry1Chr6$R2, LD_2012_Cry1Chr6$grp, FUN = mean))
r2means2012_Cry1Chr7.1 <- with(LD_2012_Cry1Chr7.1, tapply(LD_2012_Cry1Chr7.1$R2, LD_2012_Cry1Chr7.1$grp, FUN = mean))
r2means2012_Cry1Chr7.2 <- with(LD_2012_Cry1Chr7.2, tapply(LD_2012_Cry1Chr7.2$R2, LD_2012_Cry1Chr7.2$grp, FUN = mean))
#r2means2012_Cry1Chr11 <- with(LD_2012_Cry1Chr11, tapply(LD_2012_Cry1Chr11$R2, LD_2012_Cry1Chr11$grp, FUN = mean))

r2means2012_Cry2Chr7 <- with(LD_2012_Cry2Chr7, tapply(LD_2012_Cry2Chr7$R2, LD_2012_Cry2Chr7$grp, FUN = mean))

r2means2012_notSel <- with(LD_2012_notSel, tapply(LD_2012_notSel$R2, LD_2012_notSel$grp, FUN = mean))
head(r2means2012_notSel) #sanity check


#2017
r2means2017_in3Mb.1 <- with(LD_2017_in3Mb.1, tapply(LD_2017_in3Mb.1$R2, LD_2017_in3Mb.1$grp, FUN = mean))
head(r2means2017_in3Mb.1) #sanity check

r2means2017_in3Mb.2 <- with(LD_2017_in3Mb.2, tapply(LD_2017_in3Mb.2$R2, LD_2017_in3Mb.2$grp, FUN = mean))
r2means2017_Cry1Chr6 <- with(LD_2017_Cry1Chr6, tapply(LD_2017_Cry1Chr6$R2, LD_2017_Cry1Chr6$grp, FUN = mean))
r2means2017_Cry1Chr7.1 <- with(LD_2017_Cry1Chr7.1, tapply(LD_2017_Cry1Chr7.1$R2, LD_2017_Cry1Chr7.1$grp, FUN = mean))
r2means2017_Cry1Chr7.2 <- with(LD_2017_Cry1Chr7.2, tapply(LD_2017_Cry1Chr7.2$R2, LD_2017_Cry1Chr7.2$grp, FUN = mean))
#r2means2017_Cry1Chr11 <- with(LD_2017_Cry1Chr11, tapply(LD_2017_Cry1Chr11$R2, LD_2017_Cry1Chr11$grp, FUN = mean))

r2means2017_Cry2Chr7 <- with(LD_2017_Cry2Chr7, tapply(LD_2017_Cry2Chr7$R2, LD_2017_Cry2Chr7$grp, FUN = mean))

r2means2017_notSel <- with(LD_2017_notSel, tapply(LD_2017_notSel$R2, LD_2017_notSel$grp, FUN = mean))
head(r2means2017_notSel) #sanity check



##### Figures showing LD decay by year ######

png(filename = "/media/megan/easystore/Fig4_2002_LD_decay.png", units = "px", height = 700, width = 800)
par(mar = c(5,6,4,2))
plot(r2means2002_notSel, ylab = "Mean R-squared", xlab = "Distance between SNPs (Kb)", ylim = c(0.08,0.4), col = "black", 
     type = "l", lty = 1, lwd = 4, cex.lab = 2, cex.axis = 1.8)
lines(r2means2002_in3Mb.1, type = "l", col = "red", lty = 1, lwd = 2.5)
lines(r2means2002_in3Mb.2, type = "l", col = "red", lty = 1, lwd = 2.5)
lines(r2means2002_Cry1Chr7.1, type = "l", col = "blue", lty = 1, lwd = 2.5)
lines(r2means2002_Cry1Chr7.2, type = "l", col = "blue", lty = 1, lwd = 2.5)
lines(r2means2002_Cry2Chr7, type = "l", col = "blue", lty = 1, lwd = 2.5)
lines(r2means2002_Cry1Chr6, type = "l", col = "blue", lty = 1, lwd = 2.5)
#lines(r2means2002_Cry1Chr11, type = "l", col = "blue", lty = 1, lwd = 2.5)



dev.off()

png(filename = "/media/megan/easystore/Fig4_2012_LD_decay.png", units = "px", height = 700, width = 800)
par(mar = c(5,6,4,2))
plot(r2means2012_notSel, ylab = "Mean R-squared", xlab = "Distance between SNPs (Kb)", ylim = c(0.08,0.4), col = "black", 
     type = "l", lty = 1, lwd = 4, cex.lab = 2, cex.axis = 1.8)
lines(r2means2012_in3Mb.1, type = "l", col = "red", lty = 1, lwd = 2.5)
lines(r2means2012_in3Mb.2, type = "l", col = "red", lty = 1, lwd = 2.5)
lines(r2means2012_Cry1Chr7.1, type = "l", col = "blue", lty = 1, lwd = 2.5)
lines(r2means2012_Cry1Chr7.2, type = "l", col = "blue", lty = 1, lwd = 2.5)
lines(r2means2012_Cry2Chr7, type = "l", col = "blue", lty = 1, lwd = 2.5)
lines(r2means2012_Cry1Chr6, type = "l", col = "blue", lty = 1, lwd = 2.5)
#lines(r2means2012_Cry1Chr11, type = "l", col = "blue", lty = 1, lwd = 2.5)


dev.off()

png(filename = "/media/megan/easystore/Fig4_2017_LD_decay.png", units = "px", height = 700, width = 800)
par(mar = c(5,6,4,2))
plot(r2means2017_notSel, ylab = "Mean R-squared", xlab = "Distance between SNPs (Kb)", ylim = c(0.08,0.4), col = "black", 
     type = "l", lty = 1, lwd = 4, cex.lab = 2, cex.axis = 1.8)
lines(r2means2017_in3Mb.1, type = "l", col = "red", lty = 1, lwd = 2.5)
lines(r2means2017_in3Mb.2, type = "l", col = "red", lty = 1, lwd = 2.5)
lines(r2means2017_Cry1Chr7.1, type = "l", col = "blue", lty = 1, lwd = 2.5)
lines(r2means2017_Cry1Chr7.2, type = "l", col = "blue", lty = 1, lwd = 2.5)
lines(r2means2017_Cry2Chr7, type = "l", col = "blue", lty = 1, lwd = 2.5)
lines(r2means2017_Cry1Chr6, type = "l", col = "blue", lty = 1, lwd = 2.5)
#lines(r2means2017_Cry1Chr11, type = "l", col = "blue", lty = 1, lwd = 2.5)


dev.off()


##### Figures showing LD decay by Chr ######

png(filename = "/media/megan/easystore/Fig4_Chr5_LD_decay.png", units = "px", height = 700, width = 800)
par(mar = c(5,6,4,2))
plot(r2means2002_notSel, ylab = "Mean R-squared", xlab = "Distance between SNPs (Kb)", ylim = c(0.08,0.4), col = "black", 
     type = "l", lty = 1, lwd = 2, cex.lab = 1.8, cex.axis = 1.5)
lines(r2means2002_in3Mb.1, type = "l", col = "red", lty = 1, lwd = 2)
lines(r2means2002_in3Mb.2, type = "l", col = "blue", lty = 1, lwd = 2)


lines(r2means2012_in3Mb.1, type = "l", col = "red", lty = 2, lwd = 2)
lines(r2means2012_in3Mb.2, type = "l", col = "blue", lty = 2, lwd = 2)
lines(r2means2012_notSel, type = "l", col = "black", lty = 2, lwd = 1.5)


lines(r2means2017_in3Mb.1, type = "l", col = "red", lty = 3, lwd = 2)
lines(r2means2017_in3Mb.2, type = "l", col = "blue", lty = 3, lwd = 2)
lines(r2means2017_notSel, type = "l", col = "black", lty = 3, lwd = 1.5)


dev.off()


png(filename = "/media/megan/easystore/Fig4_Cry1Chr7_LD_decay.png", units = "px", height = 700, width = 800)
par(mar = c(5,6,4,2))
plot(r2means2002_notSel, ylab = "Mean R-squared", xlab = "Distance between SNPs (Kb)", ylim = c(0.08,0.26), col = "black", 
     type = "l", lty = 1, lwd = 2, cex.lab = 1.8, cex.axis = 1.5)
lines(r2means2002_Cry1Chr7.1, type = "l", col = "red", lty = 1, lwd = 2)
lines(r2means2002_Cry1Chr7.2, type = "l", col = "blue", lty = 1, lwd = 2)

lines(r2means2012_Cry1Chr7.1, type = "l", col = "red", lty = 2, lwd = 2)
lines(r2means2012_Cry1Chr7.2, type = "l", col = "blue", lty = 2, lwd = 2)
lines(r2means2012_notSel, type = "l", col = "black", lty = 2, lwd = 1.5)

lines(r2means2017_Cry1Chr7.1, type = "l", col = "red", lty = 3, lwd = 2)
lines(r2means2017_Cry1Chr7.2, type = "l", col = "blue", lty = 3, lwd = 2)
lines(r2means2017_notSel, type = "l", col = "black", lty = 3, lwd = 1.5)

dev.off()


png(filename = "/media/megan/easystore/Fig4_Cry2Chr7_LD_decay.png", units = "px", height = 700, width = 800)
par(mar = c(5,6,4,2))
plot(r2means2002_notSel, ylab = "Mean R-squared", xlab = "Distance between SNPs (Kb)", ylim = c(0.08,0.4), col = "black", 
     type = "l", lty = 1, lwd = 2, cex.lab = 1.8, cex.axis = 1.5)
lines(r2means2002_Cry2Chr7, type = "l", col = "red", lty = 1, lwd = 2)

lines(r2means2012_Cry2Chr7, type = "l", col = "red", lty = 2, lwd = 2)
lines(r2means2012_notSel, type = "l", col = "black", lty = 2, lwd = 1.5)

lines(r2means2017_Cry2Chr7, type = "l", col = "red", lty = 3, lwd = 2)
lines(r2means2017_notSel, type = "l", col = "black", lty = 3, lwd = 1.5)

dev.off()



png(filename = "/media/megan/easystore/Fig4_Cry1Chr6_LD_decay.png", units = "px", height = 700, width = 800)
par(mar = c(5,6,4,2))
plot(r2means2002_notSel, ylab = "Mean R-squared", xlab = "Distance between SNPs (Kb)", ylim = c(0.08,0.26), col = "black", 
     type = "l", lty = 1, lwd = 2, cex.lab = 1.8, cex.axis = 1.5)
lines(r2means2002_Cry1Chr6, type = "l", col = "red", lty = 1, lwd = 2)

lines(r2means2012_Cry1Chr6, type = "l", col = "red", lty = 2, lwd = 2)
lines(r2means2012_notSel, type = "l", col = "black", lty = 2, lwd = 1.5)

lines(r2means2017_Cry1Chr6, type = "l", col = "red", lty = 3, lwd = 2)
lines(r2means2017_notSel, type = "l", col = "black", lty = 3, lwd = 1.5)

dev.off()


##### Comparing LD Decay in whole genome and average within the scaffold with Window #####

Cry1_KZ118067.1_2002 <- subset(LD_2002_Cry1Chr6, BP_A > (144000) & BP_A < (146000+1000))
Cry1_KZ118067.1_2012 <- subset(LD_2012_Cry1Chr6, BP_A > (144000) & BP_A < (146000+1000))
Cry1_KZ118067.1_2017 <- subset(LD_2017_Cry1Chr6, BP_A > (144000) & BP_A < (146000+1000))
r2meansCry1_KZ118067.1_2002 <- with(Cry1_KZ118067.1_2002, tapply(Cry1_KZ118067.1_2002$R2, Cry1_KZ118067.1_2002$grp, FUN = mean))
r2meansCry1_KZ118067.1_2012 <- with(Cry1_KZ118067.1_2012, tapply(Cry1_KZ118067.1_2012$R2, Cry1_KZ118067.1_2012$grp, FUN = mean))
r2meansCry1_KZ118067.1_2017 <- with(Cry1_KZ118067.1_2017, tapply(Cry1_KZ118067.1_2017$R2, Cry1_KZ118067.1_2017$grp, FUN = mean))


Cry1_KZ116099.1_2002 <- subset(LD_2002_Cry1Chr7.1, BP_A > (2000) & BP_A < (7000+1000))
Cry1_KZ116099.1_2012 <- subset(LD_2012_Cry1Chr7.1, BP_A > (2000) & BP_A < (7000+1000))
Cry1_KZ116099.1_2017 <- subset(LD_2017_Cry1Chr7.1, BP_A > (2000) & BP_A < (7000+1000))
r2meansCry1_KZ116099.1_2002 <- with(Cry1_KZ116099.1_2002, tapply(Cry1_KZ116099.1_2002$R2, Cry1_KZ116099.1_2002$grp, FUN = mean))
r2meansCry1_KZ116099.1_2012 <- with(Cry1_KZ116099.1_2012, tapply(Cry1_KZ116099.1_2012$R2, Cry1_KZ116099.1_2012$grp, FUN = mean))
r2meansCry1_KZ116099.1_2017 <- with(Cry1_KZ116099.1_2017, tapply(Cry1_KZ116099.1_2017$R2, Cry1_KZ116099.1_2017$grp, FUN = mean))


Cry1_KZ117975.1_2002 <- subset(LD_2002_Cry1Chr7.2, BP_A > (43000) & BP_A < (46000+1000))
Cry1_KZ117975.1_2012 <- subset(LD_2012_Cry1Chr7.2, BP_A > (43000) & BP_A < (46000+1000))
Cry1_KZ117975.1_2017 <- subset(LD_2017_Cry1Chr7.2, BP_A > (43000) & BP_A < (46000+1000))
r2meansCry1_KZ117975.1_2002 <- with(Cry1_KZ117975.1_2002, tapply(Cry1_KZ117975.1_2002$R2, Cry1_KZ117975.1_2002$grp, FUN = mean))
r2meansCry1_KZ117975.1_2012 <- with(Cry1_KZ117975.1_2012, tapply(Cry1_KZ117975.1_2012$R2, Cry1_KZ117975.1_2012$grp, FUN = mean))
r2meansCry1_KZ117975.1_2017 <- with(Cry1_KZ117975.1_2017, tapply(Cry1_KZ117975.1_2017$R2, Cry1_KZ117975.1_2017$grp, FUN = mean))


Cry1_KZ118133.1_2002 <- subset(LD_2002_Cry1Chr7.2, BP_A > (34000) & BP_A < (35000+1000))
Cry1_KZ118133.1_2012 <- subset(LD_2012_Cry1Chr7.2, BP_A > (34000) & BP_A < (35000+1000))
Cry1_KZ118133.1_2017 <- subset(LD_2017_Cry1Chr7.2, BP_A > (34000) & BP_A < (35000+1000))
r2meansCry1_KZ118133.1_2002 <- with(Cry1_KZ118133.1_2002, tapply(Cry1_KZ118133.1_2002$R2, Cry1_KZ118133.1_2002$grp, FUN = mean))
r2meansCry1_KZ118133.1_2012 <- with(Cry1_KZ118133.1_2012, tapply(Cry1_KZ118133.1_2012$R2, Cry1_KZ118133.1_2012$grp, FUN = mean))
r2meansCry1_KZ118133.1_2017 <- with(Cry1_KZ118133.1_2017, tapply(Cry1_KZ118133.1_2017$R2, Cry1_KZ118133.1_2017$grp, FUN = mean))


