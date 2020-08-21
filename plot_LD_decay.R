#plotting ld decay

library(ggplot2)

# import the data
LD_2002 <- data.frame(read.table("/media/megan/easystore/Hz2002_plink.ld", header=T))
LD_2017 <- data.frame(read.table("/media/megan/easystore/Hz2017_plink.ld", header=T))

#subsetting data to include scaffolds within and outside of 3Mb region for 2002
LD_2002_in3Mb <- subset(LD_2002, CHR_A == "KZ118395.1" | CHR_A == "KZ118483.1" | CHR_A == "KZ117237.1" | CHR_A == "KZ117131.1")
head(LD_2002_in3Mb)

LD_2002_not3Mb <- subset(LD_2002, CHR_A != "KZ118395.1" & CHR_A != "KZ118483.1" & CHR_A != "KZ117237.1" & CHR_A != "KZ117131.1")
head(LD_2002_in3Mb)

#pulling out top 6 scaffolds over 1 Mb for 2002
LD_2002_scaf0 <- subset(LD_2002, CHR_A == "KZ116521.1")
head(LD_2002_scaf0)
LD_2002_scaf1 <- subset(LD_2002, CHR_A == "KZ116522.1")
head(LD_2002_scaf1)
LD_2002_scaf2 <- subset(LD_2002, CHR_A == "KZ117462.1")
head(LD_2002_scaf2)
LD_2002_scaf3 <- subset(LD_2002, CHR_A == "KZ118089.1")
head(LD_2002_scaf3)
LD_2002_scaf4 <- subset(LD_2002, CHR_A == "KZ118207.1")
head(LD_2002_scaf4)
LD_2002_scaf5 <- subset(LD_2002, CHR_A == "KZ118318.1")
head(LD_2002_scaf5)


#subsetting data to include scaffolds within and outside of 3Mb region for 2017
LD_2017_in3Mb <- subset(LD_2017, CHR_A == "KZ118395.1" | CHR_A == "KZ118483.1" | CHR_A == "KZ117237.1" | CHR_A == "KZ117131.1")
head(LD_2017_in3Mb)

LD_2017_not3Mb <- subset(LD_2017, CHR_A != "KZ118395.1" & CHR_A != "KZ118483.1" & CHR_A != "KZ117237.1" & CHR_A != "KZ117131.1")
head(LD_2017_in3Mb)

#pulling out top 6 scaffolds over 1 Mb for 2017
LD_2017_scaf0 <- subset(LD_2017, CHR_A == "KZ116521.1")
head(LD_2017_scaf0)
LD_2017_scaf1 <- subset(LD_2017, CHR_A == "KZ116522.1")
head(LD_2017_scaf1)
LD_2017_scaf2 <- subset(LD_2017, CHR_A == "KZ117462.1")
head(LD_2017_scaf2)
LD_2017_scaf3 <- subset(LD_2017, CHR_A == "KZ118089.1")
head(LD_2017_scaf3)
LD_2017_scaf4 <- subset(LD_2017, CHR_A == "KZ118207.1")
head(LD_2017_scaf4)
LD_2017_scaf5 <- subset(LD_2017, CHR_A == "KZ118318.1")
head(LD_2017_scaf5)

###getting distance between SNPs in kbps
#2002
LD_2002_in3Mb$distancekb <- with(LD_2002_in3Mb, LD_2002_in3Mb$BP_B-LD_2002_in3Mb$BP_A)/1000 ## the distance between snp1 and snp2 in kb
LD_2002_not3Mb$distancekb <- with(LD_2002_not3Mb, LD_2002_not3Mb$BP_B-LD_2002_not3Mb$BP_A)/1000 

LD_2002_scaf0$distancekb <- with(LD_2002_scaf0, LD_2002_scaf0$BP_B - LD_2002_scaf0$BP_A)/1000
LD_2002_scaf1$distancekb <- with(LD_2002_scaf1, LD_2002_scaf1$BP_B - LD_2002_scaf1$BP_A)/1000
LD_2002_scaf2$distancekb <- with(LD_2002_scaf2, LD_2002_scaf2$BP_B - LD_2002_scaf2$BP_A)/1000
LD_2002_scaf3$distancekb <- with(LD_2002_scaf3, LD_2002_scaf3$BP_B - LD_2002_scaf3$BP_A)/1000
LD_2002_scaf4$distancekb <- with(LD_2002_scaf4, LD_2002_scaf4$BP_B - LD_2002_scaf4$BP_A)/1000
LD_2002_scaf5$distancekb <- with(LD_2002_scaf5, LD_2002_scaf5$BP_B - LD_2002_scaf5$BP_A)/1000

#2017
LD_2017_in3Mb$distancekb <- with(LD_2017_in3Mb, LD_2017_in3Mb$BP_B-LD_2017_in3Mb$BP_A)/1000 ## the distance between snp1 and snp2 in kb
LD_2017_not3Mb$distancekb <- with(LD_2017_not3Mb, LD_2017_not3Mb$BP_B-LD_2017_not3Mb$BP_A)/1000 

LD_2017_scaf0$distancekb <- with(LD_2017_scaf0, LD_2017_scaf0$BP_B - LD_2017_scaf0$BP_A)/1000
LD_2017_scaf1$distancekb <- with(LD_2017_scaf1, LD_2017_scaf1$BP_B - LD_2017_scaf1$BP_A)/1000
LD_2017_scaf2$distancekb <- with(LD_2017_scaf2, LD_2017_scaf2$BP_B - LD_2017_scaf2$BP_A)/1000
LD_2017_scaf3$distancekb <- with(LD_2017_scaf3, LD_2017_scaf3$BP_B - LD_2017_scaf3$BP_A)/1000
LD_2017_scaf4$distancekb <- with(LD_2017_scaf4, LD_2017_scaf4$BP_B - LD_2017_scaf4$BP_A)/1000
LD_2017_scaf5$distancekb <- with(LD_2017_scaf5, LD_2017_scaf5$BP_B - LD_2017_scaf5$BP_A)/1000


##making 40kb bins
#2002
LD_2002_in3Mb$grp <- cut(LD_2002_in3Mb$distancekb, 0:40) ## bin 40kb
LD_2002_not3Mb$grp <- cut(LD_2002_not3Mb$distancekb, 0:40)
LD_2002_scaf0$grp <- cut(LD_2002_scaf0$distancekb, 0:40)
LD_2002_scaf1$grp <- cut(LD_2002_scaf1$distancekb, 0:40)
LD_2002_scaf2$grp <- cut(LD_2002_scaf2$distancekb, 0:40)
LD_2002_scaf3$grp <- cut(LD_2002_scaf3$distancekb, 0:40)
LD_2002_scaf4$grp <- cut(LD_2002_scaf4$distancekb, 0:40)
LD_2002_scaf5$grp <- cut(LD_2002_scaf5$distancekb, 0:40)

#2017
LD_2017_in3Mb$grp <- cut(LD_2017_in3Mb$distancekb, 0:40) ## bin 40kb
LD_2017_not3Mb$grp <- cut(LD_2017_not3Mb$distancekb, 0:40) ## bin 40kb
LD_2017_scaf0$grp <- cut(LD_2017_scaf0$distancekb, 0:40)
LD_2017_scaf1$grp <- cut(LD_2017_scaf1$distancekb, 0:40)
LD_2017_scaf2$grp <- cut(LD_2017_scaf2$distancekb, 0:40)
LD_2017_scaf3$grp <- cut(LD_2017_scaf3$distancekb, 0:40)
LD_2017_scaf4$grp <- cut(LD_2017_scaf4$distancekb, 0:40)
LD_2017_scaf5$grp <- cut(LD_2017_scaf5$distancekb, 0:40)

#A quick check before R-squared
head(LD_2002_in3Mb)
head(LD_2002_not3Mb)

head(LD_2017_in3Mb)
head(LD_2017_not3Mb)

##getting mean R-squared values per kb
#2002
r2means2002_in3Mb <- with(LD_2002_in3Mb, tapply(LD_2002_in3Mb$R2, LD_2002_in3Mb$grp, FUN = mean))
head(r2means2002_in3Mb) #sanity check

r2means2002_not3Mb <- with(LD_2002_not3Mb, tapply(LD_2002_not3Mb$R2, LD_2002_not3Mb$grp, FUN = mean))
head(r2means2002_not3Mb) #sanity check

r2means2002_scaf0 <-with(LD_2002_scaf0, tapply(LD_2002_scaf0$R2, LD_2002_scaf0$grp, FUN = mean))
r2means2002_scaf1 <-with(LD_2002_scaf1, tapply(LD_2002_scaf1$R2, LD_2002_scaf1$grp, FUN = mean))
r2means2002_scaf2 <-with(LD_2002_scaf2, tapply(LD_2002_scaf2$R2, LD_2002_scaf2$grp, FUN = mean))
r2means2002_scaf3 <-with(LD_2002_scaf3, tapply(LD_2002_scaf3$R2, LD_2002_scaf3$grp, FUN = mean))
r2means2002_scaf4 <-with(LD_2002_scaf4, tapply(LD_2002_scaf4$R2, LD_2002_scaf4$grp, FUN = mean))
r2means2002_scaf5 <-with(LD_2002_scaf5, tapply(LD_2002_scaf5$R2, LD_2002_scaf5$grp, FUN = mean))


#2017
r2means2017_in3Mb <- with(LD_2017_in3Mb, tapply(LD_2017_in3Mb$R2, LD_2017_in3Mb$grp, FUN = mean))
head(r2means2017_in3Mb) #sanity check

r2means2017_not3Mb <- with(LD_2017_not3Mb, tapply(LD_2017_not3Mb$R2, LD_2017_not3Mb$grp, FUN = mean))
head(r2means2017_not3Mb) #sanity check

r2means2017_scaf0 <-with(LD_2017_scaf0, tapply(LD_2017_scaf0$R2, LD_2017_scaf0$grp, FUN = mean))
r2means2017_scaf1 <-with(LD_2017_scaf1, tapply(LD_2017_scaf1$R2, LD_2017_scaf1$grp, FUN = mean))
r2means2017_scaf2 <-with(LD_2017_scaf2, tapply(LD_2017_scaf2$R2, LD_2017_scaf2$grp, FUN = mean))
r2means2017_scaf3 <-with(LD_2017_scaf3, tapply(LD_2017_scaf3$R2, LD_2017_scaf3$grp, FUN = mean))
r2means2017_scaf4 <-with(LD_2017_scaf4, tapply(LD_2017_scaf4$R2, LD_2017_scaf4$grp, FUN = mean))
r2means2017_scaf5 <-with(LD_2017_scaf5, tapply(LD_2017_scaf5$R2, LD_2017_scaf5$grp, FUN = mean))


png(filename = "/media/megan/easystore/Fig4_LD_decay.png", units = "px", height = 700, width = 800)
par(mar = c(5,6,4,2))
plot(r2means2002_not3Mb, ylab = "Mean R-squared", xlab = "Distance between SNPs (Kb)", ylim = c(0.05,0.23), col = "black", 
     type = "l", lty = 1, lwd = 2, cex.lab = 1.8, cex.axis = 1.5)
lines(r2means2002_in3Mb, type = "l", col = "red", lty = 1, lwd = 2)
lines(r2means2002_scaf0, type = "l", col = "grey", lty = 1, lwd = 1.5)
lines(r2means2002_scaf1, type = "l", col = "grey", lty = 1, lwd = 1.5)
lines(r2means2002_scaf2, type = "l", col = "grey", lty = 1, lwd = 1.5)
lines(r2means2002_scaf3, type = "l", col = "grey", lty = 1, lwd = 1.5)
lines(r2means2002_scaf4, type = "l", col = "grey", lty = 1, lwd = 1.5)
lines(r2means2002_scaf5, type = "l", col = "grey", lty = 1, lwd = 1.5)

lines(r2means2017_not3Mb, type = "l", col = "black", lty = 2, lwd = 2)
lines(r2means2017_in3Mb, type = "l", col = "red", lty = 2, lwd = 2)
lines(r2means2017_scaf0, type = "l", col = "grey", lty = 2, lwd = 1.5)
lines(r2means2017_scaf1, type = "l", col = "grey", lty = 2, lwd = 1.5)
lines(r2means2017_scaf2, type = "l", col = "grey", lty = 2, lwd = 1.5)
lines(r2means2017_scaf3, type = "l", col = "grey", lty = 2, lwd = 1.5)
lines(r2means2017_scaf4, type = "l", col = "grey", lty = 2, lwd = 1.5)
lines(r2means2017_scaf5, type = "l", col = "grey", lty = 2, lwd = 1.5)
legend(27, 0.23, legend=c("2002 without 3Mb sweep", "2002 in 3Mb sweep", "2002 other scafs > 1Mb",
                        "2017 without 3Mb sweep", "2017 in 3Mb sweep", "2017 other scafs > 1Mb"),
       col=c("black", "red", "grey"), lty= c(1,1,1,2,2,2), cex=1.2)

dev.off()

