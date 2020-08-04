#plotting ld decay

library(ggplot2)

# import the data
LD_2002 <- data.frame(read.table("/media/megan/easystore/Hz2002_plink.ld", header=T))
LD_2017 <- data.frame(read.table("/media/megan/easystore/Hz2017_plink.ld", header=T))

#subsetting 2017 data to include scaffolds within and outside of 3Mb region

LD_2017_in3Mb <- subset(LD_2017, CHR_A == "KZ118395.1" | CHR_A == "KZ118483.1" | CHR_A == "KZ117237.1" | CHR_A == "KZ117131.1")
head(LD_2017_in3Mb)

LD_2017_not3Mb <- subset(LD_2017, CHR_A != "KZ118395.1" & CHR_A != "KZ118483.1" & CHR_A != "KZ117237.1" & CHR_A != "KZ117131.1")
head(LD_2017_in3Mb)

LD_2002$distancekb <- with(LD_2002, LD_2002$BP_B-LD_2002$BP_A)/1000 ## the distance between snp1 and snp2 in kb
LD_2017_in3Mb$distancekb <- with(LD_2017_in3Mb, LD_2017_in3Mb$BP_B-LD_2017_in3Mb$BP_A)/1000 ## the distance between snp1 and snp2 in kb
LD_2017_not3Mb$distancekb <- with(LD_2017_not3Mb, LD_2017_not3Mb$BP_B-LD_2017_not3Mb$BP_A)/1000 ## the distance between snp1 and snp2 in kb


LD_2002$grp <- cut(LD_2002$distancekb, 0:40) ## bin 40kb
LD_2017_in3Mb$grp <- cut(LD_2017_in3Mb$distancekb, 0:40) ## bin 40kb
LD_2017_not3Mb$grp <- cut(LD_2017_not3Mb$distancekb, 0:40) ## bin 40kb

head(LD_2002)
head(LD_2017_in3Mb)
head(LD_2017_not3Mb)

#getting mean R2 values per kb
r2means2002 <- with(LD_2002, tapply(LD_2002$R2, LD_2002$grp, FUN = mean)) ##r2 mean every 1kb
head(r2means2002)

r2means2017_in3Mb <- with(LD_2017_in3Mb, tapply(LD_2017_in3Mb$R2, LD_2017_in3Mb$grp, FUN = mean))
head(r2means2017_in3Mb)

r2means2017_not3Mb <- with(LD_2017_not3Mb, tapply(LD_2017_not3Mb$R2, LD_2017_not3Mb$grp, FUN = mean))
head(r2means2017_not3Mb)

plot(r2means2002, ylab = "Mean R-squared", xlab = "Distance between SNPs (Kb)", ylim = c(0,0.25), col = "black")
lines(r2means2017_in3Mb, type = "p", col = "red")
lines(r2means2017_not3Mb, type = "p", col = "blue")

