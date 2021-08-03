#getting changes in haplotype frequencies

library(vcfR); library(haplo.stats); library(tidyr); library(plyr)

setwd("/media/fritzlab/EE9C16C89C168AEB/Reanalysis_PNAS")

##### Chr7 #####
#this became LG9

#KZ116099.1 - loading datasets
KZ116099.1_pos <- read.csv("Cry1_KZ116099.1_pos.csv", header = T)

KZ116099.1 <- read.vcfR("./KZ116099.1_only/KZ116099.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf")

KZ116099.1_loc <- vcfR2loci(KZ116099.1, return.alleles = TRUE)
str(KZ116099.1_loc)

#fixing colnames and subsetting by column
names(KZ116099.1_loc) <- gsub(x = names(KZ116099.1_loc), pattern = "KZ116099.1_", replacement = "")

importLoc <- as.character(names(KZ116099.1_loc))

idx <- match(as.character(KZ116099.1_pos$Pos), importLoc)
idx <- sort(idx)

New_KZ116099.1_loc <- KZ116099.1_loc[,idx] 

#too many NAs after the first SNPs so filtering the latter SNPs out.
New_KZ116099.1_loc <- New_KZ116099.1_loc[,c(1:5)]
names <- as.numeric(colnames(New_KZ116099.1_loc))


#splitting alleles into separate columns
KZ116099.1_spl <- data.frame(rownames(New_KZ116099.1_loc))
names(KZ116099.1_spl) <- "sample"

for (i in 1:ncol(New_KZ116099.1_loc)){
  vec <- as.vector(New_KZ116099.1_loc[,i])
  temp <- revalue(vec, c("."="NA/NA"))
  temp <- data.frame(temp)
  temp$temp <- as.character(temp$temp)
  print(head(temp))
  name1 <- paste0(i,".","1")
  name2 <- paste0(i,".","2")
  alleles <- temp %>% separate(col = "temp", into = c(name1,name2), sep = "/")
  print(head(alleles))
  KZ116099.1_spl <- cbind(KZ116099.1_spl,alleles)
}



#2002
KZ116099.1_spl.desc <- summaryGeno(KZ116099.1_spl[1:13,-1], miss.val=c(0,NA))
save.em.2002.KZ116099.1 <- haplo.em(geno=KZ116099.1_spl[1:13,-1], locus.label=names, miss.val=c(0,NA), 
                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2002.KZ116099.1)
print(save.em.2002.KZ116099.1, nlines = 12)

#2012
KZ116099.1_spl.desc <- summaryGeno(KZ116099.1_spl[14:24,-1], miss.val=c(0,NA))
save.em.2012.KZ116099.1 <- haplo.em(geno=KZ116099.1_spl[14:24,-1], locus.label=names, miss.val=c(0,NA), 
                         control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2012.KZ116099.1)
print(save.em.2012.KZ116099.1, nlines = 12)

#2017
KZ116099.1_spl.desc <- summaryGeno(KZ116099.1_spl[25:35,-1], miss.val=c(0,NA))
save.em.2017.KZ116099.1 <- haplo.em(geno=KZ116099.1_spl[24:35,-1], locus.label=names, miss.val=c(0,NA), 
                         control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2017.KZ116099.1)
print(save.em.2017.KZ116099.1, nlines = 12)


#KZ117975.1 - loading datasets
KZ117975.1_pos <- read.csv("Cry1_KZ117975.1_pos.csv", header = T)
KZ117975.1 <- read.vcfR("./KZ117975.1_only/KZ117975.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf")
KZ117975.1_loc <- vcfR2loci(KZ117975.1, return.alleles = TRUE)
str(KZ117975.1_loc)

#fixing colnames and subsetting by column
names(KZ117975.1_loc) <- gsub(x = names(KZ117975.1_loc), pattern = "KZ117975.1_", replacement = "")

importLoc <- as.character(names(KZ117975.1_loc))

idx <- match(as.character(KZ117975.1_pos$Pos), importLoc)
idx <- sort(idx)

New_KZ117975.1_loc <- KZ117975.1_loc[,idx] 

#too many NAs in last bp so filtering.
New_KZ117975.1_loc <- New_KZ117975.1_loc[,-7]
names <- as.numeric(colnames(New_KZ117975.1_loc))

#splitting alleles into separate columns
KZ117975.1_spl <- data.frame(rownames(New_KZ117975.1_loc))
names(KZ117975.1_spl) <- "sample"

for (i in 1:ncol(New_KZ117975.1_loc)){
  vec <- as.vector(New_KZ117975.1_loc[,i])
  temp <- revalue(vec, c("."="NA/NA"))
  temp <- data.frame(temp)
  temp$temp <- as.character(temp$temp)
  print(head(temp))
  name1 <- paste0(i,".","1")
  name2 <- paste0(i,".","2")
  alleles <- temp %>% separate(col = "temp", into = c(name1,name2), sep = "/")
  print(head(alleles))
  KZ117975.1_spl <- cbind(KZ117975.1_spl,alleles)
}


#2002
KZ117975.1_spl.desc <- summaryGeno(KZ117975.1_spl[1:13,-1], miss.val=c(0,NA))
save.em.2002.KZ117975.1 <- haplo.em(geno=KZ117975.1_spl[1:13,-1], locus.label=names, miss.val=c(0,NA), 
                         control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2002.KZ117975.1)
print(save.em.2002.KZ117975.1, nlines = 12)

#2012
KZ117975.1_spl.desc <- summaryGeno(KZ117975.1_spl[14:24,-1], miss.val=c(0,NA))
save.em.2012.KZ117975.1 <- haplo.em(geno=KZ117975.1_spl[14:24,-1], locus.label=names, miss.val=c(0,NA), 
                         control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2012.KZ117975.1)
print(save.em.2012.KZ117975.1, nlines = 12)

#2017
KZ117975.1_spl.desc <- summaryGeno(KZ117975.1_spl[25:35,-1], miss.val=c(0,NA))
save.em.2017.KZ117975.1 <- haplo.em(geno=KZ117975.1_spl[24:35,-1], locus.label=names, miss.val=c(0,NA), 
                         control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2017.KZ117975.1)
print(save.em.2017.KZ117975.1, nlines = 12)


#KZ117108.1 - loading datasets THIS IS CRY2!!!!
KZ117108.1_pos <- read.csv("Cry2_KZ117108.1_pos.csv", header = T)
max(KZ117108.1_pos$Pos)
min(KZ117108.1_pos$Pos)
KZ117108.1 <- read.vcfR("./KZ117108.1_only/KZ117108.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf")
KZ117108.1_loc <- vcfR2loci(KZ117108.1, return.alleles = TRUE)
str(KZ117108.1_loc)

KZ117108.1_loc2 <- Filter(function(x) !any(x=="."), KZ117108.1_loc)#removing some regions with a lot of NAs


#fixing colnames and subsetting by column
names(KZ117108.1_loc2) <- gsub(x = names(KZ117108.1_loc2), pattern = "KZ117108.1_", replacement = "")

importLoc <- as.character(names(KZ117108.1_loc2))

idx <- match(as.character(KZ117108.1_pos$Pos), importLoc)
idx <- sort(idx)

New_KZ117108.1_loc <- KZ117108.1_loc2[,idx] 

#splitting alleles into separate columns and removing columns with NAs
KZ117108.1_spl <- data.frame(rownames(New_KZ117108.1_loc))
names(KZ117108.1_spl) <- "sample"


for (i in 1:ncol(New_KZ117108.1_loc)){
  vec <- as.vector(New_KZ117108.1_loc[,i])
  temp <- revalue(vec, c("."="NA/NA"))
  temp <- data.frame(temp)
  temp$temp <- as.character(temp$temp)
  print(head(temp))
  name1 <- paste0(i,".","1")
  name2 <- paste0(i,".","2")
  alleles <- temp %>% separate(col = "temp", into = c(name1,name2), sep = "/")
  print(head(alleles))
  KZ117108.1_spl <- cbind(KZ117108.1_spl,alleles)
}




#For this locus, I removed a region with a lot of missing genotype calls, which impacted haplotypes calling
#2002
KZ117108.1_spl.desc <- summaryGeno(KZ117108.1_spl[1:13,-1], miss.val=c(0,NA))
save.em.2002.KZ117108.1 <- haplo.em(geno=KZ117108.1_spl[1:13,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2002.KZ117108.1)
print(save.em.2002.KZ117108.1, nlines = 12)

#2012
KZ117108.1_spl.desc <- summaryGeno(KZ117108.1_spl[14:24,-1], miss.val=c(0,NA))
save.em.2012.KZ117108.1 <- haplo.em(geno=KZ117108.1_spl[14:24,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2012.KZ117108.1)
print(save.em.2012.KZ117108.1, nlines = 12)

#2017
KZ117108.1_spl.desc <- summaryGeno(KZ117108.1_spl[25:35,-1], miss.val=c(0,NA))
save.em.2017.KZ117108.1 <- haplo.em(geno=KZ117108.1_spl[24:35,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2017.KZ117108.1)
print(save.em.2017.KZ117108.1, nlines = 12)



###### Chr6 ######
#this became LG18

#KZ118067.1 - loading datasets
KZ118067.1_pos <- read.csv("Cry1_KZ118067.1_pos.csv", header = T)
KZ118067.1 <- read.vcfR("./KZ118067.1_only/KZ118067.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf")
KZ118067.1_loc <- vcfR2loci(KZ118067.1, return.alleles = TRUE)
str(KZ118067.1_loc)

#fixing colnames and subsetting by column
names(KZ118067.1_loc) <- gsub(x = names(KZ118067.1_loc), pattern = "KZ118067.1_", replacement = "")

importLoc <- as.character(names(KZ118067.1_loc))

idx <- match(as.character(KZ118067.1_pos$Pos), importLoc)
idx <- sort(idx)

New_KZ118067.1_loc <- KZ118067.1_loc[,idx] 

#splitting alleles into separate columns
KZ118067.1_spl <- data.frame(rownames(New_KZ118067.1_loc))
names(KZ118067.1_spl) <- "sample"

for (i in 1:ncol(New_KZ118067.1_loc)){
  vec <- as.vector(New_KZ118067.1_loc[,i])
  temp <- revalue(vec, c("."="NA/NA"))
  temp <- data.frame(temp)
  temp$temp <- as.character(temp$temp)
  print(head(temp))
  name1 <- paste0(i,".","1")
  name2 <- paste0(i,".","2")
  alleles <- temp %>% separate(col = "temp", into = c(name1,name2), sep = "/")
  print(head(alleles))
  KZ118067.1_spl <- cbind(KZ118067.1_spl,alleles)
}


#2002
KZ118067.1_spl.desc <- summaryGeno(KZ118067.1_spl[1:13,-1], miss.val=c(0,NA))
save.em.2002.KZ118067.1 <- haplo.em(geno=KZ118067.1_spl[1:13,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2002.KZ118067.1)
print(save.em.2002.KZ118067.1, nlines = 12)

#2012
KZ118067.1_spl.desc <- summaryGeno(KZ118067.1_spl[14:24,-1], miss.val=c(0,NA))
save.em.2012.KZ118067.1 <- haplo.em(geno=KZ118067.1_spl[14:24,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2012.KZ118067.1)
print(save.em.2012.KZ118067.1, nlines = 12)

#2017
KZ118067.1_spl.desc <- summaryGeno(KZ118067.1_spl[25:35,-1], miss.val=c(0,NA))
save.em.2017.KZ118067.1 <- haplo.em(geno=KZ118067.1_spl[24:35,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2017.KZ118067.1)
print(save.em.2017.KZ118067.1, nlines = 12)


##### Chr5 #####
#this became LG13

KZ118395.1 <- read.vcfR("./KZ118395.1_only/KZ118395.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf")
KZ118395.1_loc <- vcfR2loci(KZ118395.1, return.alleles = TRUE)
str(KZ118395.1_loc)

#fixing colnames and subsetting by column
names(KZ118395.1_loc) <- gsub(x = names(KZ118395.1_loc), pattern = "KZ118395.1_", replacement = "")

importLoc <- as.numeric(as.character(names(KZ118395.1_loc)))
sub_importLoc <- importLoc[importLoc > 5000 & importLoc < 15000 & importLoc!=6990 & importLoc!=7549 & importLoc!=8700 & importLoc!=11488 & importLoc !=14317] 
#just picking one 10Kb region because whole scaffold looks to have sweep, plus had to remove 5 SNPs 
#in 10kb window that created ambiguous haplotypes b/c of missing data.

idx <- as.character(sub_importLoc)

New_KZ118395.1_loc <- data.frame(KZ118395.1_loc[,idx])

#splitting alleles into separate columns
KZ118395.1_spl <- data.frame(rownames(KZ118395.1_loc))
names(KZ118395.1_spl) <- "sample"

for (i in 1:ncol(New_KZ118395.1_loc)){
  vec <- as.vector(New_KZ118395.1_loc[,i])
  temp <- revalue(vec, c("."="NA/NA"))
  temp <- data.frame(temp)
  temp$temp <- as.character(temp$temp)
  print(head(temp))
  name1 <- paste0(i,".","1")
  name2 <- paste0(i,".","2")
  alleles <- temp %>% separate(col = "temp", into = c(name1,name2), sep = "/")
  print(head(alleles))
  KZ118395.1_spl <- cbind(KZ118395.1_spl,alleles)
}


#2002
KZ118395.1_spl.desc <- summaryGeno(KZ118395.1_spl[1:13,-1], miss.val=c(0,NA))
save.em.2002.KZ118395.1 <- haplo.em(geno=KZ118395.1_spl[1:13,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2002.KZ118395.1)
print(save.em.2002.KZ118395.1, nlines = 12)

#2012
KZ118395.1_spl.desc <- summaryGeno(KZ118395.1_spl[14:24,-1], miss.val=c(0,NA))
save.em.2012.KZ118395.1 <- haplo.em(geno=KZ118395.1_spl[14:24,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2012.KZ118395.1)
print(save.em.2012.KZ118395.1, nlines = 12)

#2017
KZ118395.1_spl.desc <- summaryGeno(KZ118395.1_spl[25:35,-1], miss.val=c(0,NA))
save.em.2017.KZ118395.1 <- haplo.em(geno=KZ118395.1_spl[24:35,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2017.KZ118395.1)
print(save.em.2017.KZ118395.1, nlines = 12)



#KZ117131.1
KZ117131.1 <- read.vcfR("./KZ117131.1_only/KZ117131.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf")
KZ117131.1_loc <- vcfR2loci(KZ117131.1, return.alleles = TRUE)
str(KZ117131.1_loc)

#fixing colnames and subsetting by column
names(KZ117131.1_loc) <- gsub(x = names(KZ117131.1_loc), pattern = "KZ117131.1_", replacement = "")

importLoc <- as.numeric(as.character(names(KZ117131.1_loc)))
sub_importLoc <- importLoc[importLoc > 27000 & importLoc < 37000 & importLoc != 28766 & importLoc != 30592 & importLoc != 34050 & importLoc != 35474] 
#just picking one 10Kb region because whole scaffold looks to have sweep, plus had to remove 4 SNPs that created ambiguous haplotypes b/c of missing data.

idx <- as.character(sub_importLoc)

New_KZ117131.1_loc <- data.frame(KZ117131.1_loc[,idx])

#splitting alleles into separate columns
KZ117131.1_spl <- data.frame(rownames(KZ117131.1_loc))
names(KZ117131.1_spl) <- "sample"

for (i in 1:ncol(New_KZ117131.1_loc)){
  vec <- as.vector(New_KZ117131.1_loc[,i])
  temp <- revalue(vec, c("."="NA/NA"))
  temp <- data.frame(temp)
  temp$temp <- as.character(temp$temp)
  print(head(temp))
  name1 <- paste0(i,".","1")
  name2 <- paste0(i,".","2")
  alleles <- temp %>% separate(col = "temp", into = c(name1,name2), sep = "/")
  print(head(alleles))
  KZ117131.1_spl <- cbind(KZ117131.1_spl,alleles)
}


#2002
KZ117131.1_spl.desc <- summaryGeno(KZ117131.1_spl[1:13,-1], miss.val=c(0,NA))
save.em.2002.KZ117131.1 <- haplo.em(geno=KZ117131.1_spl[1:13,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2002.KZ117131.1)
print(save.em.2002.KZ117131.1, nlines = 12)

#2012
KZ117131.1_spl.desc <- summaryGeno(KZ117131.1_spl[14:24,-1], miss.val=c(0,NA))
save.em.2012.KZ117131.1 <- haplo.em(geno=KZ117131.1_spl[14:24,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2012.KZ117131.1)
print(save.em.2012.KZ117131.1, nlines = 12)

#2017
KZ117131.1_spl.desc <- summaryGeno(KZ117131.1_spl[25:35,-1], miss.val=c(0,NA))
save.em.2017.KZ117131.1 <- haplo.em(geno=KZ117131.1_spl[24:35,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2017.KZ117131.1)
print(save.em.2017.KZ117131.1, nlines = 12)
