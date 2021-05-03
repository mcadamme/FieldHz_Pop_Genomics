#getting changes in haplotype frequencies

library(vcfR); library(haplo.stats); library(tidyr); library(plyr)

setwd("/media/fritzlab/EE9C16C89C168AEB/Reanalysis_PNAS")

##### Chr7 #####

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



##### Chr1 #####

#KZ118241.1 - loading datasets
KZ118241.1_pos <- read.csv("Cry1_KZ118241.1_pos.csv", header = T)
KZ118241.1 <- read.vcfR("./KZ118241.1_only/KZ118241.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf")
KZ118241.1_loc <- vcfR2loci(KZ118241.1, return.alleles = TRUE)
str(KZ118241.1_loc)

#fixing colnames and subsetting by column
names(KZ118241.1_loc) <- gsub(x = names(KZ118241.1_loc), pattern = "KZ118241.1_", replacement = "")

importLoc <- as.character(names(KZ118241.1_loc))

idx <- match(as.character(KZ118241.1_pos$Pos), importLoc)
idx <- sort(idx)

New_KZ118241.1_loc <- KZ118241.1_loc[,idx] 

#splitting alleles into separate columns
KZ118241.1_spl <- data.frame(rownames(New_KZ118241.1_loc))
names(KZ118241.1_spl) <- "sample"

for (i in 1:ncol(New_KZ118241.1_loc)){
  vec <- as.vector(New_KZ118241.1_loc[,i])
  temp <- revalue(vec, c("."="NA/NA"))
  temp <- data.frame(temp)
  temp$temp <- as.character(temp$temp)
  print(head(temp))
  name1 <- paste0(i,".","1")
  name2 <- paste0(i,".","2")
  alleles <- temp %>% separate(col = "temp", into = c(name1,name2), sep = "/")
  print(head(alleles))
  KZ118241.1_spl <- cbind(KZ118241.1_spl,alleles)
}


#2002
KZ118241.1_spl.desc <- summaryGeno(KZ118241.1_spl[1:13,-1], miss.val=c(0,NA))
save.em.2002.KZ118241.1 <- haplo.em(geno=KZ118241.1_spl[1:13,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2002.KZ118241.1)
print(save.em.2002.KZ118241.1, nlines = 12)

#2012
KZ118241.1_spl.desc <- summaryGeno(KZ118241.1_spl[14:24,-1], miss.val=c(0,NA))
save.em.2012.KZ118241.1 <- haplo.em(geno=KZ118241.1_spl[14:24,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2012.KZ118241.1)
print(save.em.2012.KZ118241.1, nlines = 12)

#2017
KZ118241.1_spl.desc <- summaryGeno(KZ118241.1_spl[25:35,-1], miss.val=c(0,NA))
save.em.2017.KZ118241.1 <- haplo.em(geno=KZ118241.1_spl[24:35,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2017.KZ118241.1)
print(save.em.2017.KZ118241.1, nlines = 12)


##### Chr11 #####

#KZ118133.1 - loading datasets
KZ118133.1_pos <- read.csv("Cry1_KZ118133.1_pos.csv", header = T)
max(KZ118133.1_pos$Pos)
min(KZ118133.1_pos$Pos)

KZ118133.1 <- read.vcfR("./KZ118133.1_only/KZ118133.1_thinned_FieldHzea_variantsonly.vcf.recode.vcf")
KZ118133.1_loc <- vcfR2loci(KZ118133.1, return.alleles = TRUE)
str(KZ118133.1_loc)

#fixing colnames and subsetting by column
names(KZ118133.1_loc) <- gsub(x = names(KZ118133.1_loc), pattern = "KZ118133.1_", replacement = "")

importLoc <- as.character(names(KZ118133.1_loc))

idx <- match(as.character(KZ118133.1_pos$Pos), importLoc)
idx <- sort(idx)

New_KZ118133.1_loc <- data.frame(KZ118133.1_loc[,idx])
#names(New_KZ118133.1_loc) <- "19"

table(New_KZ118133.1_loc[1:13,])#2002 freq
table(New_KZ118133.1_loc[14:24,])#2012 freq
table(New_KZ118133.1_loc[25:35,])#2017 freq

#splitting alleles into separate columns
KZ118133.1_spl <- data.frame(rownames(KZ118133.1_loc))
names(KZ118133.1_spl) <- "sample"

for (i in 1:ncol(New_KZ118133.1_loc)){
  vec <- as.vector(New_KZ118133.1_loc[,i])
  temp <- revalue(vec, c("."="NA/NA"))
  temp <- data.frame(temp)
  temp$temp <- as.character(temp$temp)
  print(head(temp))
  name1 <- paste0(i,".","1")
  name2 <- paste0(i,".","2")
  alleles <- temp %>% separate(col = "temp", into = c(name1,name2), sep = "/")
  print(head(alleles))
  KZ118133.1_spl <- cbind(KZ118133.1_spl,alleles)
}


#2002
KZ118133.1_spl.desc <- summaryGeno(KZ118133.1_spl[1:13,-1], miss.val=c(0,NA))
save.em.2002.KZ118133.1 <- haplo.em(geno=KZ118133.1_spl[1:13,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2002.KZ118133.1)
print(save.em.2002.KZ118133.1, nlines = 12)

#2012
KZ118133.1_spl.desc <- summaryGeno(KZ118133.1_spl[14:24,-1], miss.val=c(0,NA))
save.em.2012.KZ118133.1 <- haplo.em(geno=KZ118133.1_spl[14:24,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2012.KZ118133.1)
print(save.em.2012.KZ118133.1, nlines = 12)

#2017
KZ118133.1_spl.desc <- summaryGeno(KZ118133.1_spl[25:35,-1], miss.val=c(0,NA))
save.em.2017.KZ118133.1 <- haplo.em(geno=KZ118133.1_spl[24:35,-1], locus.label=idx, miss.val=c(0,NA), 
                                    control = haplo.em.control(n.try = 20, insert.batch.size=2))
names(save.em.2017.KZ118133.1)
print(save.em.2017.KZ118133.1, nlines = 12)