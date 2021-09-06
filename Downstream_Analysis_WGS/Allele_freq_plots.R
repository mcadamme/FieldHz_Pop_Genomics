#This is the script I used to plot allele frequency change at Cry1 and Cry2 associated loci over time
#04122021 MF

library(tidyr); library(ggplot2)

setwd("~/Desktop/Hz_fieldColl_pop_gen/Reanalysis_PNAS")

#loading resistance-associated loci with signatures of temporal change in the field
Cry1 <- read.table("./hiCry1FST2002to2017.txt", header = T)
Cry2 <- read.table("./hiCry2FST2002to2017.txt", header = T)

#loading allele frequency data
freqs_2002 <- read.table("./thinned_FieldHzea_variantsonly_2002.freq.frq", header = F, sep="\t", stringsAsFactors = F)
freqs_2012 <- read.table("./thinned_FieldHzea_variantsonly_2012.freq.frq", header = F, sep="\t", stringsAsFactors = F)
freqs_2017 <- read.table("./thinned_FieldHzea_variantsonly_2017.freq.frq", header = F, sep="\t", stringsAsFactors = F)

#filtering based on n Chr - must have at least 6 individs per year
freqs_2002_filt <- subset(freqs_2002, V4 > 12)
freqs_2012_filt <- subset(freqs_2012, V4 > 12)
freqs_2017_filt <- subset(freqs_2017, V4 > 12)

#prepping for merge with Cry1 and Cry2 datasets
#Cry1
uniq_Cry1_scafs <- data.frame(unique(Cry1$Scaf))
names(uniq_Cry1_scafs) <- "scaf"

min_HiCry1ByScaf <- as.numeric(as.character(as.vector(tapply(Cry1$WinStart, Cry1$Scaf, min))))
max_HiCry1ByScaf <- as.numeric(as.character(as.vector(tapply(Cry1$WinStart, Cry1$Scaf, max))))

minMax_cry1 <- cbind(min_HiCry1ByScaf, max_HiCry1ByScaf)
str(minMax_cry1)
names(minMax_cry1) <- c("Min", "Max")

#Cry2
uniq_Cry2_scafs <- data.frame(unique(Cry2$Scaf))
names(uniq_Cry2_scafs) <- "scaf"

min_HiCry2ByScaf <- as.numeric(as.character(as.vector(tapply(Cry2$WinStart, Cry2$Scaf, min))))
max_HiCry2ByScaf <- as.numeric(as.character(as.vector(tapply(Cry2$WinStart, Cry2$Scaf, max))))

minMax_cry2 <- cbind(min_HiCry2ByScaf, max_HiCry2ByScaf)
str(minMax_cry2)
names(minMax_cry2) <- c("Min", "Max")


#Split allele from freqs
freqs_2002_spl <- freqs_2002_filt %>% separate(col = "V5", into = c("Allele1", "Freq1"), sep = ":") %>% separate(col = "V6", into = c("Allele2", "Freq2"), sep = ":")
freqs_2012_spl <- freqs_2012_filt %>% separate(col = "V5", into = c("Allele1", "Freq1"), sep = ":") %>% separate(col = "V6", into = c("Allele2", "Freq2"), sep = ":")
freqs_2017_spl <- freqs_2017_filt %>% separate(col = "V5", into = c("Allele1", "Freq1"), sep = ":") %>% separate(col = "V6", into = c("Allele2", "Freq2"), sep = ":")

#adding year for eventual rbind
freqs_2002_spl$Year <- rep("2002", times = nrow(freqs_2002_spl))
freqs_2012_spl$Year <- rep("2012", times = nrow(freqs_2012_spl))
freqs_2017_spl$Year <- rep("2017", times = nrow(freqs_2017_spl))

#adding ScafPos for rbind
freqs_2002_spl$ScafPos <- paste0(freqs_2002_spl$V1,"_",freqs_2002_spl$V2)
freqs_2012_spl$ScafPos <- paste0(freqs_2012_spl$V1,"_",freqs_2012_spl$V2)
freqs_2017_spl$ScafPos <- paste0(freqs_2017_spl$V1,"_",freqs_2017_spl$V2)


##### Cry1 First #####
#merging Cry1 with freq datasets by scaffold
freqs_2002_Cry1 <- merge(freqs_2002_spl, uniq_Cry1_scafs, by.x = "V1", by.y = "scaf")
freqs_2012_Cry1 <- merge(freqs_2012_spl, uniq_Cry1_scafs, by.x = "V1", by.y = "scaf")
freqs_2017_Cry1 <- merge(freqs_2017_spl, uniq_Cry1_scafs, by.x = "V1", by.y = "scaf")


merged <- rbind(freqs_2002_Cry1[,c(4:10)], freqs_2012_Cry1[,c(4:10)], freqs_2017_Cry1[,c(4:10)])#long format

#reshaping long to wide for subsetting
merged_W <- reshape(merged, idvar = c("ScafPos"), timevar = "Year", direction = "wide")
merged_W_spl <- merged_W %>% separate(col = "ScafPos", into = c("Scaf", "Pos"), sep = "_")
str(merged_W_spl)
merged_W_spl$Pos <- as.numeric(as.character(merged_W_spl$Pos))


#subsetting by Cry1 scaffolds & windows - commented out scaffolds that were not significant after we changed how we defined the QTL windows.
#Cry1_KZ118241.1 <- subset(merged_W_spl, Scaf == "KZ118241.1" & Pos > 140000 & Pos < (147000+1000))#adding 1000 to get full window
Cry1_KZ118067.1 <- subset(merged_W_spl, Scaf == "KZ118067.1" & Pos > 144000 & Pos < (146000+1000))
Cry1_KZ116099.1 <- subset(merged_W_spl, Scaf == "KZ116099.1" & Pos > 2000 & Pos < (7000+1000))
Cry1_KZ117975.1 <- subset(merged_W_spl, Scaf == "KZ117975.1" & Pos > 43000 & Pos < (46000+1000))
#Cry1_KZ118133.1 <- subset(merged_W_spl, Scaf == "KZ118133.1" & Pos > 34000 & Pos < (35000+1000))

#writing position tables
#write.csv(Cry1_KZ118241.1[,c(1,2)], file = "Cry1_KZ118241.1_pos.csv")
#write.csv(Cry1_KZ118067.1[,c(1,2)], file = "Cry1_KZ118067.1_pos.csv")
#write.csv(Cry1_KZ116099.1[,c(1,2)], file = "Cry1_KZ116099.1_pos.csv")
#write.csv(Cry1_KZ117975.1[,c(1,2)], file = "Cry1_KZ117975.1_pos.csv")
#write.csv(Cry1_KZ118133.1[,c(1,2)], file = "Cry1_KZ118133.1_pos.csv")

#All Cry1
#All_Cry1 <- rbind(Cry1_KZ118241.1, Cry1_KZ118067.1, Cry1_KZ116099.1, Cry1_KZ117975.1, Cry1_KZ118133.1)
All_Cry1 <- rbind(Cry1_KZ118067.1, Cry1_KZ116099.1, Cry1_KZ117975.1)
All_Cry1 <- na.omit(All_Cry1)
str(All_Cry1)

All_Cry1$Freq1.2002 <- as.numeric(as.character(All_Cry1$Freq1.2002))
All_Cry1$Freq2.2002 <- as.numeric(as.character(All_Cry1$Freq2.2002))
All_Cry1$Freq1.2012 <- as.numeric(as.character(All_Cry1$Freq1.2012))
All_Cry1$Freq2.2012 <- as.numeric(as.character(All_Cry1$Freq2.2012))
All_Cry1$Freq1.2017 <- as.numeric(as.character(All_Cry1$Freq1.2017))
All_Cry1$Freq2.2017 <- as.numeric(as.character(All_Cry1$Freq2.2017))


#checking for same ref allele (Allele1 across years)
all(All_Cry1$Allele1.2002 == All_Cry1$Allele1.2012)
all(All_Cry1$Allele1.2002 == All_Cry1$Allele1.2017)#good

#getting only changes > 0.1
All_Cry1$TotDiff <- abs(All_Cry1$Freq1.2002 - All_Cry1$Freq1.2017)

sub_All_Cry1 <- subset(All_Cry1, TotDiff > 0.1)

#getting highest 2002-2017 cry1 allele frequency change on each scaffold
Cry1_changes <- tapply(All_Cry1$TotDiff, All_Cry1$Scaf, max)

#plotting function
dat_for_Plot <- data.frame()

for (i in 1:nrow(sub_All_Cry1)){
  if (sub_All_Cry1[i,5] - sub_All_Cry1[i,10] < 0) {
    scaf <- rep(sub_All_Cry1[i,1], times = 3)
    row <- rep(i, times = 3)
    dat <- c(sub_All_Cry1[i,5], sub_All_Cry1[i,10], sub_All_Cry1[i,15])
    year <- c("2002", "2012", "2017")
    freq_change <- data.frame(cbind(scaf, row, dat, year))
    dat_for_Plot <- rbind(dat_for_Plot, freq_change)}
  else {
    scaf <- rep(sub_All_Cry1[i,1], times = 3)
    row <- rep(i, times = 3)
    dat <- c(sub_All_Cry1[i,7], sub_All_Cry1[i,12], sub_All_Cry1[i,17])
    year <- c("2002", "2012", "2017")
    freq_change <- data.frame(cbind(scaf,row, dat, year))
    dat_for_Plot <- rbind(dat_for_Plot, freq_change)}
}

str(dat_for_Plot)
dat_for_Plot$dat <- as.numeric(as.character(dat_for_Plot$dat))
ggplot(dat_for_Plot, aes(x = year, y = dat)) + ggtitle("Cry1_KZ118133.1" ) + geom_line(aes(color = scaf, group = row))

##### Cry2 #####
#merging Cry2 with freq datasets by scaffold
freqs_2002_Cry2 <- merge(freqs_2002_spl, uniq_Cry2_scafs, by.x = "V1", by.y = "scaf")
freqs_2012_Cry2 <- merge(freqs_2012_spl, uniq_Cry2_scafs, by.x = "V1", by.y = "scaf")
freqs_2017_Cry2 <- merge(freqs_2017_spl, uniq_Cry2_scafs, by.x = "V1", by.y = "scaf")


merged2 <- rbind(freqs_2002_Cry2[,c(4:10)], freqs_2012_Cry2[,c(4:10)], freqs_2017_Cry2[,c(4:10)])#long format

#reshaping long to wide for subsetting
merged_W2 <- reshape(merged2, idvar = c("ScafPos"), timevar = "Year", direction = "wide")
merged_W2_spl <- merged_W2 %>% separate(col = "ScafPos", into = c("Scaf", "Pos"), sep = "_")
str(merged_W2_spl)
merged_W2_spl$Pos <- as.numeric(as.character(merged_W2_spl$Pos))


#subsetting by Cry1 scaffolds & windows
Cry2_KZ117108.1 <- subset(merged_W2_spl, Pos > 257000 & Pos < (270000+1000))#adding 1000 to get full window

#All Cry2
All_Cry2 <- na.omit(Cry2_KZ117108.1)


All_Cry2$Freq1.2002 <- as.numeric(as.character(All_Cry2$Freq1.2002))
All_Cry2$Freq2.2002 <- as.numeric(as.character(All_Cry2$Freq2.2002))
All_Cry2$Freq1.2012 <- as.numeric(as.character(All_Cry2$Freq1.2012))
All_Cry2$Freq2.2012 <- as.numeric(as.character(All_Cry2$Freq2.2012))
All_Cry2$Freq1.2017 <- as.numeric(as.character(All_Cry2$Freq1.2017))
All_Cry2$Freq2.2017 <- as.numeric(as.character(All_Cry2$Freq2.2017))


#checking for same ref allele (Allele1 across years)
all(All_Cry2$Allele1.2002 == All_Cry2$Allele1.2012)
all(All_Cry2$Allele1.2002 == All_Cry2$Allele1.2017)#good

#getting only changes > 0.1
All_Cry2$TotDiff <- abs(All_Cry2$Freq1.2002 - All_Cry2$Freq1.2017)

sub_All_Cry2 <- subset(All_Cry2, TotDiff > 0.1)

#getting highest 2002-2017 cry2 allele frequency change on each scaffold
Cry2_changes <- tapply(All_Cry2$TotDiff, All_Cry2$Scaf, max)


#plotting function
dat_for_Plot <- data.frame()

for (i in 1:nrow(sub_All_Cry2)){
  if (sub_All_Cry2[i,5] - sub_All_Cry2[i,10] < 0) {
    scaf <- rep(sub_All_Cry2[i,1], times = 3)
    row <- rep(i, times = 3)
    dat <- c(sub_All_Cry2[i,5], sub_All_Cry2[i,10], sub_All_Cry2[i,15])
    year <- c("2002", "2012", "2017")
    freq_change <- data.frame(cbind(scaf, row, dat, year))
    dat_for_Plot <- rbind(dat_for_Plot, freq_change)}
  else {
    scaf <- rep(sub_All_Cry2[i,1], times = 3)
    row <- rep(i, times = 3)
    dat <- c(sub_All_Cry2[i,7], sub_All_Cry2[i,12], sub_All_Cry2[i,17])
    year <- c("2002", "2012", "2017")
    freq_change <- data.frame(cbind(scaf,row, dat, year))
    dat_for_Plot <- rbind(dat_for_Plot, freq_change)}
}

str(dat_for_Plot)
dat_for_Plot$dat <- as.numeric(as.character(dat_for_Plot$dat))
ggplot(dat_for_Plot, aes(x = year, y = dat)) + geom_line(aes(color = scaf, group = row))


##### LG5 changes #####
#freq datasets by scaffold - scaffold 569
freqs_2002_569 <- subset(freqs_2002_spl, V1 == "KZ118395.1")
freqs_2012_569 <- subset(freqs_2012_spl, V1 == "KZ118395.1")
freqs_2017_569 <- subset(freqs_2017_spl, V1 == "KZ118395.1")

merged_569 <- rbind(freqs_2002_569, freqs_2012_569, freqs_2017_569)

#reshaping long to wide for subsetting
merged_569 <- reshape(merged_569, idvar = c("ScafPos"), timevar = "Year", direction = "wide")
merged_569_spl <- merged_569 %>% separate(col = "ScafPos", into = c("Scaf", "Pos"), sep = "_")
str(merged_569_spl)
merged_569_spl$Pos <- as.numeric(as.character(merged_569_spl$Pos))


nonCry_KZ118395.1 <- subset(merged_569_spl, Pos < (151000+1000))#adding 1000 to get full window

nonCry_KZ118395.1 <- na.omit(nonCry_KZ118395.1)
nonCry_KZ118395.1$Freq1.2002 <- as.numeric(as.character(nonCry_KZ118395.1$Freq1.2002))
nonCry_KZ118395.1$Freq2.2002 <- as.numeric(as.character(nonCry_KZ118395.1$Freq2.2002))
nonCry_KZ118395.1$Freq1.2012 <- as.numeric(as.character(nonCry_KZ118395.1$Freq1.2012))
nonCry_KZ118395.1$Freq2.2012 <- as.numeric(as.character(nonCry_KZ118395.1$Freq2.2012))
nonCry_KZ118395.1$Freq1.2017 <- as.numeric(as.character(nonCry_KZ118395.1$Freq1.2017))
nonCry_KZ118395.1$Freq2.2017 <- as.numeric(as.character(nonCry_KZ118395.1$Freq2.2017))


#checking for same ref allele (Allele1 across years)
all(nonCry_KZ118395.1$Allele1.2002 == nonCry_KZ118395.1$Allele1.2012)
all(nonCry_KZ118395.1$Allele1.2002 == nonCry_KZ118395.1$Allele1.2017)#good

#getting only changes > 0.1
nonCry_KZ118395.1$TotDiff <- abs(nonCry_KZ118395.1$Freq1.2002 - nonCry_KZ118395.1$Freq1.2017)

#getting highest 2002-2017 cry2 allele frequency change on each scaffold
nonCry_KZ118395.1_changes <- tapply(nonCry_KZ118395.1$TotDiff, nonCry_KZ118395.1$Scaf, max)


#scaffold 1612
freqs_2002_1612 <- subset(freqs_2002_spl, V1 == "KZ117131.1")
freqs_2012_1612 <- subset(freqs_2012_spl, V1 == "KZ117131.1")
freqs_2017_1612 <- subset(freqs_2017_spl, V1 == "KZ117131.1")

merged_1612 <- rbind(freqs_2002_1612, freqs_2012_1612, freqs_2017_1612)


#reshaping long to wide for subsetting
merged_1612 <- reshape(merged_1612, idvar = c("ScafPos"), timevar = "Year", direction = "wide")
merged_1612_spl <- merged_1612 %>% separate(col = "ScafPos", into = c("Scaf", "Pos"), sep = "_")
str(merged_1612_spl)
merged_1612_spl$Pos <- as.numeric(as.character(merged_1612_spl$Pos))


nonCry_KZ117131.1 <- subset(merged_1612_spl, Pos < (151000+1000))#adding 1000 to get full window

nonCry_KZ117131.1 <- na.omit(nonCry_KZ117131.1)
nonCry_KZ117131.1$Freq1.2002 <- as.numeric(as.character(nonCry_KZ117131.1$Freq1.2002))
nonCry_KZ117131.1$Freq2.2002 <- as.numeric(as.character(nonCry_KZ117131.1$Freq2.2002))
nonCry_KZ117131.1$Freq1.2012 <- as.numeric(as.character(nonCry_KZ117131.1$Freq1.2012))
nonCry_KZ117131.1$Freq2.2012 <- as.numeric(as.character(nonCry_KZ117131.1$Freq2.2012))
nonCry_KZ117131.1$Freq1.2017 <- as.numeric(as.character(nonCry_KZ117131.1$Freq1.2017))
nonCry_KZ117131.1$Freq2.2017 <- as.numeric(as.character(nonCry_KZ117131.1$Freq2.2017))


#checking for same ref allele (Allele1 across years)
all(nonCry_KZ117131.1$Allele1.2002 == nonCry_KZ117131.1$Allele1.2012)
all(nonCry_KZ117131.1$Allele1.2002 == nonCry_KZ117131.1$Allele1.2017)#good

#getting only changes > 0.1
nonCry_KZ117131.1$TotDiff <- abs(nonCry_KZ117131.1$Freq1.2002 - nonCry_KZ117131.1$Freq1.2017)

#getting highest 2002-2017 cry2 allele frequency change on each scaffold
nonCry_KZ117131.1_changes <- tapply(nonCry_KZ117131.1$TotDiff, nonCry_KZ117131.1$Scaf, max)

