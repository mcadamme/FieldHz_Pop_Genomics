#This is the script I used to plot read counts for all Hz populations
#MF 10/29/2018

#To get the line counts, I used the following command in BASH:
#for file in *.fq.gz; do zcat $file | wc -l >> file_line_count; done

setwd("/media/megan/New Volume/Hz_PopGen_ddRAD_demult/AMD_libs1thru8_extended_frags/renamed_samples/")

line_counts <- read.table("file_line_count", header = F)
file_names <- read.table("file_list", header = F)

data <- cbind(file_names, line_counts)

colnames(data) <- c("file_names", "line_counts")

data$read_counts <- (data$line_counts/4)

data$year <- regmatches(data$file_names, regexpr("20[[:digit:]]+", data$file_names))
data$samp_num <- (1:265)

mean_read_count <- mean(data$read_counts)
stdev_read_count <- sd(data$read_counts)
y.low <- rep(mean_read_count*0.1, times = 265)
y.high <- rep(max(data$read_counts), times = 265)

exclude_from_analysis <- subset(data, read_counts < (mean_read_count * 0.1))
exclude_from_analysis #all have fewer than 40,000 reads 
#these individuals should be excluded from pop-level analysis with vcftools.

col_year = c("red","green","blue","black" )[as.factor(data$year)] 
         
png("Hz_fieldcoll_readcounts.png", units = "px", width = 600, height = 400)
plot(data$read_counts~data$samp_num, col = col_year, pch = 16, xaxt = "n", 
     ylab = "Read Counts", xlab = "Samples Grouped By Year", cex.lab = 1.25, cex.axis = 1.2)
lines(data$samp_num, y.low, col = 'grey')
lines(data$samp_num, y.high, col = 'grey')

polygon(c(data$samp_num, rev(data$samp_num)), c(y.high, rev(y.low)),
        col=rgb(0.6,0.6,0.7,0.2), border = NA)
abline(a = mean_read_count, b = 0)

dev.off()



