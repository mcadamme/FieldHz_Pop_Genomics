#This is the script I used to plot read counts for all Hz populations
#MF 10/29/2018

#To get the line counts, I used the following commands in BASH:
#for file in *.fq.gz; do zcat $file | wc -l >> file_line_count; done
#for file in *.fq.gz; do zcat $file | echo $file >> file_list; done

setwd("/media/megan/New Volume/Hz_PopGen_ddRAD_demult/AMD_libs1thru8_extended_frags/renamed_samples/")

line_counts <- read.table("file_line_count", header = F)
file_names <- read.table("file_list", header = F)

data <- cbind(file_names, line_counts)

colnames(data) <- c("file_names", "line_counts")

data$read_counts <- (data$line_counts/4)

data$year <- regmatches(data$file_names, regexpr("20[[:digit:]]+", data$file_names))
data$samp_num <- (1:nrow(data))

mean_read_count <- mean(data$read_counts)
sd_read_count <- sd(data$read_counts)
sum_read_count <- sum(data$read_counts)
y.low <- rep(mean_read_count*0.01, times = nrow(data))
y.high <- rep(max(data$read_counts), times = nrow(data))

exclude_from_analysis <- subset(data, read_counts < (mean_read_count*0.01))
exclude_from_analysis #all have fewer than 4,000 reads 
#these individuals should be excluded from pop-level analysis with vcftools.

col_year = c("red","green","blue","black" )[as.factor(data$year)] 
         
png("Hz_fieldcoll_readcounts.png", units = "px", width = 600, height = 400)
plot(log10(data$read_counts)~data$samp_num, col = col_year, pch = 16, xaxt = "n", 
     ylab = "Log10 Read Counts", xlab = "Samples Grouped By Year", ylim = c(2,8), cex.lab = 1.25, cex.axis = 1.2)
lines(data$samp_num, log10(y.low), col = 'grey')
lines(data$samp_num, log10(y.high), col = 'grey')

polygon(c(data$samp_num, rev(data$samp_num)), c((log10(y.high)), rev(log10(y.low))),
        col=rgb(0.6,0.6,0.7,0.2), border = NA)
abline(a = log10(mean_read_count), b = 0)

dev.off()



