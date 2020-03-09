#Script to get distribution of FST values under neutral model using ms simulated genos
#03052020


setwd("~/ms_sims/FieldPops/Ne_twentythou_NoSize/")

#Loading packages
x <- c("adegenet", "hierfstat")

lapply(x, FUN = function(X) {
  do.call("library", list(X)) 
})

#function to paste alleles 1 and 2 together on same line.
f <- function(x) {
  s <- seq(2, length(x), 2)
  paste(x[s-1], x[s], sep="|")
}

out.file <- ""
dim.data <- ""

file.names1 <- dir("./", pattern ="output.txt")

for(i in 1:length(file.names1)){

  df <- read.csv(file.names1[i], header = F)
  
  # run algorithm for each column
  df2 <- as.data.frame(lapply(df, f), stringsAsFactors=FALSE)
  dim_data <- dim(df2)
  dim.data <- append(dim.data, dim_data[2])
  
  if(dim_data[2] > 9){
    
    pop <- rep(c("P", "A"), times = 1, each = 12)
    
    df2_genind <- df2genind(df2, pop=pop, sep = "\\|", ploidy = 2)
    df2_genind_sum <- summary(df2_genind)
    df2_hier <- genind2hierfstat(df2_genind, pop = pop)
    
    FST_mat <- pairwise.WCfst(df2_hier, diploid = T)
    print(FST_mat[1,2])
    out.file <- rbind(out.file, round(FST_mat[1,2], digits = 6))
  }
}


write.table(out.file[1:10001], "FST_msSims4.out", row.names = F, col.names = F)