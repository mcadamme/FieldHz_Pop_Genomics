#This is the script that I used to examine changes in allele frequencies for significantly diverged SNPs.

#used equations on pg 28 of Falconer and Mackay and solved for s.
library(scatterplot3d)

#function for selection coefficient assuming dominance of p
Dom.Sel.Coef = function (q.1, q.2, q.3, g.1, g.2) {
  delta.q1 = (q.2-q.1)/g.1
  delta.q2 = (q.3-q.2)/g.2
  s1 = (delta.q1)/((q.1^2)*(q.1+delta.q1-1))
  s2 = (delta.q2)/((q.2^2)*(q.2+delta.q2-1))
  x <- data.frame(cbind(s1,s2))
print(x)
}

#function for selection coefficient assuming incomplete dominance (or no dominance) of p
NoDom.Sel.Coef = function (q.1, q.2, q.3, g.1, g.2) {
  delta.q1 = (q.2-q.1)/g.1
  delta.q2 = (q.3-q.2)/g.2
  s1 = (delta.q1)/(q.1 *(0.5*q.1+delta.q1-0.5))
  s2 = (delta.q2)/(q.2 *(0.5*q.2+delta.q2-0.5))
  y <- data.frame(cbind(s1,s2))
print(y)
}

#function for selection coefficient assuming recessiveness of p
Rec.Sel.Coef = function (q.1, q.2, q.3, g.1, g.2) {
  delta.q1 = (q.2-q.1)/g.1
  delta.q2 = (q.3-q.2)/g.2
  s1 = (delta.q1)/(((-(1-q.1)^2)*q.1)+(2*q.1*(1-q.1)*delta.q1)+(q.1^2*delta.q1))
  s2 = (delta.q2)/(((-(1-q.2)^2)*q.2)+(2*q.2*(1-q.2)*delta.q2)+(q.2^2*delta.q2))
  z <- data.frame(cbind(s1,s2))
print(z)
}

data <- read.csv("~/Desktop/Hz_ProtCoding_deltaq.csv", header = T)
freqs_only <- data[,c(4:6)]

first_gens <- c(60,60,60,60,60,60)
second_gens <- c(30,30,30,30,30,30)

freqs_gen <- as.matrix(cbind(freqs_only, first_gens, second_gens))
str(freqs_gen)

Dom_df <- data.frame(s1=numeric(),s2=numeric())
                 
for (i in 1:nrow(freqs_gen)){
  x <- Dom.Sel.Coef(freqs_gen[i,1], freqs_gen[i,2], freqs_gen[i,3], freqs_gen[i,4], freqs_gen[i,5])
  Dom_df <- rbind(Dom_df,x)
}
  
NoDom_df <- data.frame(s1=numeric(),s2=numeric())

for (i in 1:nrow(freqs_gen)){
  x <- NoDom.Sel.Coef(freqs_gen[i,1], freqs_gen[i,2], freqs_gen[i,3], freqs_gen[i,4], freqs_gen[i,5])
  NoDom_df <- rbind(NoDom_df,x)
}

Rec_df <- data.frame(s1=numeric(),s2=numeric())

for (i in 1:nrow(freqs_gen)){
  x <- Rec.Sel.Coef(freqs_gen[i,1], freqs_gen[i,2], freqs_gen[i,3], freqs_gen[i,4], freqs_gen[i,5])
  Rec_df <- rbind(Rec_df,x)
}


#To plot selection coefficients
plot.Rec.Sel.Coef = function(q.1, q.2, g, s, dataframe) {
  delta.q = (q.2-q.1)/g
  df.name <- deparse(substitute(dataframe))
  png(file = print(paste0(df.name,"Rec.png")), units = "px", height = 600, width = 900)
  scatterplot3d(abs(delta.q), q.1, abs(s),highlight.3d = TRUE, col.axis = "blue", 
                cex = 2.5, cex.axis = 1.5, angle = 15, cex.lab = 2, cex.main = 2, col.grid = "lightblue", 
                main = "Recessiveness of p", xlab = "Delta q", ylab = "",
                zlab = "Selection Coefficient", pch = 20, zlim = c(0,.8),xlim = c(0,0.025))
  dev.off()
}
