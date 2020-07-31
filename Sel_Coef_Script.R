#This is the script that I used to examine changes in allele frequencies for significantly diverged SNPs.

library(tidyr); library(ggplot2)

#used equations on pg 28 of Falconer and Mackay and solved for s.

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

data <- read.csv("~/Desktop/Hz_ProtCoding_deltaq.csv", header = T) #top 5 sel ancestral, bottom sel alt.
freqs_only <- data[,c(4:6)]

first_gens <- c(50,50,50,50,50,50)
second_gens <- c(25,25,25,25,25,25)

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

All_data <- rbind(Dom_df, NoDom_df, Rec_df)
dom <- rep(c("Dom", "NoDom", "Rec"), each = 6, times = 1)
All_data2 <- cbind(All_data, dom)
colnames(All_data2) <- c("2002", "2012", "dom")

All_data_long <- gather(All_data2, sel_per, sel_coef, "2002":"2012", factor_key=TRUE)

#without standard errors
interaction.plot(All_data_long$sel_per, All_data_long$dom, All_data_long$sel_coef,
                 fun = mean,  # summary statistic to be plotted for response variable
                 type = "l",     # type of plot, here "l" for lines
                 ylab = "Selection coefficient",
                 xlab = "Selection period (Year)",
                 col = c("blue4", "red4", "black"),
                 lty = 1,  # line type
                 lwd = 2,  # line width
                 trace.label = "Degree of Dom",  # label for legend
                 xpd = FALSE) #,  # 'clip' legend at border)

#with standard errors
df <- with(All_data_long, aggregate(sel_coef, list(dom=dom, sel_per=sel_per), mean))
df$se <- with(All_data_long, aggregate(sel_coef, list(dom=dom, sel_per=sel_per), function(x) sd(x)/sqrt(6)))[,3]

ggplot(df, aes(x=sel_per, y=x, colour=dom, group=dom)) + geom_line(aes(linetype=dom), size=.6) + 
  geom_point(aes(shape=dom), size=3) + geom_errorbar(aes(ymax=x+se, ymin=x-se), width=.1) +
  labs(x=expression("Selection Period"), y=expression("Selection coefficient (s)")) +
  scale_color_manual(values=c("#868686FF","#0073C2FF","#ff1f39")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"), panel.border = element_rect(colour = "black", fill=NA, size=2))

