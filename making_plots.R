library(latex2exp)

load(".Rdata")

setwd("./plots")


#-------------------------------------------------------------------------------

pdf("reported_cases_to_N1_PMMoV.pdf",width = 6.5,height = 8)
par(mar = c(3,3.1,1.9,3.1), mfrow = c(5,3), oma = c(0.1,0,0,0))
for (i in 1:J) {
  plot(X2[[i]], pch = 21, bg = "blue", cex = 0.9,
       main = paste("Facility",LETTERS[i]),xlab="",ylab="")
  mtext("time (in weeks)",side=1,line=2,cex=0.68)
  mtext(TeX("N1/PMMoV$\\,\\times 10^3$"),side=2,line=2,cex=0.55)
  par(new=T)
  plot(Z_ts[[i]],xlab="",ylab="",cex=0.85,pch=21,bg="grey",axes=F)
  axis(4, at = pretty(range(Z_ts[[i]])))
  mtext(TeX("Reported Cases"),side=4,cex=0.62,line=2)
}

plot(X2[[1]],xlab = "",ylab="",axes=F,col="white")
legend("topleft", pch = 21, pt.bg = c("blue","grey"), bty="n",cex=1.15,
       legend = TeX(c("N1/PMMoV$\\,\\times 10^3$","Reported Cases")))
dev.off()

#-------------------------------------------------------------------------------
# Beta Plots

pdf("beta1_plots.pdf",height=10,width = 4)
par(mar = c(2.4,3.2,1.8,1), mfrow = c(7,2), oma = c(0.5,0,0,0))
for (i in 1:7) {
  plot(beta_list[[i]][1,], ylab = "", xlab = "",cex=0.75)
  mtext("Iterations",side=1,line=1.9,cex=0.62)
  mtext(TeX("$\\beta_{1s}$"),side=2,line=2,cex = 0.75)
  
  plot(beta_list[[i]][2,], ylab = "", xlab = "",cex=0.75)     
  mtext("Iterations",side=1,line=1.9,cex=0.62)
  mtext(TeX("$\\beta_{2s}$"),side=2,line=2,cex = 0.75)
  par(new=T)
  mtext(paste("Facility",LETTERS[i]),side=3,adj=-0.6,line=0.6,cex=0.8)
}
dev.off()

pdf("beta2_plots.pdf", height=10,width = 4)
par(mar = c(2.4,3.2,1.8,1), mfrow = c(7,2), oma = c(0.5,0,0,0))
for (i in 8:J) {
  plot(beta_list[[i]][1,], ylab = "", xlab = "",cex=0.75)
  mtext("Iterations",side=1,line=1.9,cex=0.62)
  mtext(TeX("$\\beta_{1s}$"),side=2,line=2,cex = 0.75)
  
  plot(beta_list[[i]][2,], ylab = "", xlab = "",cex=0.75)     
  mtext("Iterations",side=1,line=1.9,cex=0.62)
  mtext(TeX("$\\beta_{2s}$"),side=2,line=2,cex = 0.75)
  par(new=T)
  mtext(paste("Facility",LETTERS[i]),side=3,adj=-0.6,line=0.6,cex=0.8)
}
dev.off()


#-------------------------------------------------------------------------------
# a_s Plots

pdf("a_s_plots.pdf", width = 5, height = 7)
par(mar = c(2.9,3,1.8,0.5), mfrow = c(5,3), oma = c(0,0,0,0))
for (i in 1:J) {
  plot(a_s[[i]],main = "",ylab = "", xlab = "",cex=0.78)
  mtext(paste("Facility",LETTERS[i]), side = 3,line = 0.45,cex=0.74)
  mtext("Iterations",side=1,line=1.9,cex=0.58)
  mtext(TeX("$a_{s}$"),side=2,line=1.9,cex = 0.75)
}

dev.off()


#-------------------------------------------------------------------------------
# Reporting probs plot

pdf("report_probs.pdf",width=6.4,height=8)
par(mar = c(3.05,3.05,2,3.1), mfrow = c(5,3), oma = c(0.1,0,0,0))

for (i in 1:J) {
  plot(reportProbs_est[[i]][,1],xlab="",ylab="",cex=0.8,bg="seagreen",
       pch = 21, main = paste("Facility",LETTERS[i]), lwd=0.5,
       ylim = c(0,1))
  mtext(TeX("$\\pi_{ts}$"), side = 2, line = 2, cex = 0.6)
  mtext("time (in weeks)", side  = 1, line = 2, cex = 0.5)
  par(new = T)
  plot(Z_ts[[i]], xlab="",ylab="",axes=F,bg="grey",cex=0.77,pch=21,lwd=0.5)
  axis(4, at = pretty(range(Z_ts[[i]]))) 
  mtext("Reported Cases",side = 4, line=2,cex = 0.56)
}

plot(reportProbs_est[[14]][,1], xlab="", ylab="", col="white", axes = F)
legend("topleft", pch = 21, pt.bg = c("seagreen","grey"), bty = "n",
       legend = TeX(c("Mean $\\pi_{ts}$","Reported Cases"))
       ,cex = 1.2)
dev.off()

#-------------------------------------------------------------------------------
# Negative Binomial Probs

pdf("negBin_probs.pdf",width=6.2,height=8)
par(mar = c(3,3.2,2,3.1), mfrow = c(5,3), oma = c(0.1,0,0,0))

for (i in 1:J) {
  plot(negBinProbs_est[[i]][,1], xlab = "", ylab = "", pch = 21,
       main = paste("Facility",LETTERS[i]), cex=0.82, bg=2,lwd=0.5,
       ylim = c(0,1))
  mtext(TeX("$p_{ts}$"),side=2,line=2,cex=0.7)
  mtext("time (in weeks)",side=1,line=2,cex=0.6)
  par(new = T)
  plot(X2[[i]],xlab = "",ylab = "",axes = F,col="blue",
       cex=0.8,lwd=1,pch=21)
  axis(4, at = pretty(range(X2[[i]]))) 
  mtext(TeX("N1/PMMoV$\\times 10^3$"),side = 4, line=2.1,cex = 0.55)
}

plot(negBinProbs_est[[14]][,1], xlab="", ylab="", col="white", axes = F)
legend(-3,1, pch = c(21,1), col = c(1,"blue"),
       pt.bg = 2 ,bty = "n",cex=1.3, 
       legend = TeX(c("Mean $p_{ts}$","N1/PMMoV$\\times 10^3$")))

dev.off()


#-------------------------------------------------------------------------------
# Estimated Cases

pdf("est_N_ts_versus_Z_ts.pdf",width=5.6,height=8)
par(mar = c(2.9,3.1,1.9,0.7), mfrow = c(5,3), oma = c(0.1,0,0,0))

for (i in 1:J) {
  plot(N_t_est[[i]][,1], xlab = "", ylab = "", pch = 21, bg = 2,cex=0.95,
       main = paste("Facility",LETTERS[i]), ylim = c(0,max(N_t_est[[i]])))
  mtext(TeX("Number of Cases"),side=2,line=2,cex=0.65)
  mtext("time (in weeks)",side=1,line=2,cex=0.6)
  points(Z_ts[[i]], bg ="grey", pch = 21,cex=0.9)
}

plot(N_t_est[[1]][,1], xlab = "", ylab = "", col = "white", axes = F)
legend("topleft", pch = 21, pt.bg = c(2,"grey"), bty = "n", cex = 1.3, 
       legend = TeX(c("Mean $N_{ts}$","Reported $Z_{ts}$")))
dev.off()

#-------------------------------------------------------------------------------
# Finding confidence interval of N_ts 

quantiles <- c(0.05,0.95)
N_ts_conf <- list()

for (s in unique(location)) {
  N_ts_conf[[s]] <- apply(N_ts[[s]],1,quantile,quantiles)
}

pdf("estimated_N_ts_with_CI.pdf", width = 7.2, height = 8)
par(mar = c(2.9,3.1,2,0.8), mfrow = c(5,3), oma = c(0,0,0,0))
for(s in 1:J) {
  plot(N_t_est[[s]][,"mean"],xlab="",ylab="",pch = 21, bg = 2, 
       cex=0.95, main = paste("Facility",LETTERS[s]), lwd = 0.3,
       ylim = c(0,max(N_ts_conf[[s]])))
  mtext("time (in weeks)", side  = 1, line = 2, cex = 0.68)
  mtext(TeX("Estimated $N_{ts}$"),side=2,line=2,cex=0.65)
  arrows(x0=1:results$time_points[[s]],y0=N_ts_conf[[s]][1,], 
         x1=1:results$time_points[[s]],y1=N_ts_conf[[s]][2,], 
         code=3,angle=90, length=0.05,lty=1,lwd=0.8)
  points(N_t_est[[s]][,"mean"],pch = 21, bg = 2, cex = 0.95)
}
plot(N_t_est[[1]][,1], xlab = "", ylab = "", col = "white", axes = F)
legend("topleft", pch = c(21,21), pt.bg = c(2,"grey"),bty = "n", cex = 1.3, 
       lty = c(1,0),
       legend = TeX(c("Mean $N_{ts}$ with 90% CI")))
dev.off()

#-------------------------------------------------------------------------------
# Method 1: Plotting N1/PMMoV Values

pdf("N1_PMMoV_vals_m1.pdf",width = 7,height = 8.8)
par(mar = c(2.9,3.1,1.9,3.1), mfrow = c(5,3), oma = c(0.1,0,0,0))
for (s in 1:J) {
  plot(X2[[s]], xlab ="",ylab="",main = paste("Facility",LETTERS[s]),
      lwd=0.9,col="steelblue",pch=21)
  points(above_thresh[[s]][[4]],N1_PMMoV_ts[[s]][[4]]*1e3,bg="blue",pch=21,lwd=0.5)
  mtext(TeX("N1/PMMoV$\\times 10^3$"),side=2,line=1.9,cex=0.58)
  mtext("time (in weeks)",side=1,line=1.9,cex=0.62)
  par(new=T)
  plot(above_thresh[[s]][[4]],negBinProbs_est[[s]][above_thresh[[s]][[4]],1],
       pch = 21, col = 2, xlab = "", ylab = "", axes = F,
       xlim = c(1,results$time_points[[s]]))
  axis(4, at = pretty(range(negBinProbs_est[[s]][above_thresh[[s]][[4]],1])))
  mtext(TeX("Value of $p_{ts}$"), side=4, line=2,cex=0.6)
}
plot(X2[[1]],xlab = "",ylab="",axes=F,col="white")
legend("topleft", pch = c(21,21,21), col = c("steelblue",1,2),
       pt.bg = c("white","blue","white"),cex = 1.2,  bty = "n", 
       legend = TeX(c("N1/PMMoV", "N1/PMMoV confirming at least one case",
                      "$p_{ts} \\geq 0.8$")))
dev.off()

#-------------------------------------------------------------------------------
# Mins of N1/PMMoV at different thresholds

pdf("N1_PMMoV_mins.pdf",width = 6.2,height = 7.7)
options(digits = 2, scipen = -2)
par(mar = c(3,3.2,2,0.9), mfrow = c(5,3), oma = c(0.1,0.1,0,0))
for (i in 1:J) {
  plot(c(5:9/10,0.95),N1_PMMoV_mins[[i]], pch = 21, bg = "blue", cex=0.9,
       main = paste("Facility",LETTERS[i]),xlim=c(0.5,0.95),ylab="")
  lines(c(5:9/10,0.95),N1_PMMoV_mins[[i]],col="blue")
  mtext(TeX("$p_{ts}$"),side=1,line=2.2,cex=0.8)
  mtext("N1/PMMoV",side=2,line=2.2,cex=0.7)
}

dev.off()


#-------------------------------------------------------------------------------
# Method 2: Plotting Signal Value

pdf("N1_PMMoV_vals_m2.pdf",width = 6.5,height = 8)
par(mar = c(2.9,3.1,2,3.1), mfrow = c(5,3), oma = c(0.1,0,0,0))

for (s in (1:J)[-4]) {
  plot(signal_val[[s]][,1], xlab ="",ylab="",
       main = paste("Facility",LETTERS[s]),
       pch=21,bg=2,cex=0.9,lwd=0.3)
  mtext(TeX("Signal Value"),side=2,line=1.9,cex=0.6)
  mtext("time (in weeks)",side=1,line=1.9,cex=0.6)
  par(new=T)
  plot(X2[[s]],cex=0.8, pch = 21, xlab = "", col="blue",
       ylab = "", axes = F,lwd=0.8)
  axis(4, at = pretty(range((X2[[s]]))))
  mtext(TeX("N1/PMMoV$\\times 10^3$"),side=4,line=1.9,cex=0.51)
}
plot(X2[[1]],xlab = "",ylab="",axes=F,col="white")
legend(-5,9, pch = 21,pt.bg = c(2,"white"), col = c(1,"blue"),
       bty="n",cex=1.3, inset = c(2,0),
       legend = c("Virus Concetration signal",TeX("N1/PMMoV$\\times 10^3$")))
dev.off()

#-------------------------------------------------------------------------------

pdf("negBin_probs_to_Z_ts.pdf",width=6,height=8)
par(mar = c(2.8,3.1,1.9,3), mfrow = c(5,3), oma = c(0.05,0,0,0))

for (i in 1:J) {
  plot(negBinProbs_est[[i]][,1], xlab = "", ylab = "", pch = 21,
       main = paste("Facility",LETTERS[i]), cex=0.9, bg=2,lwd=0.3,
       ylim = c(0,1))
  mtext(TeX("$p_{ts}$"),side=2,line=1.9,cex=0.8)
  mtext("time (in weeks)",side=1,line=1.9,cex=0.58)
  par(new = T)
  plot(Z_ts[[i]],xlab = "",ylab = "", axes = F,bg = "grey",
       cex=0.85,lwd=0.7,pch=21)
  axis(4, at = pretty(range(Z_ts[[i]]))) 
  mtext(TeX("Reported Cases"),side = 4, line=1.9,cex = 0.52)
}

plot(negBinProbs_est[[14]][,1], xlab="", ylab="", col="white", axes = F)
legend(-3,1, pch = c(21,21), col = 1, pt.bg = c(2,"grey"),bty = "n",
       cex=1.3, legend = TeX(c("Mean $p_{ts}$","Reported Cases")))

dev.off()

#-------------------------------------------------------------------------------

pdf("reported_cases_to_N1_PMMoV.pdf",width = 6.5,height = 8)
par(mar = c(3,3.1,1.9,3.1), mfrow = c(5,3), oma = c(0.1,0,0,0))
for (i in 1:J) {
  plot(X2[[i]], pch = 21, bg = "blue", cex = 0.9,
       main = paste("Facility",LETTERS[i]),xlab="",ylab="")
  mtext("time (in weeks)",side=1,line=2,cex=0.68)
  mtext(TeX("N1/PMMoV$\\,\\times 10^3$"),side=2,line=2,cex=0.55)
  par(new=T)
  plot(Z_ts[[i]],xlab="",ylab="",cex=0.87,pch=21,bg="grey",axes=F)
  axis(4, at = pretty(range(Z_ts[[i]])))
  mtext(TeX("Reported Cases"),side=4,cex=0.62,line=2)
}

plot(X2[[1]],xlab = "",ylab="",axes=F,col="white")
legend("topleft", pch = 21, pt.bg = c("blue","grey"), bty="n",cex=1.15,
       legend = TeX(c("N1/PMMoV$\\,\\times 10^3$","Reported Cases")))
dev.off()


