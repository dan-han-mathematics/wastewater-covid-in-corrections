

#---------------------------------------------------------------------------------------------
# N1/PMMoV values at various thresholds

thresh <- c(0.5,0.6,0.7,0.8,0.9,0.95)

above_thresh <- list()

N1_PMMoV_ts <- list()
N1_PMMoV_s  <- list()

for (s in unique(location)) {
  
  above_thresh[[s]] <- list()
  N1_PMMoV_ts[[s]] <- list()
  N1_PMMoV_s[[s]] <- list()
  
  for (i in 1:length(thresh)){
    above_thresh[[s]][[i]] <- sort(unique(c(which(negBinProbs_est[[s]][,"median"] >= thresh[i]),
                                            which(negBinProbs_est[[s]][,"mean"] >= thresh[i]))))
    
    if (any(X2[[s]][above_thresh[[s]][[i]]] == 0)) {
      above_thresh[[s]][[i]] <- 
        above_thresh[[s]][[i]][-which(X2[[s]][above_thresh[[s]][[i]]] <= 0)]
    }
    
    N1_PMMoV_ts[[s]][[i]] <- X2[[s]][above_thresh[[s]][[i]]]*1e3
    
    N1_PMMoV_s[[s]][[i]] <- c(min=min(N1_PMMoV_ts[[s]][[i]]),
                              median=median(N1_PMMoV_ts[[s]][[i]]),
                              mean=mean(N1_PMMoV_ts[[s]][[i]]),
                              p_ts=mean(negBinProbs_est[[s]][above_thresh[[s]][[i]],1]))
  }
}

#---------------------------------------------------------------------------------------------
# Sensitivity tables

sens_tab <- list()
sens_sum <- list()

library(dplyr)

for (s in unique(location)) {
  
  sens_tab[[s]] <- list()
  sens_sum[[s]] <- list()
  
  
  for (i in 1:length(thresh)) {
    
    min <- N1_PMMoV_s[[s]][[i]]["min"]
    sens_tab[[s]][[i]] <- data.frame(Z_ts_ge_0 = ifelse(results$Z_ts[[s]]==0,0,1),
                                     N1_PMMoV_ge_min = ifelse(X2[[s]]*1e3 < min, 0, 1),
                                     TP = 0, FP = 0, TN = 0, FN = 0)
    
  
    for (k in 1:nrow(sens_tab[[s]][[i]])) {
      if (sens_tab[[s]][[i]]$Z_ts_ge_0[k]==1 & sens_tab[[s]][[i]]$N1_PMMoV_ge_min[k]==1) {
        sens_tab[[s]][[i]]$TP[k] = 1
      } else if (sens_tab[[s]][[i]]$Z_ts_ge_0[k]==0 & sens_tab[[s]][[i]]$N1_PMMoV_ge_min[k]==1) {
        sens_tab[[s]][[i]]$FP[k] = 1
      } else if (sens_tab[[s]][[i]]$Z_ts_ge_0[k]==1 & sens_tab[[s]][[i]]$N1_PMMoV_ge_min[k]==0) {
        sens_tab[[s]][[i]]$FN[k] = 1
      } else if (sens_tab[[s]][[i]]$Z_ts_ge_0[k]==0 & sens_tab[[s]][[i]]$N1_PMMoV_ge_min[k]==0) {
        sens_tab[[s]][[i]]$TN[k] = 1
      }
    }
  }
  
  

  sens_sum[[s]] <- cbind(thresh,
                         rbind(colSums(sens_tab[[s]][[1]]),
                         colSums(sens_tab[[s]][[2]]),
                         colSums(sens_tab[[s]][[3]]),
                         colSums(sens_tab[[s]][[4]]),
                         colSums(sens_tab[[s]][[5]]),
                         colSums(sens_tab[[s]][[6]])))
  sens_sum[[s]] <- cbind(sens_sum[[s]],
                         "ratio"=sens_sum[[s]][,"TP"]/sens_sum[[s]][,"Z_ts_ge_0"])
}


#---------------------------------------------------------------------------------------

sink("ROC_method_N1_PMMoV.txt")

cat("\n# Overall count at the 6 different thresholds for the 14 location #\n")

cat("> sens_sum\n")
sens_sum

sink()

df.two.locations=data.frame(Threshold=sens_sum$BCFC[,1], A=sens_sum$BCFC[,8],B=sens_sum$BCC[,8],C=sens_sum$EKCC[,8],D=sens_sum$GRCC[,8],E=sens_sum$KCIW[,8],F=sens_sum$KSP[,8],G=sens_sum$KSR[,8],H=sens_sum$LAC[,8],I=sens_sum$LSCC[,8],
                            J=sens_sum$LLCC[,8],K=sens_sum$NTC[,8], L=sens_sum$RCC[,8],M=sens_sum$SSCC[,8],N=sens_sum$WKCC[,8])
df.two.locations=data.frame(Threshold=sens_sum$BCFC[,1], A=sens_sum$BCFC[,8],B=sens_sum$BCC[,8],D=sens_sum$GRCC[,8],E=sens_sum$KCIW[,8],L=sens_sum$RCC[,8])

df.more.than.80=data.frame(Threshold=sens_sum$BCFC[,1], A=sens_sum$BCFC[,8],C=sens_sum$EKCC[,8],E=sens_sum$KCIW[,8],F=sens_sum$KSP[,8],G=sens_sum$KSR[,8],I=sens_sum$LSCC[,8],H=sens_sum$LAC[,8], J=sens_sum$LLCC[,8])
df.less.than.80=data.frame(Threshold=sens_sum$BCFC[,1], B=sens_sum$BCC[,8],D=sens_sum$GRCC[,8],K=sens_sum$NTC[,8], L=sens_sum$RCC[,8],M=sens_sum$SSCC[,8],N=sens_sum$WKCC[,8])





library(MASS) 
library(reshape2) 
df.more.than.80_melt=melt(df.more.than.80,id=c("Threshold"))
df.less.than.80_melt=melt(df.less.than.80,id=c("Threshold"))
colnames(df.more.than.80_melt) <- c('Threshold','Facility','Ratio') 
# install.packages('ggsci')
library(ggsci)  
library(ggplot2)  
library(latex2exp)
p1 = ggplot(data = df.more.than.80_melt, aes(x = Threshold, y =Ratio,color=Facility,shape = Facility)) + 
  geom_point(size=6) +
  geom_line(size=0.5)+
  xlab(TeX('Thresholds of $p_{ts}$'))+
  ylab('Capture Ratio')+
  theme_classic()+
  xlim(0.5, 0.8)+ylim(0.4, 1)+
  scale_color_lancet()+scale_fill_lancet()

print(p1)

p2 = ggplot(data = df.less.than.80_melt, aes(x = Threshold, y =Ratio,color=Facility,shape = Facility)) + 
  geom_point(size=2) +
  geom_line()+
  xlab(TeX('Thresholds of $p_{ts}$'))+
  ylab('Capture Ratio')+
  theme_classic()+
  xlim(0.5, 0.8)+ylim(0.4, 1)+
  scale_color_lancet()+scale_fill_lancet()

print(p2)


above0.8 <- list()

week_s1 <- list()
N1_PMMoV_ts1 <- list()
N1_PMMoV_s1  <- list()

for (s in unique(location)) {
  
  above0.8[[s]] <- sort(unique(c(which(negBinProbs_est[[s]][,"median"] >= 0.8),
                                 which(negBinProbs_est[[s]][,"mean"] >= 0.8))))
  
  if (any(X2[[s]][above0.8[[s]]]== 0)) {
    above0.8[[s]] <- above0.8[[s]][- which(X2[[s]][above0.8[[s]]]<= 0)]     }
  
  week_s1[[s]] <- full_data$week[full_data$id == s][above0.8[[s]]]
  N1_PMMoV_ts1[[s]] <- X2[[s]][above0.8[[s]]]
  N1_PMMoV_s1[[s]] <- c(min=min(N1_PMMoV_ts1[[s]]/1e3),
                        median=median(N1_PMMoV_ts1[[s]]/1e3),
                        mean=mean(N1_PMMoV_ts1[[s]]/1e3))
  mean(negBinProbs_est[[s]][above0.8[[s]],1])
}



pdf("N1_PMMoV_vals_comparison0.8.pdf",width = 7,height = 8.8)
par(mar = c(2,3.1,4,3.1), mfrow = c(2,2), oma = c(0.1,0,0,0))
for (s in c(7,12)) {
  plot(X2[[s]], xlab ="",ylab="",main = paste("Facility",LETTERS[s]),
       lwd=0.9,col="steelblue",pch=21)
  points(above0.8[[s]],N1_PMMoV_ts1[[s]],bg="blue",pch=21,lwd=0.5)
  mtext(TeX("N1/PMMoV$\\times 10^3$"),side=2,line=1.9,cex=0.58)
  mtext("time (in weeks)",side=1,line=1.9,cex=0.62)
  par(new=T)
  plot(above0.8[[s]],negBinProbs_est[[s]][above0.8[[s]],1],
       pch = 21, col = 2, xlab = "", ylab = "", axes = F,
       xlim = c(1,results$time_points[[s]]),ylim=c(0,1))
  axis(4)
  mtext(TeX("$p_{ts}$"), side=4, line=2,cex=0.75)
}
mtext("A", side = 3, line =-1.5, outer = TRUE)
par(mai=c(0,0,0,0))
plot.new()
legend("topleft", pch = c(21,21,21), col = c("steelblue",1,2),
       pt.bg = c("white","blue","white"),cex = 1.2,  bty = "n", 
       legend = TeX(c("N1/PMMoV", "N1/PMMoV confirming at least one case",
                      "$p_{ts} \\geq 0.8$")))
dev.off()



####################################################################



above0.5 <- list()

week_s1 <- list()
N1_PMMoV_ts1 <- list()
N1_PMMoV_s1  <- list()

for (s in unique(location)) {
  
  above0.5[[s]] <- sort(unique(c(which(negBinProbs_est[[s]][,"median"] >= 0.5),
                                 which(negBinProbs_est[[s]][,"mean"] >= 0.5))))
  
  if (any(X2[[s]][above0.5[[s]]]== 0)) {
    above0.5[[s]] <- above0.5[[s]][- which(X2[[s]][above0.5[[s]]]<= 0)]     }
  
  week_s1[[s]] <- full_data$week[full_data$id == s][above0.5[[s]]]
  N1_PMMoV_ts1[[s]] <- X2[[s]][above0.5[[s]]]
  N1_PMMoV_s1[[s]] <- c(min=min(N1_PMMoV_ts1[[s]]/1e3),
                        median=median(N1_PMMoV_ts1[[s]]/1e3),
                        mean=mean(N1_PMMoV_ts1[[s]]/1e3))
  mean(negBinProbs_est[[s]][above0.5[[s]],1])
}



pdf("N1_PMMoV_vals_comparison0.5.pdf",width = 7,height = 8.8)
par(mar = c(2,3.1,4,3.1), mfrow = c(2,2), oma = c(0.1,0,0,0))
for (s in c(7,12)) {
  plot(X2[[s]], xlab ="",ylab="",main = paste("Facility",LETTERS[s]),
       lwd=0.9,col="steelblue",pch=21)
  points(above0.5[[s]],N1_PMMoV_ts1[[s]],bg="blue",pch=21,lwd=0.5)
  mtext(TeX("N1/PMMoV$\\times 10^3$"),side=2,line=1.9,cex=0.58)
  mtext("time (in weeks)",side=1,line=1.9,cex=0.62)
  par(new=T)
  plot(above0.5[[s]],negBinProbs_est[[s]][above0.5[[s]],1],
       pch = 21, col = 2, xlab = "", ylab = "", axes = F,
       xlim = c(1,results$time_points[[s]]),ylim=c(0,1))
  axis(4)
  mtext(TeX("Value of $p_{ts}$"), side=4, line=2,cex=0.6)
}
mtext("B", side = 3, line =-1.5, outer = TRUE)
par(mai=c(0,0,0,0))
plot.new()
legend("topleft", pch = c(21,21,21), col = c("steelblue",1,2),
       pt.bg = c("white","blue","white"),cex = 1.2,  bty = "n", 
       legend = TeX(c("N1/PMMoV", "N1/PMMoV confirming at least one case",
                      "$p_{ts} \\geq 0.5$")))
dev.off()


