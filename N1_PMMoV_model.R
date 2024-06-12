
source("underreport_final.R")

#-------------------------------------------------------------------------------
# Loading Data 

full_data <- read.csv(paste0("full_wkly_data.csv"))
full_data <- na.omit(full_data) 

location <- full_data$id
J <- length(unique(location))

population <- full_data$in_pop
capacity_ratio <- full_data$capacity_ratio
N1_PMMoV <- full_data$N1_PMMoV
em_cases_ratio <- full_data$em_cases/full_data$capacity

#-------------------------------------------------------------------------------
# Model

# Design matrix containing factors that influence the report probability
X_1 <- data.frame(location,capacity_ratio,em_cases_ratio)
p_1 <- ncol(X_1)-1

# Design matrix containing factors that influence the true positive covid cases
X_2 <- data.frame(location,N1_PMMoV=N1_PMMoV*10^3)
p_2 <- ncol(X_2)-1

# Matrix containing the inmate population of each location at time t
X_pop <- data.frame(location,population)

# Reported positive cases
Z_t <- data.frame(location, reported = full_data$in_cases)

# Mean and variance for epsilon's prior, for location s=1,2,3,...J
mu_epsilon_0 <- rep(0,J)
sigma_epsilon_0 <- rep(1,J)

# variance for e's prior, for location s=1,2,3,...J
sigma_e_0 <- rep(1,J)

# Mean and variance of a's prior for location s=1,2,3,...J
mu_a_0 <- rep(0,J)
sigma_a_0 <- rep(1000,J)

h <- 1

# m0 is the mean vector of the prior of beta_s
# P_0 is the covariance matrix of the prior of beta_s
m0 <- rep(0,p_1)
P_0 <- diag(1,p_1)

main.iters <- 10000
burn <- 30000
total <- main.iters+burn

set.seed(313)
results <- underreport.R(Z_t, X_1 = X_1 ,X_2 = X_2,
                         X_population = X_pop, 
                         location = location, 
                         h = h, m0 = m0, P_0 = P_0,
                         mu_epsilon_0 = mu_epsilon_0,
                         sigma_epsilon_0 = sigma_epsilon_0,
                         mu_a_0 = mu_a_0,
                         sigma_a_0 = sigma_a_0,
                         sigma_e_0 = sigma_e_0,
                         main.iters = main.iters,
                         burn = burn, verbose=500)

rm(mu_a_0,mu_epsilon_0,sigma_a_0,Z_t,sigma_e_0,m0,
   P_0,em_cases_ratio,sigma_epsilon_0,main.iters,h,
   X_1,X_2,X_pop,population,capacity_ratio,p_1,p_2)

gc()


################################################################################
# a_s and beta estimates

library(matrixStats)

burn <- burn+1
quantiles <- c(0.05, 0.95)

beta_list <- list()
beta1_conf <- list()   ;   beta2_conf <- list()
beta1_est <- list()    ;   beta2_est <- list()

a_s <- list()
a_est <- list()   ;   a_conf <- list()

for (i in 1:J) {
  s <- unique(location)[i]
  beta_list[[s]] <- results$beta[,i,burn:total]
  
  beta1_conf[[s]] <- quantile(beta_list[[s]],quantiles)
  beta2_conf[[s]] <- quantile(beta_list[[s]][2,],quantiles)
  
  beta1_est[[s]] <- c("mean"=mean(beta_list[[s]]),
                      "med"=median(beta_list[[s]]))
  
  beta2_est[[s]] <- c("mean"=mean(beta_list[[s]][2,]),
                      "med"=median(beta_list[[s]][2,]))
  
  a_s[[s]] <- results$a[[s]][burn:total]
  a_est[[s]] <- c("mean"=mean(a_s[[s]]),"median"=median(a_s[[s]]))
  a_conf[[s]] <- quantile(a_s[[s]],quantiles)
}

beta1_est
beta2_est
a_est

#-------------------------------------------------------------------------------
# All other estimates

Z_ts <- results$Z_t
X2 <- results$X_2

N_ts <- list()
N_t_est <- list()
N_ts_conf <- list()

gamma <- list()
gamma_est <- list()

report_probs <- list()
negBin_probs <- list()
reportProbs_est <- list()
negBinProbs_est <- list()

for (s in unique(location)) {
  N_ts[[s]] <- results$N_ts[[s]][, burn:total]
  gamma[[s]] <- results$gamma[[s]][, burn:total]
  
  report_probs[[s]] <- results$reporting_probs[[s]][ ,burn:total]
  negBin_probs[[s]] <- results$negative_binomial_probs[[s]][ ,burn:total]
  
  N_t_est[[s]] <- cbind("mean"=apply(N_ts[[s]],1, mean),
                        "median"=apply(N_ts[[s]],1, median))
  
  N_ts_conf[[s]] <- apply(N_ts[[s]],1,quantile,quantiles)
  
  gamma_est[[s]] <- cbind("mean"=apply(gamma[[s]], 1, mean),
                          "median"=apply(gamma[[s]], 1, median))
  
  reportProbs_est[[s]] <- cbind("mean"=apply(report_probs[[s]], 1, mean),
                                "median"=apply(report_probs[[s]], 1, median))
  
  negBinProbs_est[[s]] <- cbind("mean"=apply(negBin_probs[[s]], 1, mean),
                                "median"=apply(negBin_probs[[s]], 1, median))
}


################################################################################
####          Finding values of N1/PMMoV that predict a covid case          ####
################################################################################

# Method 1: Find p_ts > 0.8

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
        above_thresh[[s]][[i]][-which(X2[[s]][above_thresh[[s]][[i]]] == 0)]
    }
    
    N1_PMMoV_ts[[s]][[i]] <- X2[[s]][above_thresh[[s]][[i]]]/1e3
    
    N1_PMMoV_s[[s]][[i]] <- c(min=min(N1_PMMoV_ts[[s]][[i]]),
                              median=median(N1_PMMoV_ts[[s]][[i]]),
                              mean=mean(N1_PMMoV_ts[[s]][[i]]),
                              p_ts=mean(negBinProbs_est[[s]][above_thresh[[s]][[i]],1]))
  }
}


# quantiles <- c(0.025,0.25,0.5,0.75,0.975)
# N1_PMMoV_vals <- unlist(N1_PMMoV_ts)
# N1_PMMoV_state <- quantile(N1_PMMoV_vals, quantiles)
# (N1_PMMoV_state <- c(N1_PMMoV_state,"mean"=mean(N1_PMMoV_vals)))

N1_PMMoV_vals <- list()
for (i in 1:6 ) {
  N1_PMMoV_vals[[i]] <- data.frame(Min = 0,Median = 0, Mean = 0, p_ts = 0)
  for (j in 1:14) {
    N1_PMMoV_vals[[i]][j,] <- N1_PMMoV_s[[j]][[i]]
  }
  rownames(N1_PMMoV_vals[[i]]) <- unique(location)
}

# Saving Results

sink("N1_PMMoV_vals_diff_thresh.txt")

cat("# Values of N1/PMMoV that can predict one case with different probability #\n")
options(digits = 4,scipen = -2)

cat("\n# Values of N1/PMMoV such that p_ts >=",paste0(thresh[1]*100,"%"),"\n\n")
N1_PMMoV_vals[[1]]
cat("\n")

cat("\n# Values of N1/PMMoV such that p_ts >=",paste0(thresh[2]*100,"%"),"\n\n")
N1_PMMoV_vals[[2]]
cat("\n")

cat("\n# Values of N1/PMMoV such that p_ts >=",paste0(thresh[3]*100,"%"),"\n\n")
N1_PMMoV_vals[[3]]
cat("\n")

cat("\n# Values of N1/PMMoV such that p_ts >=",paste0(thresh[4]*100,"%"),"\n\n")
N1_PMMoV_vals[[4]]
cat("\n")

cat("\n# Values of N1/PMMoV such that p_ts >=",paste0(thresh[5]*100,"%"),"\n\n")
N1_PMMoV_vals[[5]]
cat("\n")

cat("\n# Values of N1/PMMoV such that p_ts >=",paste0(thresh[6]*100,"%"),"\n\n")
N1_PMMoV_vals[[6]]
cat("\n\n")

sink()

#-------------------------------------------------------------------------------

# Finding the min N1/PMMoV under all thresholds of p_ts
N1_PMMoV_mins <- list()

for (i in unique(location)) {
  N1_PMMoV_mins[[i]] <- N1_PMMoV_s[[i]][[1]]["min"]
  for (k in 2:6) {
    N1_PMMoV_mins[[i]] <- append(N1_PMMoV_mins[[i]],N1_PMMoV_s[[i]][[k]]["min"])
  }
}

#-------------------------------------------------------------------------------
# Finding the capture ratio under all p_ts thresholds 

sens_tab <- list()
sens_sum <- list()

library(dplyr)

Z_ts <- list()

for (s in unique(location)) {
  Z_ts[[s]] <- full_data$in_cases[full_data$id == s]
}

for (s in unique(location)) {
  
  sens_tab[[s]] <- list()
  sens_sum[[s]] <- list()
  
  
  for (i in 1:length(thresh)) {
    
    min <- N1_PMMoV_s[[s]][[i]]["min"]
    sens_tab[[s]][[i]] <- data.frame(Z_ts_ge_0 = ifelse(Z_ts[[s]]==0,0,1),
                                     N1_PMMoV_ge_min = ifelse(X2[[s]] < min*1e3, 0, 1),
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
  sens_sum[[s]] <- cbind(thresh,
                         "ratio"=sens_sum[[s]][,"TP"]/sens_sum[[s]][,"Z_ts_ge_0"])
}

#---------------------------------------------------------------------------------------

sink("ROC_table.txt")

cat("\n# Overall count at the 6 different thresholds for the 14 location #\n")

cat("> sens_sum\n")
sens_sum
cat("\n")

sink()
#-------------------------------------------------------------------------------
# Method 2: Using the derivation formula

signal_val <- list()

for (s in unique((location))) {
  signal_val[[s]] <- matrix(NA,ncol=2,nrow=results$time_points[[s]])
  colnames(signal_val[[s]]) <- colnames(gamma_est[[s]])
  
  for (i in 1:results$time_points[[s]]) {
    signal_val[[s]][i,1] <- max(0,log(4) - gamma_est[[s]][i,1])
    signal_val[[s]][i,2] <- max(0,log(4) - gamma_est[[s]][i,2])
  }
}

save.image()

