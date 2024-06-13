
###############################################################################
#  Projections of wastewater as an indicator of COVID-19 cases in corrections facilities: a modelling study          #
#																	      	  #
#          Copyright(2021): Dan Han, Pamela Linares, University of Louisville, KY, USA
#																			  #	
#                         Version 1 - June 13, 2024

# In this code,
# ln(p_ts/(1-p_ts))=Virus+gamma(t,s) where gamma(t,s)=a1s+a2s*gamma(t-1,s)+e(t,s)
# that is ln(p_ts/(1-p_ts))=Virus+a1s+a2s*gamma(t-1,s)

###############################################################################





library(truncnorm)
library(BayesLogit)

#------------------------------------------------------------------------------
# P0 is the covariance matrix of the prior of beta
# m0 is the mean vector of the prior of beta
# mu_alpha_0 is the vector of means for the prior of alpha_s
# mu_alpha_0[s] is the mean of the prior of alpha_s
# Sigma_alpha_0[s] is the variance of the prior of alpha_s for s locations 

# mu_a_0 is a vector of the means of the prior of a_s for s locations 
# mu_a_0[s] is the mean of the prior of a_s
# sigma_a_0 is a vector of the variances of the prior of a_s for s locations 
# sigma_a_0[s] is the variance of the prior of a_s

# mu_epsilon_0 is a vector of means for the prior of epsilon_s for s locations 
# mu_epsilon_0[s] is the mean of the prior of epsilon_s
# Sigma_epsilon_0[s]=sigma_epsilon_0[s]*diag(N) the variance of the prior of epsilon_s for s locations 

model.R <- function(Z_t, X_1 = X_1,
                          X_2, X_population,
                          location = location,
                          h = h, m0 = m0, P_0 = P0,
                          sigma_epsilon_0, mu_epsilon_0,
                          sigma_e_0, mu_a_0 = mu_a_0,
                          sigma_a_0 = sigma_a_0,
                          main.iters = 1000, burn = 500,
                          verbose = 500) {
  
  # the number of categories/locations
  J <- length(unique(location))
  
  # the number of parameters in X_1 after excluding the location/id
  p_1 <- ncol(X_1)-1
  
  # the number of parameters in X_2 after excluding the location/id
  p_2 <- ncol(X_2)-1
  
  total.iters <- main.iters + burn
  beta <- array(0.0, dim=c(p_1,J,total.iters))
  
  # Timing
  start.time <- Sys.time()
  
  X_1_list <- list()
  X_2_list <- list()
  X_population_list <- list()
  
  Z_t_list <- list()
  Ns_list <- list()
  N_ts_list <- list()
  
  epsilon_list <- list()
  gamma_list <- list()
  gamma_s_list <- list()
  
  lambda_list <- list()
  
  e_list <- list()
  a_list <- list()
  
  reporting_probs <- list()
  negative_binomial_probs <- list()

  for (s in unique(location)){
    
    # Design matrix containing factors that influence the report prob in sth location
    X_1_list[[s]] <- as.matrix(X_1[X_1$location == s,-1])
    
    # Design matrix containing factors that influence true covid cases in sth location
    X_2_list[[s]] <- as.matrix(X_2[X_2$location == s,-1])
    
    # The population at time t, location s
    X_population_list[[s]] <- as.matrix(X_population[X_population$location == s,-1])
    
    # The recorded inmate positive cases number in sth location
    Z_t_list[[s]] <- as.matrix(Z_t[Z_t$location == s,-1])
    
    # Ns is the total number of samples/time points at location s
    Ns_list[[s]] <- dim(X_1_list[[s]])[1]
    
    # N_ts is the true positive cases number for time t at location s
    N_ts_0 <- array(0, dim=c(Ns_list[[s]],total.iters))
    N_ts_0[,1] <- Z_t_list[[s]]
    N_ts_list[[s]] <- N_ts_0
    
    # create empty matrix for other parameters
    gamma_list[[s]] <- array(0.0, dim=c(Ns_list[[s]],total.iters))
    gamma_s_list[[s]] <- array(0.0, dim=c(Ns_list[[s]],total.iters))
    
    lambda_list[[s]] <- array(0.0, dim=c(Ns_list[[s]],total.iters))
    
    e_list[[s]] <- array(0.0, dim=c(Ns_list[[s]],total.iters))
    a_list[[s]] <- array(0.0, total.iters)
    
    epsilon_list[[s]] <- array(0.0, dim=c(Ns_list[[s]],total.iters))
    
    reporting_probs[[s]] <- array(0.0, dim=c(Ns_list[[s]],total.iters))
    negative_binomial_probs[[s]] <- array(0.0, dim=c(Ns_list[[s]],total.iters))
  }
  
  # First iteration of epsilon, e_ts, and gamma
  for (s in 1:J){
    
    for (i in 1:Ns_list[[s]]) {
      epsilon_list[[s]][i,1] <- mu_epsilon_0[s]+sqrt(sigma_epsilon_0[s])%*%rnorm(1)
      e_list[[s]][i,1] <- rnorm(1, mean = 0, sd = sqrt(sigma_e_0[s]))
    }
    
    for (i in 2:Ns_list[[s]]) {
      gamma_list[[s]][i,1] <- a_list[[s]][1]*gamma_list[[s]][i-1,1]+e_list[[s]][i,1]
      gamma_s_list[[s]][i,1] <- gamma_list[[s]][i-1,1] 
    }
  }
  
  
  message("> Begin Sampling")
  
  for (k in 2:total.iters) {
    
    if (k == burn + 1) start.ess <- Sys.time()
    if (k %% 1000 == 0 & k < total.iters) cat("\niterations",k,"to",999+k)
    
    for (s in 1:J){
      
      Z_ts <- as.matrix(Z_t_list[[s]])
      N_ts <- as.matrix(N_ts_list[[s]][,k-1])
      X_1s <- as.matrix(X_1_list[[s]])
      X_2s <- as.matrix(X_2_list[[s]])
      Ns <- Ns_list[[s]]
      kappa1 <- Z_ts-N_ts/2
      
      # Drawing w_1s
      psi_ts <- X_1s%*%beta[,s,k-1]+epsilon_list[[s]][,k-1]
      w_s1 <- rpg.devroye(Ns,N_ts, psi_ts)
      
      # Get the invert matrix of V_t.
      V_t_invert <- t(X_1s)%*%(diag(w_s1)+sigma_epsilon_0[s]^{-1}*diag(Ns))%*%X_1s+solve(P_0)
    
      # Get V_j, the covariance matrix of beta_j
      V_t <- chol2inv(chol(V_t_invert))
      
      # Get m_j, the mean of beta_t
      m_t <- V_t%*%(t(X_1s)%*%(kappa1-diag(w_s1)%*%epsilon_list[[s]][,k-1]) +
                      sigma_epsilon_0[s]^{-1}*t(X_1s)%*%
                      (psi_ts-mu_epsilon_0[s]*rep(1,Ns))+chol2inv(chol(P_0))%*%m0)
      
      # Sample the vector beta_t for location s
      beta_t <-  m_t + t(chol(V_t))%*%rnorm(p_1)
      beta[,s,k] <- beta_t
      
      # Sample epsilon_s
      for (i in 1:Ns) {
        sigma_epsilon_is <- 1/(w_s1[i]+1/sigma_epsilon_0[s])
        mu_epsilon_is <- (kappa1[i] - w_s1[i]*X_1s[i,]%*%beta_t +
                            mu_epsilon_0[s]/sigma_epsilon_0[s])*sigma_epsilon_is
        
        # Sample the error vector epsilon_s for s locations
        epsilon_list[[s]][i,k] <- rnorm(1, mu_epsilon_is, sqrt(sigma_epsilon_is))
      }
      
      # Calculate the reporting probabilities
      log_p_ts <- X_1s%*%beta_t+epsilon_list[[s]][,k]-log(1+exp(X_1s%*%beta_t+epsilon_list[[s]][,k]))
      p_ts <- exp(log_p_ts)
      reporting_probs[[s]][,k] <- p_ts
      
      ################################################################################################
      
      kappa2 <- (N_ts-h)/2
      
      # Drawing w_2
      phi_ts <- X_2s+gamma_list[[s]][,k-1]
      w_s2 <- rpg.devroye(Ns, N_ts+h, phi_ts)
      
      
      # Sample a_s
      V_a_s_invert <- t(gamma_s_list[[s]][,k-1])%*%diag(w_s2)%*%gamma_s_list[[s]][,k-1]+
                      t(gamma_s_list[[s]][,k-1])%*%diag(sigma_e_0[s]^{-1},Ns)%*%
                      gamma_s_list[[s]][,k-1]+sigma_a_0[s]^{-1}
      V_a_s <- chol2inv(chol(V_a_s_invert))
      
      mu_a_s <- V_a_s%*%(t(gamma_s_list[[s]][,k-1])%*%diag(w_s2)%*%(kappa2/w_s2-X_2s-e_list[[s]][,k-1])+
                           t(gamma_s_list[[s]][,k-1])%*%diag(sigma_e_0[s]^{-1},Ns)%*%(phi_ts-X_2s)+
                           sigma_a_0[s]^{-1}*mu_a_0[s])
      
      a_list[[s]][k] <- rtruncnorm(1, a=-1, b=1, mean = mu_a_s, sd = sqrt(V_a_s))
      a_s <- a_list[[s]][k]
      
      
      # sample e_s
      for (i in 1:Ns) {
        sigma_e_is <- 1/(w_s2[i]+1/sigma_e_0[s])
        mu_e_is <- (kappa2[i]-w_s2[i]*X_2s[i]-w_s2[i]*gamma_s_list[[s]][i,k-1]%*%a_s)*sigma_e_is
        
        # sample the error vector e_s for the sth location
        e_list[[s]][i,k] <- rnorm(1, mean=mu_e_is, sd=sqrt(sigma_e_is))
      }
      
      for (i in 2:Ns_list[[s]]) {
        gamma_list[[s]][i,k] <- a_s*gamma_list[[s]][i-1,k-1]+e_list[[s]][i,k]
        gamma_s_list[[s]][i,k] <- gamma_list[[s]][i-1,k]
      }
      
      
      # Calculate N_ts negative binomial probabilities
      log_p_ts_nb <- X_2s+gamma_list[[s]][,k]-log(1+exp(X_2s+gamma_list[[s]][,k]))
      p_ts_nb <- exp(log_p_ts_nb)
      
      if (any(is.infinite(exp(X_2s+gamma_list[[s]][,k])))) {
        p_ts_nb[which(is.infinite(exp(X_2s+gamma_list[[s]][,k])))] <- 1
      }
      
      negative_binomial_probs[[s]][,k] <- p_ts_nb
      
      # Sampling lambda_ts
      lambda <- rgamma(Ns, shape = N_ts+h, scale = p_ts_nb)
      lambda_list[[s]][,k] <- lambda

      
      # Sampling N_ts, the true number of positive cases
      for (i in 1:Ns) {
        N_ts_list[[s]][i,k] <- min(rpois(1,lambda=(1-p_ts[i])*lambda[i])+Z_ts[i],
                                   X_population_list[[s]][i]) 
        } 
    }
  }
  
  end.time <- Sys.time()
  
  total.time <- difftime(end.time, start.time)
  ess.time <- difftime(end.time, start.ess)
  
  output <- list(Z_ts = Z_t_list, N_ts = N_ts_list, 
                 reporting_probs = reporting_probs,
                 negative_binomial_probs = negative_binomial_probs,
                 beta = beta, a = a_list, gamma = gamma_s_list, 
                 lambda = lambda_list, epsilon = epsilon_list, e_ts = e_list, 
                 X_1 = X_1_list, X_2 = X_2_list, time_points = Ns_list,
                 total.time = total.time, ess.time = ess.time)
  
  return(output) 
  
}
