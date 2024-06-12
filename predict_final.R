predict.unReport=function(object, new_X1 = NULL,
         new_X2 = NULL, new_pop = NULL, 
         location = location,
         total.iters=10000,
         burn.used = 5000, seed = 123) {
  
  x <- object
  J <- length(unique(location))
  h <- 1
  
  clock.start1 <- Sys.time()
  
  # Saving results parameter estimates
  
  total <- dim(x$beta)[3]
  burn <- burn.used + 1
  
  a_s <- NULL
  beta_est <- list()
  lambda0 <- list()
  gamma <- list()
  
  for (i in 1:J) {
    s <- unique(location)[i]
    beta_est[[s]] <- rowMeans(x$beta[,i,burn:total])

    lambda0[[s]] <- rowMeans(x$lambda[[s]][, burn:total])

    gamma[[s]] <- rowMeans(x$gamma[[s]][, burn:total])


    a_s[s] <- mean(x$a[[s]][burn:total])
  }
  
  
  gamma_list <- list()
  lambda_ts <- list()
  
  Nts_list <- list()
  Z_t_list <- list()
  Ns_list <- list()
  
  X1_list <- list()
  X2_list <- list()
  X_pop_list <- list()
  
  e_list <- list()
  epsilon_list <- list()
  
  negBin_probs <- list()
  report_probs <- list()
  
  total.iters <-total.iters
  
  for (s in unique(location)) {
    
    X1_list[[s]] <- as.matrix(new_X1[new_X1$location == s, -1])
    X2_list[[s]] <- as.matrix(new_X2[new_X2$location == s, -1])
    X_pop_list[[s]] <- as.matrix(new_pop[new_pop$location == s,-1])
    
    Ns_list[[s]] <- dim(X1_list[[s]])[1]
    Nts_list[[s]] <- array(0, dim=c(Ns_list[[s]],total.iters))
    Z_t_list[[s]] <- array(0, dim=c(Ns_list[[s]],total.iters))
    
    e_list[[s]] <- array(0, dim=c(Ns_list[[s]],total.iters))
    epsilon_list[[s]] <- array(0, dim=c(Ns_list[[s]],total.iters))
    
    gamma_list[[s]] <- array(0, dim=c(Ns_list[[s]],total.iters))
    lambda_ts[[s]] <- array(0, dim=c(Ns_list[[s]],total.iters))
    
    report_probs[[s]] <- array(0, dim=c(Ns_list[[s]],total.iters))
    negBin_probs[[s]] <- array(0, dim=c(Ns_list[[s]],total.iters))
  }
  
  set.seed(seed)
  
  # Sample
  message(" Begin Sampling")
  
  for (k in 1:total.iters) {
    
    # if (k == burn.in + 1) 
    clock.start2 <- Sys.time()

    for (s in 1:J) {

      
      Ns <- Ns_list[[s]]
      n <- length(gamma[[s]])

      X_1s <- as.matrix(X1_list[[s]])
      X_2s <- as.matrix(X2_list[[s]])
 
      # Sample epsilon
      epsilon_list[[s]][,k] <- rnorm(Ns)

      # Calculate the reporting probabilities
      log_p_ts <- (X_1s%*%beta_est[[s]]+epsilon_list[[s]][,k]-
                     log(1+exp(X_1s%*%beta_est[[s]]+epsilon_list[[s]][,k])))
      p_ts <- exp(log_p_ts)
      report_probs[[s]][,k] <- p_ts
      
      gamma_list[[s]][1,k] <- a_s[s]*gamma[[s]][n] + rnorm(1,sd=0.1)
      
      for (i in 2:Ns_list[[s]]) {
        gamma_list[[s]][i,k] <- a_s[s]*gamma_list[[s]][i-1,k]+rnorm(1,sd=0.1)
      }
      
      
      # Calculate N_ts's negative binomial probabilities
      log_p_ts_nb <- X_2s+gamma_list[[s]][,k]-log(1+exp(X_2s+gamma_list[[s]][,k]))
      p_ts_nb <- exp(log_p_ts_nb)
      
      negBin_probs[[s]][,k] <- p_ts_nb
      
      lambda <- rgamma(Ns, shape = h, 
                       scale = p_ts_nb/(1-p_ts_nb))
      
      # lambda_ts[[s]][,k] <- c(max(lambda0[[s]][n],lambda[1]), lambda[1:(Ns-1)])
      lambda_ts[[s]][,k] <- c(max(lambda0[[s]][n],lambda[1]), lambda[2:(Ns)])
      
      Z_ts=rpois(Ns,lambda=p_ts*lambda_ts[[s]][,k])

   
      Nts_list[[s]][,k] <- rnbinom(Ns,size=h,prob = (1-p_ts_nb)/(1-p_ts_nb*p_ts) )+Z_ts

      
      
    }
    
  }
 
  
  clock.end <- Sys.time()
  
  total.time <- difftime(clock.end, clock.start1)
  ess.time <- difftime(clock.end,clock.start2)
  
  out <- list(N_ts = Nts_list, Z_ts = Z_t_list, X_2 = X2_list,
              report_probs = report_probs,negBin_probs = negBin_probs,
              beta = beta_est, a_s = a_s, gamma = gamma_list, 
              e_ts = e_list, epsilon = epsilon_list, lambda = lambda_ts,
              total.time = total.time, ess.time = ess.time)
  
  return(out)
  
}
