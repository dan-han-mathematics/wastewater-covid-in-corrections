# This files generate paper tables

#-------------------------------------------------------------------------------

# Table of a_s estimates
a_table <- bind_rows(a_est) %>%
  bind_cols(bind_rows(a_conf)) %>%
  mutate(location = paste(LETTERS[1:14],names(X2),sep ="_")) %>%
  select(location,everything())

write_csv(a_table,"a_table.csv")


# Table of beta coefficients
beta_table <- bind_cols(location=paste(LETTERS[1:14],names(X2),sep ="_"),
                        coef=rep("$\\beta_{1}$",14),bind_rows(beta1_est)) %>%
  bind_rows(bind_cols(location=paste(LETTERS[1:14],names(X2),sep ="_"),
                      coef=rep("$\\beta_{2}$",14),bind_rows(beta2_est))) %>%
  left_join(.,bind_cols(location=paste(rep(LETTERS[1:14],2),
                                       rep(names(beta2_est),2),sep = "_"),
                        coef=rep(c("$\\beta_{1}$","$\\beta_{2}$"),each=14),
                        bind_rows(bind_rows(beta1_conf),bind_rows(beta2_conf))),
            by = c("location", "coef")) %>%
  arrange(location,coef)

write_csv(beta_table,"beta_table.csv")


#-------------------------------------------------------------------------------

# Table of N1/PMMoV inconsistencies with Z_ts 
max_p_ts1 <- NULL
for (s in unique(location)) {
  max_p_ts1[s] <- negBinProbs_est[[s]][which.max(X2[[s]]),"mean"]
}

max_N1_PMMoV <- full_data %>%
  summarise(week = week[which.max(N1_PMMoV)],
            week_num = which.max(N1_PMMoV),
            max_N1_PMMoV = N1_PMMoV[which.max(N1_PMMoV)], 
            Z_ts = in_cases[which.max(N1_PMMoV)],.by = "id") %>%
  mutate(mean_p_ts = max_p_ts1, Facility = LETTERS[1:14]) %>%
  unite(name,c("Facility","id")) %>%
  unite(week, week:week_num, sep = " (")

write_csv(max_N1_PMMoV,"./max_N1_PMMoV_to_Z_ts_and_p_ts.csv")


# Table of N1/PMMoV inconsistencies with Z_ts 
max_p_ts2 <- NULL
for (s in unique(location)) {
  max_p_ts2[s] <- negBinProbs_est[[s]][which.max(Z_ts[[s]]),"mean"]
}

max_Z_ts <- full_data %>%
  summarise(week = week[which.max(in_cases)], 
            week_num = which.max(in_cases),
            max_Z_ts = in_cases[which.max(in_cases)],
            N1_PMMoV = N1_PMMoV[which.max(in_cases)],.by ="id") %>%
  mutate(mean_p_ts = max_p_ts2, Facility = LETTERS[1:14]) %>%
  unite(name,c("Facility","id")) %>%
  unite(week, week:week_num, sep = " (")

write_csv(max_Z_ts,"./max_Z_ts_to_N1_PMMoV_and_p_ts.csv")

  



              
