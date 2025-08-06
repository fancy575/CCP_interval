### Function to generate real time from piecewise exponential distribution
piece_exp_gen <- function(nsim=1000,beta=c(1,1,-1,0.5,-0.5,0.1,-0.1,1,-1,0.5)
                          ,seed=123, reffect = rep(0,1000),
                          interval = c(0), ## pre-defined time intervals
                          lambda = c(exp(1))## pre-defined baseline hazard function 
){
  
  
  set.seed(seed)
  X13 <- sample(x=c(1:4),size=nsim,replace = TRUE,prob = c(0.3,0.2,0.2,0.3))
  
  Z <- rbinom(nsim,1,2/3)
  
  X1 <- ifelse(X13==2,1,0)
  X2 <- ifelse(X13==3,1,0)
  X3 <- ifelse(X13==4,1,0)
  
  X4 <- rnorm(nsim,0,0.5^2)
  X5 <- rnorm(nsim,0,1)
  X6 <- rnorm(nsim,0,1.5^2)
  X7 <- rbinom(nsim,1,0.3)
  X8 <- rbinom(nsim,1,0.4)
  X9 <- rbinom(nsim,1,0.5)
  
  
  exp_const <- c(exp(matrix(beta,ncol=length(beta)) %*% 
                       rbind(matrix(Z,nrow=1),
                             matrix(X1,nrow=1),
                             matrix(X2,nrow=1),
                             matrix(X3,nrow=1),
                             matrix(X4,nrow=1),
                             matrix(X5,nrow=1),
                             matrix(X6,nrow=1),
                             matrix(X7,nrow=1),
                             matrix(X8,nrow=1),
                             matrix(X9,nrow=1)) + reffect ))
  
  J <- length(interval)
  if(J != length(lambda)){
    return ("number of intervals does not match with number of baseline hazard")
  }
  
  diff_interval <- diff(interval)
  
  tau <- runif(nsim)
  T <- data.frame(time=numeric(nsim),interval = numeric(nsim))
  
  for (i in 1:nsim) {
    tau_cond1 <- -exp_const[i] * lambda[-length(lambda)] * diff_interval
    tau_cond2 <- rev(exp(cumsum(tau_cond1)))
    loc <- J- findInterval(tau[i],tau_cond2)
    lambda_temp <- lambda[0:(loc-1)]
    lambda_minus <- lambda[0:(loc-1)]/lambda[loc]*diff_interval[0:(loc-1)]
    t <- interval[loc] - log(tau[i])/(lambda[loc]*exp_const[i]) - sum(lambda_minus)
    T$time[i] <- t
    T$interval[i] <- loc
  }
  T <- cbind(T,Z,X1,X2,X3,X4,X5,X6,X7,X8,X9,reffect)
  
  ## table(T$inteval)
  return(T)
  
}




### Generate interval censor from given realtime data and visit time, windows

interval_gen <- function(true_time_dt, ### columns, "time","interval","Z,X1 - X9"
                         intervals_1=c(3.29,6,12,24,48,60), ## low risk schedule
                         margin_1 = c(0.1,1,1,1,1,1),## low risk margin
                         multiplier_1 = 1, ## low risk multiplier
                         intervals_2=c(1,2,6,12,24,36,48,60), ## high risk schedule
                         margin_2 = c(0.1,0.1,1.5,1.5,2,2,2,2), ## high risk margin
                         multiplier_2 = 1, ## high risk multiplier
                         realtime_rate = 0.1, ### rate for true event time
                         layer_freq_1 = c(0.25,0.5,0.25), # proportion of each intervals low risk
                         layer_miss_1 = c(0.02,0.2,0.5), # missing rate for each layers low risk
                         layer_freq_2 = c(1), # proportion of each intervals high risk
                         layer_miss_2 = c(0.02), # missing rate for each layers high risk
                         T_multi = 0.001,
                         risk_cov="X9"
){
  
  ### Calculate the visit interval based on visit time and window
  
  ## low risk interval
  pred_interval_L1 <- intervals_1 - margin_1*multiplier_1
  pred_interval_U1 <- intervals_1 + margin_1*multiplier_1
  
  ## high risk interval
  pred_interval_L2 <- intervals_2 - margin_2*multiplier_2
  pred_interval_U2 <- intervals_2 + margin_2*multiplier_2

  

  ### time transformation to meet the right censor proportion
  censor_T_trans <- true_time_dt
  censor_T_trans$time_trans <- censor_T_trans$time/T_multi
  
  
  ### create a new data frame by split the data into different given visit frequency group
  layer <- unlist(mapply(rep, 1:length(layer_freq_1),each=floor(nrow(censor_T_trans)*layer_freq_1)))
  layer <- sample(layer)
  censor_T_trans$layer <- layer
  

  ### For each frequency group, generate interval censor data
  for(i in 1:nrow(censor_T_trans)){
    freq_loc <- censor_T_trans$layer[i]
    
    if(censor_T_trans[i,risk_cov]==0 ){ ## low risk
      gen_censor_interval <- mapply(function(x,y) runif(1,x,y), pred_interval_L1,pred_interval_U1)
      gen_censor_interval <- c(0,gen_censor_interval,Inf)
      
      obs_ind <- 1-rbinom(length(gen_censor_interval)-2,1,layer_miss_1[freq_loc])


      
    }else if (censor_T_trans[i,risk_cov]==1 ){ ## high risk
      gen_censor_interval <- mapply(function(x,y) runif(1,x,y), pred_interval_L2,pred_interval_U2)
      gen_censor_interval <- c(0,gen_censor_interval,Inf)
      
      obs_ind <- 1-rbinom(length(gen_censor_interval)-2,1,layer_miss_2[1])
    }
    
    new_gen_interval <- gen_censor_interval[2:(length(gen_censor_interval)-1)][which(obs_ind==1)]
    new_gen_interval <- c(0,new_gen_interval,Inf)
    obs_loc <- findInterval(censor_T_trans$time_trans[i],new_gen_interval)
    
    censor_T_trans$obs_L_trans[i] <- new_gen_interval[obs_loc]
    censor_T_trans$obs_R_trans[i] <- new_gen_interval[obs_loc+1]
  }

  
  real_sample <- sample(1:nrow(censor_T_trans), nrow(censor_T_trans)*realtime_rate,replace=F )

  censor_T_trans[real_sample,'obs_L_trans'] <- ifelse(censor_T_trans[real_sample,'time_trans'] >= max(intervals_1),max(intervals_1),censor_T_trans[real_sample,'time_trans'])
  censor_T_trans[real_sample,'obs_R_trans'] <- ifelse(censor_T_trans[real_sample,'time_trans'] >= max(intervals_1),Inf,censor_T_trans[real_sample,'time_trans'])
  
  ## Transform back to original generated time scale
  censor_T_trans$obs_L <- censor_T_trans$obs_L_trans*T_multi
  censor_T_trans$obs_R <- censor_T_trans$obs_R_trans*T_multi
  censor_T_trans_final <- censor_T_trans[,-which(names(censor_T_trans) %in% c("layer","time_trans","obs_L_trans","obs_R_trans"))]
  censor_T_trans_final$obs_R[which(is.na(censor_T_trans_final$obs_R))] <- Inf
  
  censor_T_trans_final[which(censor_T_trans_final$obs_L == 0),'censor'] <- 0
  censor_T_trans_final[which(censor_T_trans_final$obs_R == Inf),'censor'] <- 2
  censor_T_trans_final[which(censor_T_trans_final$obs_L !=0 & censor_T_trans_final$obs_R != Inf),'censor'] <- 1
  censor_T_trans_final[which(censor_T_trans_final$obs_L == censor_T_trans_final$obs_R),'censor'] <- 3 ## real event
  
  
  ### Midpoint imputation
  censor_T_trans_final[which(censor_T_trans_final$censor==0 |
                               censor_T_trans_final$censor==1),"obs_L"] <-  censor_T_trans_final[which(censor_T_trans_final$censor==0 |
                                                                                                         censor_T_trans_final$censor==1),"obs_R"] <-
    censor_T_trans_final$obs_R[which(censor_T_trans_final$censor==0 | censor_T_trans_final$censor==1)]
  censor_T_trans_final[which(censor_T_trans_final$censor==0 |
                               censor_T_trans_final$censor==1),"censor"] <- 3
  
  
  
  
  return(censor_T_trans_final)
  
  
}



gen_interval_censor_HC_CT <- function(
  nH = 600,
  nC = 200,
  est = c(1,-1,0.5,-0.5,0.1,-0.1,1,-1,0.5),
  delta = rep(0,9),seed=123,gamma = 1,
  interval = c(0), ## pre-defined time intervals
  lambda0 = exp(1), delta_lambda = exp(0),
  ct_intervals = c(0.92,1.84,5.92,12,24,36,48,60),
  ct_margin = c(0.1,0.1,1.48,1.48,1.97,1.96,1.96,1.96),
  ct_end = 60,
  ct_multi = 1,
  ct_event = 0.1,
  hc_intervals = c(3.29,6,12,24,36,48,60),
  hc_margin = c(0.1,1,1,1,1,1,1),
  hc_end = 60,
  hc_multi = 1,
  hc_event = 0.1,
  risk_cov="X9",
  kappa = c(0,1),
  prob_end = 0.5,
  hc_freq = c(0.25,0.5,0.25),
  hc_miss = c(0.02,0.2,0.5),
  ct_freq = c(1),
  ct_miss = c(0.02),
  adm_censor=F

){
  
  set.seed(seed)
  
  
  
  ### Generate CT data
  
  ### Generate covariate
  beta_ct <- est
  

  
  ### Generate random effect and match in CT and HC
  ### Generate random effect and match in CT and HC
  
  reffect_CT <- rnorm(nC,0,1)
  reffect_HC <- rep(0,nH)
  
  rt <- nH/nC
  
  for(i in 1:nC){
    reffect_HC[(((i-1)*rt)+1) : (((i-1)*rt)+rt) ] <- reffect_CT[i]
  }  
  
  
  
 
  ### generate real time for CT
  
  lambda_ct <- lambda0
  
  
  
 
  CT_T <- piece_exp_gen(nC,beta=c(gamma,beta_ct),seed=seed,reffect = reffect_CT,
                      interval=interval,lambda=lambda_ct)
  
  #T_50 <- quantile(CT_T$time,probs = c(prob_end))
  #T_multi <- T_50/ct_end
  T_multi <- 0.006946491
  
  ### Generate interval censor L, R
  
  
  CT_T_interval <- interval_gen(CT_T,intervals_1=ct_intervals,margin_1=ct_margin,
                                multiplier_1=ct_multi,
                                intervals_2 = ct_intervals,margin_2 = ct_margin,
                                multiplier_2 = ct_multi,
                                realtime_rate = ct_event,layer_freq_1 = ct_freq,layer_miss_1=ct_miss,
                                layer_freq_2 = ct_miss, layer_miss_2 = ct_miss,
                                T_multi = T_multi,risk_cov = risk_cov)
  
  CT_T_interval$obs_R[which(CT_T_interval$obs_R == Inf)] <- NA
  
  CT_T_interval$L_int <- findInterval(CT_T_interval$obs_L,kappa)
  CT_T_interval$R_int <- findInterval(CT_T_interval$obs_R,kappa)
  
  CT_T_interval$cluster <- match(CT_T_interval$reffect,(sort(unique.default(CT_T_interval$reffect))))
  
  
  CT_T_interval$HC <- 0
  
  ### Generate Historical Data HC
  
  ### Generate covariate

  
  lambda_hc <- lambda0 * delta_lambda
  beta_hc <- est + delta
  

  HC_T <- piece_exp_gen(nH,beta=c(0,beta_hc),seed=seed+1,reffect = reffect_HC,
                                  interval=interval,lambda=lambda_hc)
  
  
  ### Generate interval censor L, R
  
  HC_T_interval <-  interval_gen(HC_T,intervals_1=hc_intervals,margin_1=hc_margin,
                                 multiplier_1=hc_multi,
                                 intervals_2 = ct_intervals,margin_2 = ct_margin,
                                 multiplier_2 = ct_multi,
                                 layer_freq_1 = hc_freq,layer_miss_1 =hc_miss,
                                 layer_freq_2 = ct_freq, layer_miss_2 = ct_miss,
                                 realtime_rate = hc_event,
                                 T_multi = T_multi,risk_cov = risk_cov)
  HC_T_interval$Z <- 0
  HC_T_interval$obs_R[which(HC_T_interval$obs_R == Inf)] <- NA
  
  
  
  HC_T_interval$L_int <- findInterval(HC_T_interval$obs_L,kappa)
  HC_T_interval$R_int <- findInterval(HC_T_interval$obs_R,kappa)
  
  HC_T_interval$cluster <- match(HC_T_interval$reffect,(sort(unique.default(HC_T_interval$reffect))))
  
  HC_T_interval$HC <- 1

  
  CT_HC <- rbind(CT_T_interval,HC_T_interval)
  
  
  if(adm_censor == T){
    etime <- CT_HC$time/T_multi
    C1 <- runif(nC+nH,12,60)
    C2 <- rexp(nC+nH,0.05343455)
    C <- apply(cbind(C1,C2), 1, min)
    adm_ind <- which(etime >= C)
    CT_HC[adm_ind,"obs_L"] <- C[adm_ind]*T_multi
    CT_HC[adm_ind,"obs_R"] <- NA
    CT_HC$L_int[adm_ind] <- findInterval(CT_HC$obs_L[adm_ind],kappa)
    CT_HC$R_int[adm_ind] <- findInterval(CT_HC$obs_L[adm_ind],kappa)
    CT_HC[adm_ind,"censor"] <- 2
  }
  
  return(CT_HC)
  
}
