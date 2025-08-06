source("HC CT Data Generation Functions_skip3freq_X_dep.R")

library(R2jags)

ct_interval <- c(1,2,3,6,9,12,18,24,30,36,42,48,54,60)
ct_margin <- c(0.1,0.1,0.1,
               1,1,1,
               1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5)

hc_interval <- c(3.29,6,12,24,36,48,60)
hc_margin <- c(0.2,2,3,3,3,3,3)

interval <- c(0,3,6,9,12)
Monte_out <- list()
M <- 500
for(m in 1:M){
  
  full_dt <- gen_interval_censor_HC_CT(nH=1800,nC=600,est = c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                                       delta=rep(0,9), seed=m*153+3,interval = c(0),gamma=0,
                                       lambda0 = c(exp(2.4)),delta_lambda = c(exp(0)),
                                       ct_intervals = ct_interval,ct_margin = ct_margin,
                                       ct_multi = 1,ct_event=0.1,
                                       hc_intervals=hc_interval,hc_margin = hc_margin,
                                       hc_multi = 1,hc_event = 0.1,
                                       kappa = c(0),prob_end=0.5, hc_freq = c(1),hc_miss = c(0.2),
                                       ct_freq = c(1),ct_miss = c(0.02), adm_censor=T
  )
  
  
  K <- length(interval)
  T_multi <- 0.006946491
  kappa <- T_multi * interval
  kappa[K+1] <- max(full_dt$obs_L,full_dt$obs_R,full_dt$time,na.rm = T) + 0.0001
  # kappa <- c(0,quantile(full_dt$time,prob=c(0.5)))
  ##  kappa <- c(0)
  ##  idx.kappasE <- which(full_dt$censor==3)
  ##  kappa[K+1] <- max(c(full_dt$obs_L,na.omit(full_dt$obs_R),full_dt$time[idx.kappasE])) + 0.0001
  
  
  full_dt$L_int <- findInterval(full_dt$obs_L,kappa)
  full_dt$R_int <- findInterval(full_dt$obs_R,kappa)
  
  out.H <- full_dt[full_dt$HC==1,] ## Historical
  out.C <- full_dt[full_dt$HC==0,] ## CT
  
  
  
  ### HC data
  out.H <- out.H[order(out.H$censor),]
  idxL = max(which(out.H$censor==0))
  idxI = c(min(which(out.H$censor==1)),max(which(out.H$censor==1)))
  idxR = max(which(out.H$censor==2))
  idxE = max(which(out.H$censor==3))
  
  X <- out.H[,c(paste("X",1:9,sep=""))]
  p = ncol(X)
  L = out.H$obs_L
  R = out.H$obs_R
  lower.int.obs = out.H$L_int
  upper.int.obs = out.H$R_int
  
  ones <- rep(1,nrow(out.H))
  cluster <- out.H$cluster
  n <- length(cluster)
  
  ### CT data
  out.C <- out.C[order(out.C$censor),]
  idxL.C = max(which(out.C$censor==0))
  idxI.C = c(min(which(out.C$censor==1)),max(which(out.C$censor==1)))
  idxR.C = max(which(out.C$censor==2))
  idxE.C = max(which(out.C$censor==3))
  
  X.C <- out.C[,c(paste("X",1:9,sep=""))]
  p = ncol(X.C)
  trt.C <- out.C$Z
  L.C = out.C$obs_L
  R.C = out.C$obs_R
  
  lower.int.obs.C = out.C$L_int
  upper.int.obs.C = out.C$R_int
  ones.C <- rep(1,nrow(out.C))
  cluster.C <- out.C$cluster
  n.C <- length(cluster.C)
  
  
  idxcluster <- rep(0,n)
  for(i in 1 : n){
    idxcluster[i] <- which(cluster[i]==cluster.C)
  }
  nint = length(kappa)-1
  
  
  model.data <- list("idxL","idxI","idxR","idxE","lower.int.obs","upper.int.obs","X","p","nint","kappa","R","L","ones",
                     "idxL.C","idxI.C","idxR.C","idxE.C","lower.int.obs.C","upper.int.obs.C","X.C","R.C","L.C","ones.C","idxcluster",
                     "n","n.C","trt.C")

  model.parameters <- model.parameters <- c("beta","alpha","nu","betaC","alphaC","gammaC","zz","cc","aaa1")  #Clustering method
  model.fit.cluster <- jags(model.file= "exact_likelihood_random_effect_HC_CT.txt", data=model.data,
  parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)
  Monte_out[[m]] <- model.fit.cluster
  saveRDS(Monte_out, file = "fully_ex2.4(3).rds")
}