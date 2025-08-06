setwd("/scratch/g/skim/CT_HD_Bayesian_ones(1228)/Midpoint_X_depend_5_interval/CP")



tau_dist_fun <- function(indata){
  tau_quant_dt <- array(0,dim = c(3,14,length(indata)))
  
  for(m in 1:length(indata)){
    tau_dist  <- matrix(nrow=nrow(indata[[m]]$BUGSoutput$sims.list$cc),
                        ncol=ncol(indata[[m]]$BUGSoutput$sims.list$cc))
    for(i in 1:ncol(tau_dist)){
      if(i<=10){
        tau_dist[,i] <- ifelse(indata[[m]]$BUGSoutput$sims.list$cc[,i] == 1,
                               indata[[m]]$BUGSoutput$sims.list$zz[,1], 
                               indata[[m]]$BUGSoutput$sims.list$zz[,2])
      }else if(i>10){
        tau_dist[,i] <- ifelse(indata[[m]]$BUGSoutput$sims.list$cc[,i] == 1,
                               indata[[m]]$BUGSoutput$sims.list$aaa1, 
                               indata[[m]]$BUGSoutput$sims.list$zz[,2])
      }

    }
    tau_quant_dt[,,m] <- apply(tau_dist, 2, quantile,probs=c(0.25,0.5,0.75))
    
  }
  
  tau_25 <- apply(t(tau_quant_dt[1,,]),2,mean)
  tau_50 <- apply(t(tau_quant_dt[2,,]),2,mean)
  tau_75 <- apply(t(tau_quant_dt[3,,]),2,mean)
  tau_final <- t(rbind(tau_25,tau_50,tau_75))
  return(tau_final)
}



## fully exchangeable
fullyex <- readRDS("fully_ex1.rds")
fullyex2 <- readRDS("fully_ex1(2).rds")
fullyex <- c(fullyex,fullyex2)
fully_tau <- tau_dist_fun(fullyex)
write.csv(fully_tau,file="fully1_tau.csv")


fullyex <- readRDS("fully_ex2.4.rds")
fullyex2 <- readRDS("fully_ex2.4(2).rds")
fullyex <- c(fullyex,fullyex2)
fully_tau <- tau_dist_fun(fullyex)
write.csv(fully_tau,file="fully2.4_tau.csv")


# nonfully 0.2
non02 <- readRDS("nonfully1_02.rds")
non022 <- readRDS("nonfully1_02(2).rds")
non02 <- c(non02,non022)

non02_tau <- tau_dist_fun(non02)
write.csv(non02_tau,file="non1_02_tau.csv")


non02 <- readRDS("nonfully2.4_02.rds")
non022 <- readRDS("nonfully2.4_02(2).rds")
non02 <- c(non02,non022)
non02_tau <- tau_dist_fun(non02)
write.csv(non02_tau,file="non2.4_02_tau.csv")


# nonfully 0.4
non04 <- readRDS("nonfully1_04.rds")
non042 <- readRDS("nonfully1_04(2).rds")
non04 <- c(non04,non042)
non04_tau <- tau_dist_fun(non04)
write.csv(non04_tau,file="non1_04_tau.csv")


non04 <- readRDS("nonfully2.4_04.rds")
non042 <- readRDS("nonfully2.4_04(2).rds")
non04 <- c(non04,non042)
non04_tau <- tau_dist_fun(non04)
write.csv(non04_tau,file="non2.4_04_tau.csv")



# nonfully 0.8
non08 <- readRDS("nonfully1_08.rds")
non082 <- readRDS("nonfully1_08(2).rds")
non08 <- c(non08,non082)
non08_tau <- tau_dist_fun(non08)
write.csv(non08_tau,file="non1_08_tau.csv")


non08 <- readRDS("nonfully2.4_08.rds")
non082 <- readRDS("nonfully2.4_08(2).rds")
non08 <- c(non08,non082)
non08_tau <- tau_dist_fun(non08)
write.csv(non08_tau,file="non2.4_08_tau.csv")



# nonfully 1.2
non1.2 <- readRDS("nonfully1_1.2.rds")
non1.22 <- readRDS("nonfully1_1.2(2).rds")

non1.2 <- c(non1.2,non1.22)

non1.2_tau <- tau_dist_fun(non1.2)
write.csv(non1.2_tau,file="non1_1.2_tau.csv")



non1.2 <- readRDS("nonfully2.4_1.2.rds")
non1.22 <- readRDS("nonfully2.4_1.2(2).rds")
non1.2 <- c(non1.2,non1.22)
non1.2_tau <- tau_dist_fun(non1.2)
write.csv(non1.2_tau,file="non2.4_1.2_tau.csv")


# nonfully 1.6
non1.6 <- readRDS("nonfully1_1.6.rds")
non1.62 <- readRDS("nonfully1_1.6(2).rds")
non1.6 <- c(non1.6,non1.62)
non1.6_tau <- tau_dist_fun(non1.6)
write.csv(non1.6_tau,file="non1_1.6_tau.csv")


non1.6 <- readRDS("nonfully2.4_1.6.rds")
non1.62 <- readRDS("nonfully2.4_1.6(2).rds")
non1.6 <- c(non1.6,non1.62)
non1.6_tau <- tau_dist_fun(non1.6)
write.csv(non1.6_tau,file="non2.4_1.6_tau.csv")


# nonfully 2.0
non2.0 <- readRDS("nonfully1_2.0.rds")
non2.02 <- readRDS("nonfully1_2.0(2).rds")
non2.0 <- c(non2.0,non2.02)
non2.0_tau <- tau_dist_fun(non2.0)
write.csv(non2.0_tau,file="non1_2.0_tau.csv")


non2.0 <- readRDS("nonfully2.4_2.0.rds")
non2.02 <- readRDS("nonfully2.4_2.0(2).rds")
non2.0 <- c(non2.0,non2.02)
non2.0_tau <- tau_dist_fun(non2.0)
write.csv(non2.0_tau,file="non2.4_2.0_tau.csv")

