library(Rcpp)

library(survival)
library(MatchIt)
library(R2jags)

raw_dt <- read.csv("BMT_data.csv")
raw_dt <- raw_dt[which(raw_dt$main != 5),]

raw_dt$source <- ifelse(raw_dt$main %in% c(1,3), "CT","HD")
raw_dt$treatment <- ifelse(raw_dt$main %in% c(1,2), 0,1  ) # 1: Haplo, 0: dUCB



raw_dt$disease <- ifelse(raw_dt$disease2gp %in% c(1,2,3), 0, 1 ) # 0:Leukemia, 1:Lymphoma
raw_dt$raceb <- ifelse(raw_dt$race2gp==2, 1 ,0)
raw_dt$raceo <- ifelse(raw_dt$race2gp==3, 1 ,0)
raw_dt$eth <- ifelse(raw_dt$ethnicit == 1, 0, 1)
raw_dt$red_dis_int <- ifelse(raw_dt$drigp == 2, 1,0 )
raw_dt$red_dis_high <- ifelse(raw_dt$drigp %in% c(3,4), 1, 0 )
# raw_dt$red_dis_veryhigh <- ifelse(raw_dt$drigp==4, 1, 0 )
raw_dt$dis_trans_CR <- ifelse(raw_dt$alstatprgpnew_num %in% c(1:5),1,0) ## CR
raw_dt$dis_trans_PR <- ifelse(raw_dt$alstatprgpnew_num == 6,1,0)  ## PR
raw_dt$karno <- ifelse(raw_dt$karnofcat==1,0,1)
raw_dt$hct_re <- ifelse(raw_dt$hct==1,0,1)
raw_dt$gender <- ifelse(raw_dt$sex==1,0,1) # 1=male
raw_dt$CT <- ifelse(raw_dt$source=="CT",1,0)

ana_data <- raw_dt[,c(15,13,14,10,11,12,21,3,32,25,30,26,27,8,
                      31,22,23,24,33)]

ana_data <- as.matrix(ana_data)
ana_data[which(ana_data==99)] <- NA
ana_data <- na.omit(ana_data)
ana_data <- as.data.frame(ana_data)

ana_data <- ana_data[which(ana_data$disease==0),]

set.seed(1)
#match_form <- as.formula("CT~age+gender+eth+karno+red_dis_int+red_dis_high+red_dis_high+
#                         cmv+hct_re+disease+raceb+raceo")

match_form <- as.formula("CT~age+gender+eth+karno+red_dis_int+red_dis_high+red_dis_high+
                         cmv+hct_re+raceb+raceo")

m.out <- matchit(match_form,data=ana_data,caliper = 0.2,ratio=2)
matched.data <- match.data(m.out,drop.unmatched = F)

K = 10 #Number of intervals for alpha
kappa=c(0,quantile(matched.data$intxrel[which(matched.data$dfs==1 & matched.data$CT==1)],
                   seq(0,1,by=1/K))[-c(1,K+1)],
        max(matched.data$intxrel)+0.0001) 

hc_interval <- c(0,3.29,6,12,24,36,48,60,Inf)
ct_interval <- c(0,7,14,21,28,35,42,49,56,180,365,730, Inf)/30.4

matched.data$censor <- NA
matched.data$obs_L <- NA
matched.data$obs_R <- NA
matched.data$L_int <- NA
matched.data$R_int <- NA

for(i in 1:nrow(matched.data)){
  if(matched.data$dead[i]==1){
    matched.data$obs_L[i] <- matched.data$intxrel[i]
    matched.data$obs_R[i] <- matched.data$intxrel[i]
  }else if (matched.data$dfs[i]==1){
    if(matched.data$CT[i]==1){
      temp_ind <- findInterval(matched.data$intxrel[i],ct_interval)
      matched.data$obs_L[i] <- ct_interval[temp_ind]
      matched.data$obs_R[i] <- matched.data$intxrel[i]
    }else{
      temp_ind <- findInterval(matched.data$intxrel[i],hc_interval)
      matched.data$obs_L[i] <- hc_interval[temp_ind]
      matched.data$obs_R[i] <- matched.data$intxrel[i]
    }
    
    
  }else if(matched.data$dfs[i]==0){
    matched.data$obs_L[i] <- matched.data$intxrel[i]
    matched.data$obs_R[i] <- Inf
  }
  
  if(matched.data$obs_L[i] == 0){
    matched.data$censor[i] <- 0
  }else if(matched.data$obs_R[i] == Inf){
    matched.data$censor[i] <- 2
  }else if (matched.data$obs_L[i] != 0 & matched.data$obs_R[i] != Inf ){
    matched.data$censor[i] <- 1
  }
  
  if(matched.data$dead[i]==1){
    matched.data$censor[i] <- 3
  }
  
  
  
}


matched.data$L_int <- findInterval(matched.data$obs_L,kappa)
matched.data$R_int <- findInterval(matched.data$obs_R,kappa)

matched.data <- matched.data[-which(is.na(matched.data$subclass) & matched.data$CT==0 ),]

matched.data$subclass <- as.numeric(levels(matched.data$subclass))[matched.data$subclass]
matched.data$subclass[which(is.na(matched.data$subclass))] <- (max(as.numeric(matched.data$subclass),na.rm=T)+1) : (sum(is.na(matched.data$subclass)) + max(as.numeric(matched.data$subclass),na.rm=T))
matched.data <- matched.data[order(matched.data$subclass),]



## interval censored data
data_all <- matched.data[order(matched.data$censor),]
out.H=data_all[data_all$CT==0,];out.C=data_all[data_all$CT==1,]



## HD data
out.H <- out.H[order(out.H$censor),]
idxL = max(which(out.H$censor==0))
idxI = c(min(which(out.H$censor==1)),max(which(out.H$censor==1)))
idxR = max(which(out.H$censor==2))
idxE = max(which(out.H$censor==3))

X <- out.H[,c(7:15,17:18)]
p = ncol(X)
L = out.H$obs_L
R = out.H$obs_R
lower.int.obs = out.H$L_int
upper.int.obs = out.H$R_int

ones <- rep(1,nrow(out.H))
cluster <- out.H$subclass
n <- length(cluster)


### CT data
out.C <- out.C[order(out.C$censor),]
idxL.C = max(which(out.C$censor==0))
idxI.C = c(min(which(out.C$censor==1)),max(which(out.C$censor==1)))
idxR.C = max(which(out.C$censor==2))
idxE.C = max(which(out.C$censor==3))

X.C <- out.C[,c(7:15,17:18)]
p = ncol(X.C)
##trt.C <- out.C$Trt
L.C = out.C$obs_L
R.C = out.C$obs_R


lower.int.obs.C = out.C$L_int
upper.int.obs.C = out.C$R_int
ones.C <- rep(1,nrow(out.C))
cluster.C <- out.C$subclass
n.C <- length(cluster.C)


idxcluster <- rep(0,n)
for(i in 1 : n){
  idxcluster[i] <- which(cluster[i]==cluster.C)
}
nint = length(kappa)-1









#    JAGS Code   For Murray's method #
#model.data <- list("idxR","idxE","lower.int.obs","X","p","nint","kappa","L","ones",
#                   "idxR.C","idxE.C","lower.int.obs.C","X.C","L.C","ones.C","idxcluster",
#                   "n","n.C")


model.data <- list("idxL","idxI","idxR","idxE","lower.int.obs","upper.int.obs","X","p","nint","kappa","R","L","ones",
                   "idxL.C","idxI.C","idxR.C","idxE.C","lower.int.obs.C","upper.int.obs.C","X.C","R.C","L.C","ones.C","idxcluster",
                   "n","n.C")

model.parameters <- c("beta","alpha","nu","betaC","alphaC","t","mmm","mmm1") # Murray's Method


model.fit.cluster <- jags(model.file= "exact_likelihood_random_effect_HC_CT_murray.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)


survival_murrray_parm <- model.fit.cluster$BUGSoutput$summary[33:44,c(1,3,7)]
rownames(survival_murrray_parm) <- colnames(X.C)
survival_murrray_parm <- cbind(survival_murrray_parm,
                               paste(round(survival_murrray_parm[,1],2)," (",
                                     round(survival_murrray_parm[,2],2),",",
                                     round(survival_murrray_parm[,3],2),")",sep="")
)

write.csv(survival_murrray_parm,file="survival_murray_parm.csv")



tau_quant <- matrix(nrow=nrow(model.fit.cluster$BUGSoutput$sims.list$t),
                    ncol=ncol(model.fit.cluster$BUGSoutput$sims.list$t))

for(p in 1:ncol(model.fit.cluster$BUGSoutput$sims.list$t)){
  if (p <= 12 ){
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$t[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$mmm,200)
  }else if (p == 13){
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$t[,p]==1,
                            0.0001,200)
  }else{
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$t[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$mmm1,200)
  }
  
}


#tau_quant <- read.csv("survival_murray_tau.csv")
apply(tau_quant, 2, median)
write.csv(tau_quant,file="survival_murray_tau.csv")



model.parameters <- model.parameters <- c("beta","alpha","nu","betaC","alphaC","zz","cc","aaa1")  #Clustering method
model.fit.cluster <- jags(model.file= "exact_likelihood_random_effect_HC_CT.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_clustering_parm <- model.fit.cluster$BUGSoutput$summary[33:43,c(1,3,7)]
rownames(survival_clustering_parm) <- colnames(X.C)

survival_clustering_parm <- cbind(survival_clustering_parm,
                                  paste(round(survival_clustering_parm[,1],2)," (",
                                        round(survival_clustering_parm[,2],2),",",
                                        round(survival_clustering_parm[,3],2),")",sep="")
)


write.csv(survival_clustering_parm,file="survival_clustering_parm.csv")


tau_quant <- matrix(nrow=nrow(model.fit.cluster$BUGSoutput$sims.list$cc),
                    ncol=ncol(model.fit.cluster$BUGSoutput$sims.list$cc))

for(p in 1:ncol(model.fit.cluster$BUGSoutput$sims.list$cc)){
  if (p <= 13 ){
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$cc[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$zz[,1],
                            model.fit.cluster$BUGSoutput$sims.list$zz[,2])
  }else{
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$cc[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$aaa1,
                            model.fit.cluster$BUGSoutput$sims.list$zz[,2])
  }
  
}

apply(tau_quant, 2, median)
tau_quant <- read.csv("survival_clustering_tau.csv")
write.csv(tau_quant,file="survival_clustering_tau.csv")


## Meta analysis
model.parameters <- model.parameters <- c("beta","alpha","nu","betaC","alphaC","tau_sb","tau_sa","mu_sb","mu_sa") 
model.fit.cluster <- jags(model.file= "exact_likelihood_random_effect_shrinksame.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_meta_parm <- model.fit.cluster$BUGSoutput$summary[33:44,c(1,3,7)]
rownames(survival_meta_parm) <- colnames(X.C)
survival_meta_parm <- cbind(survival_meta_parm,
                            paste(round(survival_meta_parm[,1],2)," (",
                                  round(survival_meta_parm[,2],2),",",
                                  round(survival_meta_parm[,3],2),")",sep="")
)

write.csv(survival_meta_parm,file="survival_meta_parm.csv")


## Fit CT all separately
model.data <- list("idxL","idxI","idxR","idxE","lower.int.obs","upper.int.obs","X","p","nint","kappa","R","L","ones",
                   "idxL.C","idxI.C","idxR.C","idxE.C","lower.int.obs.C","upper.int.obs.C","X.C","R.C","L.C","ones.C","idxcluster",
                   "n","n.C")

model.parameters <- c("beta","alpha","nu","betaC","alphaC","t","mmm","mmm1") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_CT.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_sep <- model.fit.sep$BUGSoutput$summary[11:22,c(1,3,7)]
rownames(survival_all_sep) <- colnames(X.C)
survival_all_sep <- cbind(survival_all_sep,
                          paste(round(survival_all_sep[,1],2)," (",
                                round(survival_all_sep[,2],2),",",
                                round(survival_all_sep[,3],2),")",sep="")
)

write.csv(survival_all_sep,file="survival_sep_CT_parm.csv")


model.data <- list("idxL","idxI","idxR","idxE","lower.int.obs","upper.int.obs","X","p","nint","kappa","R","L","ones",
                   "idxL.C","idxI.C","idxR.C","idxE.C","lower.int.obs.C","upper.int.obs.C","X.C","R.C","L.C","ones.C","idxcluster",
                   "n","n.C")

model.parameters <- c("beta","alpha","nu","betaC","alphaC","t","mmm","mmm1") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_HD.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_sep <- model.fit.sep$BUGSoutput$summary[11:22,c(1,3,7)]
rownames(survival_all_sep) <- colnames(X)
survival_all_sep <- cbind(survival_all_sep,
                          paste(round(survival_all_sep[,1],2)," (",
                                round(survival_all_sep[,2],2),",",
                                round(survival_all_sep[,3],2),")",sep="")
)

write.csv(survival_all_sep,file="survival_sep_HD_parm.csv")





## Midpoint

data_all <- matched.data[order(matched.data$censor),]

data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )] <- 
  (data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] +
     data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )])/2

data_all$censor[which(data_all$censor==0 | data_all$censor==1 )] <- 3

data_all$L_int <- findInterval(data_all$obs_L,kappa)
data_all$R_int <- findInterval(data_all$obs_R,kappa)



out.H=data_all[data_all$CT==0,];out.C=data_all[data_all$CT==1,]



## HD data
out.H <- out.H[order(out.H$censor),]
idxL = 0
idxI = c(0,0)
idxR = max(which(out.H$censor==2))
idxE = max(which(out.H$censor==3))

X <- out.H[,c(7:15,17:18)]
p = ncol(X)
L = out.H$obs_L
R = out.H$obs_R
lower.int.obs = out.H$L_int
upper.int.obs = out.H$R_int

ones <- rep(1,nrow(out.H))
cluster <- out.H$subclass
n <- length(cluster)


### CT data
out.C <- out.C[order(out.C$censor),]
idxL.C = 0
idxI.C = c(0,0)
idxR.C = max(which(out.C$censor==2))
idxE.C = max(which(out.C$censor==3))

X.C <- out.C[,c(7:15,17:18)]
p = ncol(X.C)
##trt.C <- out.C$Trt
L.C = out.C$obs_L
R.C = out.C$obs_R


lower.int.obs.C = out.C$L_int
upper.int.obs.C = out.C$R_int
ones.C <- rep(1,nrow(out.C))
cluster.C <- out.C$subclass
n.C <- length(cluster.C)


idxcluster <- rep(0,n)
for(i in 1 : n){
  idxcluster[i] <- which(cluster[i]==cluster.C)
}
nint = length(kappa)-1

model.data <- list("idxL","idxI","idxR","idxE","lower.int.obs","upper.int.obs","X","p","nint","kappa","R","L","ones",
                   "idxL.C","idxI.C","idxR.C","idxE.C","lower.int.obs.C","upper.int.obs.C","X.C","R.C","L.C","ones.C","idxcluster",
                   "n","n.C")

model.parameters <- c("beta","alpha","nu","betaC","alphaC","t","mmm","mmm1") # Murray's Method


model.fit.cluster <- jags(model.file= "exact_likelihood_random_effect_HC_CT_murray_imp.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_murrray_parm <- model.fit.cluster$BUGSoutput$summary[32:42,c(1,3,7)]
rownames(survival_murrray_parm) <- colnames(X.C)
survival_murrray_parm <- cbind(survival_murrray_parm,
                               paste(round(survival_murrray_parm[,1],2)," (",
                                     round(survival_murrray_parm[,2],2),",",
                                     round(survival_murrray_parm[,3],2),")",sep="")
)

write.csv(survival_murrray_parm,file="survival_murray_parm_midpoint.csv")



tau_quant <- matrix(nrow=nrow(model.fit.cluster$BUGSoutput$sims.list$t),
                    ncol=ncol(model.fit.cluster$BUGSoutput$sims.list$t))

for(p in 1:ncol(model.fit.cluster$BUGSoutput$sims.list$t)){
  if (p <= 12 ){
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$t[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$mmm,200)
  }else if (p == 13){
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$t[,p]==1,
                            0.0001,200)
  }else{
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$t[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$mmm1,200)
  }
  
}


#tau_quant <- read.csv("survival_murray_tau.csv")
apply(tau_quant, 2, median)
write.csv(tau_quant,file="survival_murray_tau.csv")



model.parameters <- model.parameters <- c("beta","alpha","nu","betaC","alphaC","zz","cc","aaa1")  #Clustering method
model.fit.cluster <- jags(model.file= "exact_likelihood_random_effect_HC_CT_imp.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_clustering_parm <- model.fit.cluster$BUGSoutput$summary[33:43,c(1,3,7)]
rownames(survival_clustering_parm) <- colnames(X.C)

survival_clustering_parm <- cbind(survival_clustering_parm,
                                  paste(round(survival_clustering_parm[,1],2)," (",
                                        round(survival_clustering_parm[,2],2),",",
                                        round(survival_clustering_parm[,3],2),")",sep="")
)


write.csv(survival_clustering_parm,file="survival_clustering_parm_imp.csv")


tau_quant <- matrix(nrow=nrow(model.fit.cluster$BUGSoutput$sims.list$cc),
                    ncol=ncol(model.fit.cluster$BUGSoutput$sims.list$cc))

for(p in 1:ncol(model.fit.cluster$BUGSoutput$sims.list$cc)){
  if (p <= 13 ){
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$cc[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$zz[,1],
                            model.fit.cluster$BUGSoutput$sims.list$zz[,2])
  }else{
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$cc[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$aaa1,
                            model.fit.cluster$BUGSoutput$sims.list$zz[,2])
  }
  
}

apply(tau_quant, 2, median)
tau_quant <- read.csv("survival_clustering_tau.csv")
write.csv(tau_quant,file="survival_clustering_tau.csv")


## Meta analysis
model.parameters <- model.parameters <- c("beta","alpha","nu","betaC","alphaC","tau_sb","tau_sa","mu_sb","mu_sa") 
model.fit.cluster <- jags(model.file= "exact_likelihood_random_effect_shrinksame_imp.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_meta_parm <- model.fit.cluster$BUGSoutput$summary[32:42,c(1,3,7)]
rownames(survival_meta_parm) <- colnames(X.C)
survival_meta_parm <- cbind(survival_meta_parm,
                            paste(round(survival_meta_parm[,1],2)," (",
                                  round(survival_meta_parm[,2],2),",",
                                  round(survival_meta_parm[,3],2),")",sep="")
)

write.csv(survival_meta_parm,file="survival_meta_parm_midpoint.csv")



## Fit CT all separately
model.data <- list("idxL","idxI","idxR","idxE","lower.int.obs","upper.int.obs","X","p","nint","kappa","R","L","ones",
                   "idxL.C","idxI.C","idxR.C","idxE.C","lower.int.obs.C","upper.int.obs.C","X.C","R.C","L.C","ones.C","idxcluster",
                   "n","n.C")

model.parameters <- c("beta","alpha","nu","betaC","alphaC","t","mmm","mmm1") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_CT_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_sep <- model.fit.sep$BUGSoutput$summary[11:21,c(1,3,7)]
rownames(survival_all_sep) <- colnames(X.C)
survival_all_sep <- cbind(survival_all_sep,
                          paste(round(survival_all_sep[,1],2)," (",
                                round(survival_all_sep[,2],2),",",
                                round(survival_all_sep[,3],2),")",sep="")
)

write.csv(survival_all_sep,file="survival_sep_CT_parm_midpoint.csv")


model.data <- list("idxL","idxI","idxR","idxE","lower.int.obs","upper.int.obs","X","p","nint","kappa","R","L","ones",
                   "idxL.C","idxI.C","idxR.C","idxE.C","lower.int.obs.C","upper.int.obs.C","X.C","R.C","L.C","ones.C","idxcluster",
                   "n","n.C")

model.parameters <- c("beta","alpha","nu","betaC","alphaC","t","mmm","mmm1") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_HD_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_sep <- model.fit.sep$BUGSoutput$summary[11:21,c(1,3,7)]
rownames(survival_all_sep) <- colnames(X)
survival_all_sep <- cbind(survival_all_sep,
                          paste(round(survival_all_sep[,1],2)," (",
                                round(survival_all_sep[,2],2),",",
                                round(survival_all_sep[,3],2),")",sep="")
)

write.csv(survival_all_sep,file="survival_sep_HD_parm_midpoint.csv")






# right imputation

data_all <- matched.data[order(matched.data$censor),]

data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )] <- 
     data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )]

data_all$censor[which(data_all$censor==0 | data_all$censor==1 )] <- 3

data_all$L_int <- findInterval(data_all$obs_L,kappa)
data_all$R_int <- findInterval(data_all$obs_R,kappa)



out.H=data_all[data_all$CT==0,];out.C=data_all[data_all$CT==1,]



## HD data
out.H <- out.H[order(out.H$censor),]
idxL = 0
idxI = c(0,0)
idxR = max(which(out.H$censor==2))
idxE = max(which(out.H$censor==3))

X <- out.H[,c(7:15,17:18)]
p = ncol(X)
L = out.H$obs_L
R = out.H$obs_R
lower.int.obs = out.H$L_int
upper.int.obs = out.H$R_int

ones <- rep(1,nrow(out.H))
cluster <- out.H$subclass
n <- length(cluster)


### CT data
out.C <- out.C[order(out.C$censor),]
idxL.C = 0
idxI.C = c(0,0)
idxR.C = max(which(out.C$censor==2))
idxE.C = max(which(out.C$censor==3))

X.C <- out.C[,c(7:15,17:18)]
p = ncol(X.C)
##trt.C <- out.C$Trt
L.C = out.C$obs_L
R.C = out.C$obs_R


lower.int.obs.C = out.C$L_int
upper.int.obs.C = out.C$R_int
ones.C <- rep(1,nrow(out.C))
cluster.C <- out.C$subclass
n.C <- length(cluster.C)


idxcluster <- rep(0,n)
for(i in 1 : n){
  idxcluster[i] <- which(cluster[i]==cluster.C)
}
nint = length(kappa)-1

model.data <- list("idxL","idxI","idxR","idxE","lower.int.obs","upper.int.obs","X","p","nint","kappa","R","L","ones",
                   "idxL.C","idxI.C","idxR.C","idxE.C","lower.int.obs.C","upper.int.obs.C","X.C","R.C","L.C","ones.C","idxcluster",
                   "n","n.C")

model.parameters <- c("beta","alpha","nu","betaC","alphaC","t","mmm","mmm1") # Murray's Method


model.fit.cluster <- jags(model.file= "exact_likelihood_random_effect_HC_CT_murray_imp.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_murrray_parm <- model.fit.cluster$BUGSoutput$summary[32:42,c(1,3,7)]
rownames(survival_murrray_parm) <- colnames(X.C)
survival_murrray_parm <- cbind(survival_murrray_parm,
                               paste(round(survival_murrray_parm[,1],2)," (",
                                     round(survival_murrray_parm[,2],2),",",
                                     round(survival_murrray_parm[,3],2),")",sep="")
)

write.csv(survival_murrray_parm,file="survival_murray_parm_right.csv")



tau_quant <- matrix(nrow=nrow(model.fit.cluster$BUGSoutput$sims.list$t),
                    ncol=ncol(model.fit.cluster$BUGSoutput$sims.list$t))

for(p in 1:ncol(model.fit.cluster$BUGSoutput$sims.list$t)){
  if (p <= 12 ){
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$t[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$mmm,200)
  }else if (p == 13){
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$t[,p]==1,
                            0.0001,200)
  }else{
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$t[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$mmm1,200)
  }
  
}


#tau_quant <- read.csv("survival_murray_tau.csv")
apply(tau_quant, 2, median)
write.csv(tau_quant,file="survival_murray_tau.csv")



model.parameters <- model.parameters <- c("beta","alpha","nu","betaC","alphaC","zz","cc","aaa1")  #Clustering method
model.fit.cluster <- jags(model.file= "exact_likelihood_random_effect_HC_CT_imp.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_clustering_parm <- model.fit.cluster$BUGSoutput$summary[33:43,c(1,3,7)]
rownames(survival_clustering_parm) <- colnames(X.C)

survival_clustering_parm <- cbind(survival_clustering_parm,
                                  paste(round(survival_clustering_parm[,1],2)," (",
                                        round(survival_clustering_parm[,2],2),",",
                                        round(survival_clustering_parm[,3],2),")",sep="")
)


write.csv(survival_clustering_parm,file="survival_clustering_parm_right.csv")


tau_quant <- matrix(nrow=nrow(model.fit.cluster$BUGSoutput$sims.list$cc),
                    ncol=ncol(model.fit.cluster$BUGSoutput$sims.list$cc))

for(p in 1:ncol(model.fit.cluster$BUGSoutput$sims.list$cc)){
  if (p <= 13 ){
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$cc[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$zz[,1],
                            model.fit.cluster$BUGSoutput$sims.list$zz[,2])
  }else{
    tau_quant[,p] <- ifelse(model.fit.cluster$BUGSoutput$sims.list$cc[,p]==1,
                            model.fit.cluster$BUGSoutput$sims.list$aaa1,
                            model.fit.cluster$BUGSoutput$sims.list$zz[,2])
  }
  
}

apply(tau_quant, 2, median)
tau_quant <- read.csv("survival_clustering_tau.csv")
write.csv(tau_quant,file="survival_clustering_tau.csv")


## Meta analysis
model.parameters <- model.parameters <- c("beta","alpha","nu","betaC","alphaC","tau_sb","tau_sa","mu_sb","mu_sa") 
model.fit.cluster <- jags(model.file= "exact_likelihood_random_effect_shrinksame_imp.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_meta_parm <- model.fit.cluster$BUGSoutput$summary[32:42,c(1,3,7)]
rownames(survival_meta_parm) <- colnames(X.C)
survival_meta_parm <- cbind(survival_meta_parm,
                            paste(round(survival_meta_parm[,1],2)," (",
                                  round(survival_meta_parm[,2],2),",",
                                  round(survival_meta_parm[,3],2),")",sep="")
)

write.csv(survival_meta_parm,file="survival_meta_parm_right.csv")

## Fit CT all separately
model.data <- list("idxL","idxI","idxR","idxE","lower.int.obs","upper.int.obs","X","p","nint","kappa","R","L","ones",
                   "idxL.C","idxI.C","idxR.C","idxE.C","lower.int.obs.C","upper.int.obs.C","X.C","R.C","L.C","ones.C","idxcluster",
                   "n","n.C")

model.parameters <- c("beta","alpha","nu","betaC","alphaC","t","mmm","mmm1") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_CT_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_sep <- model.fit.sep$BUGSoutput$summary[11:21,c(1,3,7)]
rownames(survival_all_sep) <- colnames(X.C)
survival_all_sep <- cbind(survival_all_sep,
                          paste(round(survival_all_sep[,1],2)," (",
                                round(survival_all_sep[,2],2),",",
                                round(survival_all_sep[,3],2),")",sep="")
)

write.csv(survival_all_sep,file="survival_sep_CT_parm_right.csv")


model.data <- list("idxL","idxI","idxR","idxE","lower.int.obs","upper.int.obs","X","p","nint","kappa","R","L","ones",
                   "idxL.C","idxI.C","idxR.C","idxE.C","lower.int.obs.C","upper.int.obs.C","X.C","R.C","L.C","ones.C","idxcluster",
                   "n","n.C")

model.parameters <- c("beta","alpha","nu","betaC","alphaC","t","mmm","mmm1") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_HD_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_sep <- model.fit.sep$BUGSoutput$summary[11:21,c(1,3,7)]
rownames(survival_all_sep) <- colnames(X)
survival_all_sep <- cbind(survival_all_sep,
                          paste(round(survival_all_sep[,1],2)," (",
                                round(survival_all_sep[,2],2),",",
                                round(survival_all_sep[,3],2),")",sep="")
)

write.csv(survival_all_sep,file="survival_sep_HD_parm_right.csv")

ana_data[which(ana_data$CT==1),] %>% group_by(treatment) %>% summarise_at(vars(age),list(age_mean=mean,age_se=sd))

ana_data[which(ana_data$CT==1),] %>% group_by(treatment,gender) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))

ana_data[which(ana_data$CT==1),] %>% group_by(treatment,eth) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))

ana_data[which(ana_data$CT==1),] %>% group_by(treatment,karno) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))

ana_data$rdri_re <- ifelse(ana_data$red_dis_int==1,1,
                           ifelse(ana_data$red_dis_high==1,2,0))


ana_data[which(ana_data$CT==1),] %>% group_by(treatment,rdri_re) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))

ana_data[which(ana_data$CT==1),] %>% group_by(treatment,cmv) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))


ana_data[which(ana_data$CT==1),] %>% group_by(treatment,hct_re) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))

ana_data$race_re <-  ifelse(ana_data$raceb==1,1,
                            ifelse(ana_data$raceo==1,2,0))

ana_data[which(ana_data$CT==1),] %>% group_by(treatment,race_re) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))



matched.data[which(matched.data$CT==0),] %>% group_by(treatment) %>% summarise_at(vars(age),list(age_mean=mean,age_se=sd))

matched.data[which(matched.data$CT==0),] %>% group_by(treatment,gender) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))

matched.data[which(matched.data$CT==0),] %>% group_by(treatment,eth) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))

matched.data[which(matched.data$CT==0),] %>% group_by(treatment,karno) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))

matched.data$rdri_re <- ifelse(matched.data$red_dis_int==1,1,
                           ifelse(matched.data$red_dis_high==1,2,0))


matched.data[which(matched.data$CT==0),] %>% group_by(treatment,rdri_re) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))

matched.data[which(matched.data$CT==0),] %>% group_by(treatment,cmv) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))


matched.data[which(matched.data$CT==0),] %>% group_by(treatment,hct_re) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))

matched.data$race_re <-  ifelse(matched.data$raceb==1,1,
                            ifelse(matched.data$raceo==1,2,0))

matched.data[which(matched.data$CT==0),] %>% group_by(treatment,race_re) %>%
  summarise(count=n()) %>%
  mutate(freq = count/sum(count))


