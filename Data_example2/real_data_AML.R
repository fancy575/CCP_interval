library(Rcpp)

library(survival)
library(MatchIt)
library(R2jags)

raw_dt <- read.csv("BMT_data.csv")
#raw_dt <- raw_dt[which(raw_dt$disease == 50 & raw_dt$donorgp != 3 ),]

raw_dt <- raw_dt[which(raw_dt$disease==10),] 
raw_dt <- raw_dt[which(!raw_dt$distatus %in% c(14,15) ),]

raw_dt$source <- ifelse(raw_dt$strata %in% c(1,2), 1,0) # 1: CT, 0: HD
raw_dt$treatment <- ifelse(raw_dt$strata %in% c(1,4), 1,0  ) # 1: RIC, 0: MAC
table(raw_dt$source)

raw_dt$agegp_new <- ifelse(raw_dt$agegp %in% c(2,3,4),0,1 ) #0 < 50 ,  1 >=50
raw_dt$karnofcat_new <- ifelse(raw_dt$karnofcat == 1, 1,
                               ifelse(raw_dt$karnofcat == 2, 0,NA))
raw_dt$graftypecat_new1 <- ifelse(raw_dt$graftypecat==1,1,ifelse(raw_dt$graftypecat==99,NA,0)) # Bone marrow
raw_dt$graftypecat_new2 <- ifelse(raw_dt$graftypecat==2,1,ifelse(raw_dt$graftypecat==99,NA,0)) # peripheral, cord blood as reference


raw_dt$donorgp_new1 <- ifelse(raw_dt$donorgp==1,1,0) # HLA
raw_dt$donorgp_new2<- ifelse(raw_dt$donorgp==4,1,0) # well matched, 3 and 5 Other related/Partially-matched unrelated  as reference

raw_dt$atgcampathgp_new <- ifelse(raw_dt$atgcampathgp==2,1,0) # No ATG or CAMPATH as reference

raw_dt$distatus_new <- ifelse(raw_dt$distatus == 12,1,0) # CR1 as reference

raw_dt$cytogene_new1 <- ifelse(raw_dt$cytogene %in% c(1,2,3),1,0 ) 
raw_dt$cytogene_new2 <- ifelse(raw_dt$cytogene ==4,1,0 )  ## unknown as reference


raw_dt$CT <- ifelse(raw_dt$source==1,1,0)

ana_data <- raw_dt[,c(56,45:55,25,26,13)]

ana_data <- as.matrix(ana_data)
ana_data[which(ana_data==99)] <- NA
ana_data <- na.omit(ana_data)
ana_data <- as.data.frame(ana_data)


set.seed(1)
match_form <- as.formula("CT~agegp_new+karnofcat_new+graftypecat_new1+
                         donorgp_new1+donorgp_new2 + atgcampathgp_new+distatus_new+cytogene_new1 +cytogene_new1
                         ")


m.out <- matchit(match_form,data=ana_data,caliper = 0.2,ratio=2)
matched.data <- match.data(m.out,drop.unmatched = F)

K = 10 #Number of intervals for alpha
kappa=c(0,quantile(matched.data$intxrel[which(matched.data$dfs==1 & matched.data$CT==1)],
                   seq(0,1,by=1/K))[-c(1,K+1)],
        max(matched.data$intxrel)+0.0001) 

hc_interval <- c(0,3.29,6,12,24,36,48,60,Inf)
ct_interval <- c(c(seq(7,91,7),100)/30.4,6,12,18, Inf)

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


## Fit combine without interaction
data_all <- matched.data[order(matched.data$censor),]
idxL_all = max(which(data_all$censor==0))
idxI_all = c(min(which(data_all$censor==1)),max(which(data_all $censor==1)))
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))

X_all <- data_all[,c(1:5,7:12,16)]
p_all = ncol(X_all)
L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1


model.data <- list("idxL_all","idxI_all","idxR_all","idxE_all","lower.int.obs_all","upper.int.obs_all","X_all","p_all","nint_all","kappa","R_all","L_all","ones_all",
                   "n_all")

model.parameters <- c("beta","alpha","nu") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_no_int_combine.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_no_int <- model.fit.sep$BUGSoutput$summary[11:22,c(1,3,7)]
rownames(survival_all_no_int) <- colnames(X_all)
survival_all_no_int <- cbind(survival_all_no_int,
                             paste(round(survival_all_no_int[,1],2)," (",
                                   round(survival_all_no_int[,2],2),",",
                                   round(survival_all_no_int[,3],2),")",sep="")
)

write.csv(survival_all_no_int,file="survival_AML_combine_no_int.csv")



### midpoint


data_all <- matched.data[order(matched.data$censor),]

data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )] <- 
  (data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] +
     data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )])/2

data_all$censor[which(data_all$censor==0 | data_all$censor==1 )] <- 3

data_all$L_int <- findInterval(data_all$obs_L,kappa)
data_all$R_int <- findInterval(data_all$obs_R,kappa)


idxL_all = 0
idxI_all = c(0,0)
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))

X_all <- data_all[,c(1:5,7:12,16)]
p_all = ncol(X_all)
L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1


model.data <- list("idxL_all","idxI_all","idxR_all","idxE_all","lower.int.obs_all","upper.int.obs_all","X_all","p_all","nint_all","kappa","R_all","L_all","ones_all",
                   "n_all")

model.parameters <- c("beta","alpha","nu") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_no_int_combine_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_no_int <- model.fit.sep$BUGSoutput$summary[11:22,c(1,3,7)]
rownames(survival_all_no_int) <- colnames(X_all)
survival_all_no_int <- cbind(survival_all_no_int,
                             paste(round(survival_all_no_int[,1],2)," (",
                                   round(survival_all_no_int[,2],2),",",
                                   round(survival_all_no_int[,3],2),")",sep="")
)

write.csv(survival_all_no_int,file="survival_AML_combine_no_int_midpoint.csv")


cox_dt <- data_all
cox_dt$event <- ifelse(data_all$censor==2,0,1)
cox_fit <- coxph(Surv(obs_L,event)~CT+treatment + agegp_new + karnofcat_new + 
                   graftypecat_new1 +  donorgp_new1 + 
                   donorgp_new2 + atgcampathgp_new + distatus_new + 
                   cytogene_new1 + cytogene_new2,data=cox_dt)


cox_all_no_int <- data.frame(estimate =coef(cox_fit),confint = paste(round(coef(cox_fit),2)," (",
                                                                       round(confint(cox_fit)[,1],2),",",
                                                                       round(confint(cox_fit)[,2],2),")",sep="" )
)

write.csv(cox_all_no_int, file="cox_combine_no_int_midpoint.csv")

### right side


data_all <- matched.data[order(matched.data$censor),]

data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )]

data_all$censor[which(data_all$censor==0 | data_all$censor==1 )] <- 3

data_all$L_int <- findInterval(data_all$obs_L,kappa)
data_all$R_int <- findInterval(data_all$obs_R,kappa)


idxL_all = 0
idxI_all = c(0,0)
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))

X_all <- data_all[,c(1:5,7:12,16)]
p_all = ncol(X_all)
L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1


model.data <- list("idxL_all","idxI_all","idxR_all","idxE_all","lower.int.obs_all","upper.int.obs_all","X_all","p_all","nint_all","kappa","R_all","L_all","ones_all",
                   "n_all")

model.parameters <- c("beta","alpha","nu") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_no_int_combine_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_no_int <- model.fit.sep$BUGSoutput$summary[11:22,c(1,3,7)]
rownames(survival_all_no_int) <- colnames(X_all)
survival_all_no_int <- cbind(survival_all_no_int,
                             paste(round(survival_all_no_int[,1],2)," (",
                                   round(survival_all_no_int[,2],2),",",
                                   round(survival_all_no_int[,3],2),")",sep="")
)

write.csv(survival_all_no_int,file="survival_AML_combine_no_int_right.csv")


cox_dt <- data_all
cox_dt$event <- ifelse(data_all$censor==2,0,1)
cox_fit <- coxph(Surv(obs_L,event)~CT+treatment + agegp_new + karnofcat_new + 
                   graftypecat_new1  + donorgp_new1 + 
                   donorgp_new2 + atgcampathgp_new + distatus_new + 
                   cytogene_new1 + cytogene_new2,data=cox_dt)


cox_all_no_int <- data.frame(estimate =coef(cox_fit),confint = paste(round(coef(cox_fit),2)," (",
                                                                     round(confint(cox_fit)[,1],2),",",
                                                                     round(confint(cox_fit)[,2],2),")",sep="" )
)

write.csv(cox_all_no_int, file="cox_combine_no_int_right.csv")



# Fit combine with interaction

data_all <- matched.data[order(matched.data$censor),]


idxL_all = max(which(data_all$censor==0))
idxI_all = c(min(which(data_all$censor==1)),max(which(data_all $censor==1)))
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))




X_all <- data_all[,c(1:12,16)]

X_all$RIC_CT <- ifelse(X_all$CT==1 & X_all$treatment==1, 1, 0)
X_all$RIC_HD <- ifelse(X_all$CT==0 & X_all$treatment==1, 1, 0)
X_all$MAC_CT <- ifelse(X_all$CT==1 & X_all$treatment==0, 1, 0)
X_all <- X_all[c(14:16,3:5,7:12,13)]
p_all = ncol(X_all)

L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1



model.parameters <- c("beta","alpha","beta_diff","beta_ric_ct") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_int_combine.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_int <- model.fit.sep$BUGSoutput$summary[11:25,c(1,3,7)]
rownames(survival_all_int) <- c(colnames(X_all),"trt_diff","beta_ric_ct")
survival_all_int <- cbind(survival_all_int,
                          paste(round(survival_all_int[,1],2)," (",
                                round(survival_all_int[,2],2),",",
                                round(survival_all_int[,3],2),")",sep="")
)

write.csv(survival_all_int,file="survival_AML_combine_int.csv")


## midpoint
data_all <- matched.data[order(matched.data$censor),]

data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )] <- 
  (data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] +
     data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )])/2

data_all$censor[which(data_all$censor==0 | data_all$censor==1 )] <- 3

data_all$L_int <- findInterval(data_all$obs_L,kappa)
data_all$R_int <- findInterval(data_all$obs_R,kappa)


idxL_all = 0
idxI_all = c(0,0)
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))




X_all <- data_all[,c(1:12,16)]

X_all$RIC_CT <- ifelse(X_all$CT==1 & X_all$treatment==1, 1, 0)
X_all$RIC_HD <- ifelse(X_all$CT==0 & X_all$treatment==1, 1, 0)
X_all$MAC_CT <- ifelse(X_all$CT==1 & X_all$treatment==0, 1, 0)
X_all <- X_all[c(14:16,3:5,7:12,13)]
p_all = ncol(X_all)

L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1

model.parameters <- c("beta","alpha","beta_diff","beta_ric_ct") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_int_combine_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_int <- model.fit.sep$BUGSoutput$summary[11:25,c(1,3,7)]
rownames(survival_all_int) <- c(colnames(X_all),"trt_diff","beta_ric_ct")
survival_all_int <- cbind(survival_all_int,
                          paste(round(survival_all_int[,1],2)," (",
                                round(survival_all_int[,2],2),",",
                                round(survival_all_int[,3],2),")",sep="")
)

write.csv(survival_all_int,file="survival_AML_combine_int_midpoint.csv")

## right

data_all <- matched.data[order(matched.data$censor),]

data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )]

data_all$censor[which(data_all$censor==0 | data_all$censor==1 )] <- 3

data_all$L_int <- findInterval(data_all$obs_L,kappa)
data_all$R_int <- findInterval(data_all$obs_R,kappa)


idxL_all = 0
idxI_all = c(0,0)
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))




X_all <- data_all[,c(1:12,16)]

X_all$RIC_CT <- ifelse(X_all$CT==1 & X_all$treatment==1, 1, 0)
X_all$RIC_HD <- ifelse(X_all$CT==0 & X_all$treatment==1, 1, 0)
X_all$MAC_CT <- ifelse(X_all$CT==1 & X_all$treatment==0, 1, 0)
X_all <- X_all[c(14:16,3:5,7:12,13)]
p_all = ncol(X_all)

L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1

model.parameters <- c("beta","alpha","beta_diff","beta_ric_ct") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_int_combine_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_int <- model.fit.sep$BUGSoutput$summary[11:25,c(1,3,7)]
rownames(survival_all_int) <- c(colnames(X_all),"trt_diff","beta_ric_ct")
survival_all_int <- cbind(survival_all_int,
                          paste(round(survival_all_int[,1],2)," (",
                                round(survival_all_int[,2],2),",",
                                round(survival_all_int[,3],2),")",sep="")
)

write.csv(survival_all_int,file="survival_AML_combine_int_right.csv")


## Matched data 

matched.data <- matched.data[-which(is.na(matched.data$subclass) & matched.data$CT==0 ),]

matched.data$subclass <- as.numeric(levels(matched.data$subclass))[matched.data$subclass]
matched.data$subclass[which(is.na(matched.data$subclass))] <- (max(as.numeric(matched.data$subclass),na.rm=T)+1) : (sum(is.na(matched.data$subclass)) + max(as.numeric(matched.data$subclass),na.rm=T))
matched.data <- matched.data[order(matched.data$subclass),]

# without interaction
# interval censor

data_all <- matched.data[order(matched.data$censor),]
idxL_all = max(which(data_all$censor==0))
idxI_all = c(min(which(data_all$censor==1)),max(which(data_all $censor==1)))
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))

X_all <- data_all[,c(1:5,7:12,16)]
p_all = ncol(X_all)
L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1


model.data <- list("idxL_all","idxI_all","idxR_all","idxE_all","lower.int.obs_all","upper.int.obs_all","X_all","p_all","nint_all","kappa","R_all","L_all","ones_all",
                   "n_all")

model.parameters <- c("beta","alpha","nu") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_no_int_combine.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_no_int <- model.fit.sep$BUGSoutput$summary[11:22,c(1,3,7)]
rownames(survival_all_no_int) <- colnames(X_all)
survival_all_no_int <- cbind(survival_all_no_int,
                             paste(round(survival_all_no_int[,1],2)," (",
                                   round(survival_all_no_int[,2],2),",",
                                   round(survival_all_no_int[,3],2),")",sep="")
)

write.csv(survival_all_no_int,file="survival_matched_no_int.csv")




# midpoint

data_all <- matched.data[order(matched.data$censor),]

data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )] <- 
  (data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] +
     data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )])/2

data_all$censor[which(data_all$censor==0 | data_all$censor==1 )] <- 3

data_all$L_int <- findInterval(data_all$obs_L,kappa)
data_all$R_int <- findInterval(data_all$obs_R,kappa)


idxL_all = 0
idxI_all = c(0,0)
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))

X_all <- data_all[,c(1:5,7:12,16)]
p_all = ncol(X_all)
L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1



model.data <- list("idxL_all","idxI_all","idxR_all","idxE_all","lower.int.obs_all","upper.int.obs_all","X_all","p_all","nint_all","kappa","R_all","L_all","ones_all",
                   "n_all")

model.parameters <- c("beta","alpha","nu") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_no_int_combine_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_no_int <- model.fit.sep$BUGSoutput$summary[11:22,c(1,3,7)]
rownames(survival_all_no_int) <- colnames(X_all)
survival_all_no_int <- cbind(survival_all_no_int,
                             paste(round(survival_all_no_int[,1],2)," (",
                                   round(survival_all_no_int[,2],2),",",
                                   round(survival_all_no_int[,3],2),")",sep="")
)

write.csv(survival_all_no_int,file="survival_matched_no_int_midpoint.csv")

cox_dt <- data_all
cox_dt$event <- ifelse(data_all$censor==2,0,1)
cox_fit <- coxph(Surv(obs_L,event)~CT+treatment + agegp_new + karnofcat_new + 
                   graftypecat_new1  + donorgp_new1 + 
                   donorgp_new2 + atgcampathgp_new + distatus_new + 
                   cytogene_new1 + cytogene_new2+distance,data=cox_dt)


cox_all_no_int <- data.frame(estimate =coef(cox_fit),confint = paste(round(coef(cox_fit),2)," (",
                                                                     round(confint(cox_fit)[,1],2),",",
                                                                     round(confint(cox_fit)[,2],2),")",sep="" )
)

write.csv(cox_all_no_int, file="cox_matched_no_int_midpoint.csv")

# right

data_all <- matched.data[order(matched.data$censor),]

data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )]

data_all$censor[which(data_all$censor==0 | data_all$censor==1 )] <- 3

data_all$L_int <- findInterval(data_all$obs_L,kappa)
data_all$R_int <- findInterval(data_all$obs_R,kappa)


idxL_all = 0
idxI_all = c(0,0)
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))

X_all <- data_all[,c(1:5,7:12,16)]
p_all = ncol(X_all)
L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1



model.data <- list("idxL_all","idxI_all","idxR_all","idxE_all","lower.int.obs_all","upper.int.obs_all","X_all","p_all","nint_all","kappa","R_all","L_all","ones_all",
                   "n_all")

model.parameters <- c("beta","alpha","nu") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_no_int_combine_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_no_int <- model.fit.sep$BUGSoutput$summary[11:22,c(1,3,7)]
rownames(survival_all_no_int) <- colnames(X_all)
survival_all_no_int <- cbind(survival_all_no_int,
                             paste(round(survival_all_no_int[,1],2)," (",
                                   round(survival_all_no_int[,2],2),",",
                                   round(survival_all_no_int[,3],2),")",sep="")
)

write.csv(survival_all_no_int,file="survival_matched_no_int_right.csv")


cox_dt <- data_all
cox_dt$event <- ifelse(data_all$censor==2,0,1)
cox_fit <- coxph(Surv(obs_L,event)~CT+treatment + agegp_new + karnofcat_new + 
                   graftypecat_new1  + donorgp_new1 + 
                   donorgp_new2 + atgcampathgp_new + distatus_new + 
                   cytogene_new1 + cytogene_new2+distance,data=cox_dt)


cox_all_no_int <- data.frame(estimate =coef(cox_fit),confint = paste(round(coef(cox_fit),2)," (",
                                                                     round(confint(cox_fit)[,1],2),",",
                                                                     round(confint(cox_fit)[,2],2),")",sep="" )
)

write.csv(cox_all_no_int, file="cox_matched_no_int_right.csv")


# Fit matched with interaction

data_all <- matched.data[order(matched.data$censor),]


idxL_all = max(which(data_all$censor==0))
idxI_all = c(min(which(data_all$censor==1)),max(which(data_all $censor==1)))
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))

L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1




X_all <- data_all[,c(1:12,16)]

X_all$trt_CT <- ifelse(X_all$CT==1 & X_all$treatment==1, 1, 0)
X_all$trt_HD <- ifelse(X_all$CT==0 & X_all$treatment==1, 1, 0)
X_all$notrt_CT <- ifelse(X_all$CT==1 & X_all$treatment==0, 1, 0)
X_all <- X_all[c(14:16,3:5,7:12,13)]
p_all = ncol(X_all)

model.parameters <- c("beta","alpha","beta_diff","beta_ric_ct") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_int_combine.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_int <- model.fit.sep$BUGSoutput$summary[11:25,c(1,3,7)]
rownames(survival_all_int) <- c(colnames(X_all),"trt_diff","beta_ric_ct")
survival_all_int <- cbind(survival_all_int,
                          paste(round(survival_all_int[,1],2)," (",
                                round(survival_all_int[,2],2),",",
                                round(survival_all_int[,3],2),")",sep="")
)

write.csv(survival_all_int,file="survival_matched_all_int.csv")



## midpoint
data_all <- matched.data[order(matched.data$censor),]

data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )] <- 
  (data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] +
     data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )])/2

data_all$censor[which(data_all$censor==0 | data_all$censor==1 )] <- 3

data_all$L_int <- findInterval(data_all$obs_L,kappa)
data_all$R_int <- findInterval(data_all$obs_R,kappa)


idxL_all = 0
idxI_all = c(0,0)
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))

L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int




X_all <- data_all[,c(1:12,16)]

X_all$trt_CT <- ifelse(X_all$CT==1 & X_all$treatment==1, 1, 0)
X_all$trt_HD <- ifelse(X_all$CT==0 & X_all$treatment==1, 1, 0)
X_all$notrt_CT <- ifelse(X_all$CT==1 & X_all$treatment==0, 1, 0)
X_all <- X_all[c(14:16,3:5,7:12,13)]
p_all = ncol(X_all)
L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1

model.parameters <- c("beta","alpha","beta_diff","beta_ric_ct") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_int_combine_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_int <- model.fit.sep$BUGSoutput$summary[11:25,c(1,3,7)]
rownames(survival_all_int) <- c(colnames(X_all),"trt_diff","beta_ric_ct")
survival_all_int <- cbind(survival_all_int,
                          paste(round(survival_all_int[,1],2)," (",
                                round(survival_all_int[,2],2),",",
                                round(survival_all_int[,3],2),")",sep="")
)

write.csv(survival_all_int,file="survival_matched_int_midpoint.csv")

## right

data_all <- matched.data[order(matched.data$censor),]

data_all$obs_L[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )] <- 
  data_all$obs_R[which(data_all$censor==0 | data_all$censor==1 )]

data_all$censor[which(data_all$censor==0 | data_all$censor==1 )] <- 3

data_all$L_int <- findInterval(data_all$obs_L,kappa)
data_all$R_int <- findInterval(data_all$obs_R,kappa)


idxL_all = 0
idxI_all = c(0,0)
idxR_all = max(which(data_all$censor==2))
idxE_all = max(which(data_all$censor==3))




X_all <- data_all[,c(1:12,16)]

X_all$trt_CT <- ifelse(X_all$CT==1 & X_all$treatment==1, 1, 0)
X_all$trt_HD <- ifelse(X_all$CT==0 & X_all$treatment==1, 1, 0)
X_all$notrt_CT <- ifelse(X_all$CT==1 & X_all$treatment==0, 1, 0)
X_all <- X_all[c(14:16,3:5,7:12,13)]
p_all = ncol(X_all)
L_all = data_all$obs_L
R_all = data_all$obs_R
lower.int.obs_all = data_all$L_int
upper.int.obs_all = data_all$R_int

ones_all <- rep(1,nrow(data_all))
n_all <- length(ones_all)
nint_all = length(kappa)-1

model.parameters <- c("beta","alpha","beta_diff","beta_ric_ct") # Murray's Method

model.fit.sep <- jags(model.file= "exact_likelihood_no_random_int_combine_imp.txt", data=model.data,
                      parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_all_int <- model.fit.sep$BUGSoutput$summary[11:25,c(1,3,7)]
rownames(survival_all_int) <- c(colnames(X_all),"trt_diff","beta_ric_ct")
survival_all_int <- cbind(survival_all_int,
                          paste(round(survival_all_int[,1],2)," (",
                                round(survival_all_int[,2],2),",",
                                round(survival_all_int[,3],2),")",sep="")
)

write.csv(survival_all_int,file="survival_matched_int_right.csv")


## matched random effect no interaction

# interval censor
data_all <- matched.data[order(matched.data$censor),]


out.H=data_all[data_all$CT==0,];out.C=data_all[data_all$CT==1,]


## HD data
out.H <- out.H[order(out.H$censor),]
idxL = max(which(out.H$censor==0))
idxI = c(min(which(out.H$censor==1)),max(which(out.H$censor==1)))
idxR = max(which(out.H$censor==2))
idxE = max(which(out.H$censor==3))

X <- out.H[,c(1:5,7:12)]
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

X.C <- out.C[,c(1:5,7:12)]
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

model.parameters <- c("beta","alpha","nu") # Murray's Method


model.fit.cluster <- jags(model.file= "exact_likelihood_random_combine_no_int.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_murrray_parm <- model.fit.cluster$BUGSoutput$summary[11:21,c(1,3,7)]
rownames(survival_murrray_parm) <- colnames(X.C)
survival_murrray_parm <- cbind(survival_murrray_parm,
                               paste(round(survival_murrray_parm[,1],2)," (",
                                     round(survival_murrray_parm[,2],2),",",
                                     round(survival_murrray_parm[,3],2),")",sep="")
)

write.csv(survival_murrray_parm,file="survival_matched_random_no_int.csv")


# midpoint
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

X <- out.H[,c(1:5,7:12)]
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

X.C <- out.C[,c(1:5,7:12)]
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

model.parameters <- c("beta","alpha","nu") 


model.fit.cluster <- jags(model.file= "exact_likelihood_random_combine_no_int_imp.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_murrray_parm <- model.fit.cluster$BUGSoutput$summary[11:21,c(1,3,7)]
rownames(survival_murrray_parm) <- colnames(X.C)
survival_murrray_parm <- cbind(survival_murrray_parm,
                               paste(round(survival_murrray_parm[,1],2)," (",
                                     round(survival_murrray_parm[,2],2),",",
                                     round(survival_murrray_parm[,3],2),")",sep="")
)

write.csv(survival_murrray_parm,file="survival_matched_random_no_int_midpoint.csv")


# right
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

X <- out.H[,c(1:5,7:12)]
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

X.C <- out.C[,c(1:5,7:12)]
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

model.parameters <- c("beta","alpha","nu") 


model.fit.cluster <- jags(model.file= "exact_likelihood_random_combine_no_int_imp.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_murrray_parm <- model.fit.cluster$BUGSoutput$summary[11:21,c(1,3,7)]
rownames(survival_murrray_parm) <- colnames(X.C)
survival_murrray_parm <- cbind(survival_murrray_parm,
                               paste(round(survival_murrray_parm[,1],2)," (",
                                     round(survival_murrray_parm[,2],2),",",
                                     round(survival_murrray_parm[,3],2),")",sep="")
)

write.csv(survival_murrray_parm,file="survival_matched_random_no_int_right.csv")

## matched random effect interaction

# interval censor


data_all <- matched.data[order(matched.data$censor),]


out.H=data_all[data_all$CT==0,];out.C=data_all[data_all$CT==1,]


## HD data
out.H <- out.H[order(out.H$censor),]
idxL = max(which(out.H$censor==0))
idxI = c(min(which(out.H$censor==1)),max(which(out.H$censor==1)))
idxR = max(which(out.H$censor==2))
idxE = max(which(out.H$censor==3))

X <- out.H[,c(1:12)]
X$trt_CT <- ifelse(X$CT==1 & X$treatment==1, 1, 0)
X$trt_HD <- ifelse(X$CT==0 & X$treatment==1, 1, 0)
X$notrt_CT <- ifelse(X$CT==1 & X$treatment==0, 1, 0)
X <- X[c(13:15,3:5,7:12)]


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

X.C <- out.C[,c(1:12)]
X.C$trt_CT <- ifelse(X.C$CT==1 & X.C$treatment==1, 1, 0)
X.C$trt_HD <- ifelse(X.C$CT==0 & X.C$treatment==1, 1, 0)
X.C$notrt_CT <- ifelse(X.C$CT==1 & X.C$treatment==0, 1, 0)
X.C <- X.C[c(13:15,3:5,7:12)]

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

model.parameters <- c("beta","alpha","nu","beta_diff","beta_ric_ct") # Murray's Method


model.fit.cluster <- jags(model.file= "exact_likelihood_random_combine_int.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_murrray_parm <- model.fit.cluster$BUGSoutput$summary[11:24,c(1,3,7)]
rownames(survival_murrray_parm) <- c(colnames(X.C),"trt_diff","beta_ric_ct")
survival_murrray_parm <- cbind(survival_murrray_parm,
                               paste(round(survival_murrray_parm[,1],2)," (",
                                     round(survival_murrray_parm[,2],2),",",
                                     round(survival_murrray_parm[,3],2),")",sep="")
)

write.csv(survival_murrray_parm,file="survival_matched_random_int.csv")



# midpoint


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

X <- out.H[,c(1:12)]
X$trt_CT <- ifelse(X$CT==1 & X$treatment==1, 1, 0)
X$trt_HD <- ifelse(X$CT==0 & X$treatment==1, 1, 0)
X$notrt_CT <- ifelse(X$CT==1 & X$treatment==0, 1, 0)
X <- X[c(13:15,3:5,7:12)]


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

X.C <- out.C[,c(1:12)]
X.C$trt_CT <- ifelse(X.C$CT==1 & X.C$treatment==1, 1, 0)
X.C$trt_HD <- ifelse(X.C$CT==0 & X.C$treatment==1, 1, 0)
X.C$notrt_CT <- ifelse(X.C$CT==1 & X.C$treatment==0, 1, 0)
X.C <- X.C[c(13:15,3:5,7:12)]

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

model.parameters <- c("beta","alpha","nu","beta_diff","beta_ric_ct") # Murray's Method


model.fit.cluster <- jags(model.file= "exact_likelihood_random_combine_int_imp.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_murrray_parm <- model.fit.cluster$BUGSoutput$summary[11:24,c(1,3,7)]
rownames(survival_murrray_parm) <- c(colnames(X.C),"trt_diff","beta_ric_ct")
survival_murrray_parm <- cbind(survival_murrray_parm,
                               paste(round(survival_murrray_parm[,1],2)," (",
                                     round(survival_murrray_parm[,2],2),",",
                                     round(survival_murrray_parm[,3],2),")",sep="")
)

write.csv(survival_murrray_parm,file="survival_matched_random_int_midpoint.csv")





# right


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

X <- out.H[,c(1:12)]
X$trt_CT <- ifelse(X$CT==1 & X$treatment==1, 1, 0)
X$trt_HD <- ifelse(X$CT==0 & X$treatment==1, 1, 0)
X$notrt_CT <- ifelse(X$CT==1 & X$treatment==0, 1, 0)
X <- X[c(13:15,3:5,7:12)]


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

X.C <- out.C[,c(1:12)]
X.C$trt_CT <- ifelse(X.C$CT==1 & X.C$treatment==1, 1, 0)
X.C$trt_HD <- ifelse(X.C$CT==0 & X.C$treatment==1, 1, 0)
X.C$notrt_CT <- ifelse(X.C$CT==1 & X.C$treatment==0, 1, 0)
X.C <- X.C[c(13:15,3:5,7:12)]

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

model.parameters <- c("beta","alpha","nu","beta_diff","beta_ric_ct") # Murray's Method


model.fit.cluster <- jags(model.file= "exact_likelihood_random_combine_int_imp.txt", data=model.data,
                          parameters=model.parameters, n.chains = 2,n.iter=5000,n.thin=5,n.burnin=1000)

survival_murrray_parm <- model.fit.cluster$BUGSoutput$summary[11:24,c(1,3,7)]
rownames(survival_murrray_parm) <- c(colnames(X.C),"trt_diff","beta_ric_ct")
survival_murrray_parm <- cbind(survival_murrray_parm,
                               paste(round(survival_murrray_parm[,1],2)," (",
                                     round(survival_murrray_parm[,2],2),",",
                                     round(survival_murrray_parm[,3],2),")",sep="")
)

write.csv(survival_murrray_parm,file="survival_matched_random_int_right.csv")




