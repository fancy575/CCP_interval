setwd("/scratch/g/skim/CT_HD_Bayesian_ones(1228)/X_dependent_5_interval/CP/combine")
summary_func <- function(rds_data,true_value){
  raw_dt <- rds_data[[1]]$BUGSoutput$summary[,1]
  raw_25 <- rds_data[[1]]$BUGSoutput$summary[,3]
  raw_75 <- rds_data[[1]]$BUGSoutput$summary[,7]
  
  
  for (i in 2: length(rds_data)) {
    raw_dt <- rbind(raw_dt,rds_data[[i]]$BUGSoutput$summary[,1])
    raw_25 <- rbind(raw_25,rds_data[[i]]$BUGSoutput$summary[,3])
    raw_75 <- rbind(raw_75,rds_data[[i]]$BUGSoutput$summary[,7])
  }
  raw_result <- raw_dt[,c(2:29,45:48)]
  raw_ci25 <- raw_25[,c(2:29,45:48)]
  raw_ci75 <- raw_75[,c(2:29,45:48)]
  
  
  remove_id <- which(abs(raw_result[,29])>2)
  
  if(length(remove_id)>0){
    raw_result <- raw_result[-remove_id,]
    raw_ci25 <- raw_ci25[-remove_id,]
    raw_ci75 <- raw_ci75[-remove_id,]
  }
  
  
  raw_sum <- data.frame(true=c(rep(true_value[1],5),rep(true_value[2],5),true_value[3:22],NA,NA),
                        mean = apply(raw_result, 2, mean),
                        sd = apply(raw_result,2,sd))
  
  raw_mse <- sweep(raw_result[,1:30],2,c(rep(true_value[1],5),rep(true_value[2],5),true_value[3:22]),"-")
  raw_se_margin <- apply(raw_result[,1:30],2,sd)/sqrt(length(rds_data))*qnorm(0.975)
  raw_mse <- raw_mse^2
  raw_mse <- apply(raw_mse, 2, mean)
  
  raw_sum$mse <- c(raw_mse,NA,NA)
  raw_sum$se_margin <- c(raw_se_margin,NA,NA)
  
  true_mt <- t(replicate(length(rds_data)-length(remove_id),c(rep(true_value[1],5),rep(true_value[2],5),true_value[3:22])))
  ci <- ifelse(raw_ci25[,1:30] <= true_mt & true_mt <= raw_ci75[,1:30],1,0 )
  raw_sum$cr <- c(apply(ci,2,mean),NA,NA)
  
  raw_sum <- raw_sum[c(1:10,29,11:28,30:32),]
  raw_sum[c(31:32),c(1,3:5)] <- ""
  
  
  return(raw_sum)
  
  
}







### Fully exchangeable
fullyex <- readRDS("fully_ex1.rds")
fullyex_true <- c(1,1, 
                  c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                  c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                  0,1)

fullyex_sum <- summary_func(fullyex,fullyex_true)
write.csv(fullyex_sum,file="fullyex1_sum.csv")


fullyex <- readRDS("fully_ex2.4.rds")

fullyex_true <- c(2.4,2.4, 
                  c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                  c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                  0,1)

fullyex_sum <- summary_func(fullyex,fullyex_true)
write.csv(fullyex_sum,file="fullyex2.4_sum.csv")


# nonfully 0.2

non02 <- readRDS("nonfully1_02.rds")

non02_true <- c(1.2,1, 
                c(1,-1,0.5,-0.5,0,0.2,1.2,-0.8,0.7),
                c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                0,1)

non02_sum <- summary_func(non02,non02_true)
write.csv(non02_sum,file="non1_02_sum.csv")


non02 <- readRDS("nonfully2.4_02.rds")


non02_true <- c(2.6,2.4, 
                c(1,-1,0.5,-0.5,0,0.2,1.2,-0.8,0.7),
                c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                0,1)

non02_sum <- summary_func(non02,non02_true)
write.csv(non02_sum,file="non2.4_02_sum.csv")




# nonfully 0.2n

non02n <- readRDS("nonfully1_02n.rds")


non02n_true <- c(0.8,1, 
                c(1,-1,0.5,-0.5,0,-0.2,0.8,-1.2,0.3),
                c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                0,1)

non02n_sum <- summary_func(non02n,non02n_true)
write.csv(non02n_sum,file="non1_02n_sum.csv")


non02n <- readRDS("nonfully2.4_02n.rds")


non02n_true <- c(2.2,2.4, 
                c(1,-1,0.5,-0.5,0,-0.2,0.8,-1.2,0.3),
                c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                0,1)

non02n_sum <- summary_func(non02n,non02n_true)
write.csv(non02n_sum,file="non2.4_02n_sum.csv")



# nonfully 0.4

non04 <- readRDS("nonfully1_04.rds")


non04_true <- c(1.4,1, 
                c(1,-1,0.5,-0.5,0,0.4,1.4,-0.6,0.9),
                c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                0,1)

non04_sum <- summary_func(non04,non04_true)
write.csv(non04_sum,file="non1_04_sum.csv")


non04 <- readRDS("nonfully2.4_04.rds")



non04_true <- c(2.8,2.4, 
                c(1,-1,0.5,-0.5,0,0.4,1.4,-0.6,0.9),
                c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                0,1)

non04_sum <- summary_func(non04,non04_true)
write.csv(non04_sum,file="non2.4_04_sum.csv")




# nonfully 0.8

non08 <- readRDS("nonfully1_08.rds")



non08_true <- c(1.8,1, 
                c(1,-1,0.5,-0.5,0,0.8,1.8,-0.2,1.3),
                c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                0,1)

non08_sum <- summary_func(non08,non08_true)
write.csv(non08_sum,file="non1_08_sum.csv")


non08 <- readRDS("nonfully2.4_08.rds")




non08_true <- c(3.2,2.4, 
                c(1,-1,0.5,-0.5,0,0.8,1.8,-0.2,1.3),
                c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                0,1)

non08_sum <- summary_func(non08,non08_true)
write.csv(non08_sum,file="non2.4_08_sum.csv")


# nonfully 1.2

non1.2 <- readRDS("nonfully1_1.2.rds")


non1.2_true <- c(2.2,1, 
                c(1,-1,0.5,-0.5,0,1.2,2.2,0.2,1.7),
                c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                0,1)

non1.2_sum <- summary_func(non1.2,non1.2_true)
write.csv(non1.2_sum,file="non1_1.2_sum.csv")


non1.2 <- readRDS("nonfully2.4_1.2.rds")



non1.2_true <- c(3.6,2.4, 
                c(1,-1,0.5,-0.5,0,1.2,2.2,0.2,1.7),
                c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                0,1)

non1.2_sum <- summary_func(non1.2,non1.2_true)
write.csv(non1.2_sum,file="non2.4_1.2_sum.csv")



# nonfully 1.6

non1.6 <- readRDS("nonfully1_1.6.rds")



non1.6_true <- c(2.6,1, 
                 c(1,-1,0.5,-0.5,0,1.6,2.6,0.6,2.1),
                 c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                 0,1)

non1.6_sum <- summary_func(non1.6,non1.6_true)
write.csv(non1.6_sum,file="non1_1.6_sum.csv")


non1.6 <- readRDS("nonfully2.4_1.6.rds")


non1.6_true <- c(4.0,2.4, 
                 c(1,-1,0.5,-0.5,0,1.6,2.6,0.6,2.1),
                 c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                 0,1)

non1.6_sum <- summary_func(non1.6,non1.6_true)
write.csv(non1.6_sum,file="non2.4_1.6_sum.csv")


# nonfully 2.0

non2.0 <- readRDS("nonfully1_2.0.rds")



non2.0_true <- c(3.0,1, 
                 c(1,-1,0.5,-0.5,0,2.0,3.0,1.0,2.5),
                 c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                 0,1)

non2.0_sum <- summary_func(non2.0,non2.0_true)
write.csv(non2.0_sum,file="non1_2.0_sum.csv")


non2.0 <- readRDS("nonfully2.4_2.0.rds")


non2.0_true <- c(4.4,2.4, 
                 c(1,-1,0.5,-0.5,0,2.0,3.0,1.0,2.5),
                 c(1,-1,0.5,-0.5,0,0,1,-1,0.5),
                 0,1)

non2.0_sum <- summary_func(non2.0,non2.0_true)
write.csv(non2.0_sum,file="non2.4_2.0_sum.csv")

