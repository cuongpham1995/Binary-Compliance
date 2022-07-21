####This file contains the functions that use in the binary compliance simulation
####In line 42, we can find the function "gen.compl.data" This is the function that generate the binary compliance data
#########################################################################


expit<-function(x){
  exp(x)/(1+exp(x))
}

#this function w, generate the w(alpha1, alpha0, y)
w<-function(y, alpha0, alpha1){
  expit(alpha0 +alpha1*y)
}

w2<-function(alpha0, y, alpha1, gamma){
  mean(w(y, alpha0, alpha1)/gamma)
}


#this function generate the density function of the compliers
fIV = function(y, alpha0, alpha1, mu, sigma, gamma){
  w(y, alpha0, alpha1)*dnorm(y,mu,sigma)/gamma
}

#generate the covariate and unmeasured confounder
cov.gen = function(n = 500, seed){
  require(bridgedist)
  set.seed(seed)
  L = matrix(rep(NA, n*2), nrow = n)
  for(i in 1:2){
      L[,i] = runif(n,-1,1)
  }
  
  U = rbridge(n)/4 #unmeasured confounder
  
  dat = data.frame(L, U)
  return(dat)
}

#generate the outcome Y, the treatment Z, the compliance level A

#Inputs
#seed.l controls the randomness of the covariates
#seed.y controls the randomness of the outcomes
#sample size controls the number of sample, the default same size is set at 500
#alpha_n1 and alpha_p1 are alpha^-1 and alpha^+1 respectively

gen.compl.data = function(seed.l, seed.y, sample.size = 500, alpha_n1 = 0, alpha_p1 = 0){
  require(dplyr)
  
  n = sample.size
  
  #generate the two covariates from uniform distribution [-1,1] and unmeasured confounders with n samples
  L = cov.gen(n = n, seed = seed.l)[,1:2] 
  U = cov.gen(n = n, seed = seed.l)[,3]
  
  
  #generate the treatment Z. 
  #Z = -1 denotes taking treatment negative, Z = 1 denotes taking treatment positive 
  Z = rbinom(n,1,expit(0.3-2*L[,1]+2*L[,2])) 
  Z = ifelse(Z == 0,1,-1) 
  
  
  #generate the 6 PI strata based on L1, Z, and U
  #denominator of the probability of each PI
  denom = 1 + exp(0.5*L[,1] + 0.5*Z  - 1*U) + exp(-0.5*L[,1] + 0.5*Z + U) + exp(-0.5*L[,1] + 0.5*Z - U) + 
    exp(0.5*L[,1] + 0.5*Z - U) + exp(0.5*L[,1] + 0.5*Z + U) 
  
  #generate the probability of belonging to each PI based on the covariate L1 and unmeasured confounder U.
  vProb = cbind(exp(0.5*L[,1] + 0.5*Z + U)/denom,exp(-0.5*L[,1] + 0.5*Z + U)/denom, 1/denom ,  
                exp(-0.5*L[,1] + 0.5*Z - U)/denom, 
                exp(0.5*L[,1] + 0.5*Z - U)/denom,exp(0.5*L[,1] + 0.5*Z  - 1*U)/denom)
  
  
  mChoices = t(apply(vProb, 1, rmultinom, n = 1, size = 1)) #obtained the observed compliance value A
  
  
  dat = cbind.data.frame(PI = apply(mChoices, 1, function(x) which(x==1)), Z) 
  dat$A = rep(NA, nrow(dat)) #dat is a data frame with 2 column A and Z. 
  
  #Determine the compliance level A based on their PI strata and the treatment (Z)
  
  dat$A[dat$PI == 4] = dat$Z[dat$PI == 4]
  dat$A[dat$PI == 1] = -1
  dat$A[dat$PI == 2] = 1
  dat$A[dat$PI == 3] = 0
  dat$A[dat$PI == 5 & dat$Z == -1] = -1
  dat$A[dat$PI == 5 & dat$Z == 1] = 0 
  dat$A[dat$PI == 6 & dat$Z == -1] = 0
  dat$A[dat$PI == 6 & dat$Z == 1] = 1 
  dat$outcome = rep(NA, n)
  
  #################### Generate the density of Y ##################################
  
  n_n1 = nrow(dat[dat$Z == -1 & dat$A == -1,])
  n_p1 = nrow(dat[dat$Z == 1 & dat$A == 1,])
  
  L1_n1 = L[dat$Z == -1 & dat$A == -1,1]
  L2_n1 = L[dat$Z == -1 & dat$A == -1,2]
  
  L1_p1 = L[dat$Z == 1 & dat$A == 1, 1]
  L2_p1 = L[dat$Z == 1 & dat$A == 1, 2]
  
  #Density of the potential compliers: f_(-1) and f(1)
  set.seed(seed.y)
  y_n1 = rnorm(n_n1, 1+2*(L1_n1 + L2_n1 ), 0.5)
  y_p1 = rnorm(n_p1, 1-0*(L1_p1 + L2_p1 ), 0.5)
  
  
  n1 = nrow(dat[dat$Z == -1 & dat$A == 1,])
  L1n1 =  L[dat$Z == -1 & dat$A == 1,1]
  L2n1 =  L[dat$Z == -1 & dat$A == 1,2]
  
  y_n12 = rnorm(n1, 3 +(L1n1 + L2n1 ) , 0.5)
 
  # the density of the non-compliers
  n2 = nrow(dat[dat$Z == 1 & dat$A == -1,])
  L1n2 = L[dat$Z == 1 & dat$A == -1,1]
  L2n2 = L[dat$Z == 1 & dat$A == -1,2]
  
  y_p12 = rnorm(n2, -1+L1n2 + L2n2 , 0.5)
 
  
  dat$outcome[dat$Z == -1 & dat$A == -1] = y_n1
  dat$outcome[dat$Z == 1 & dat$A == 1] = y_p1
  dat$outcome[dat$Z == -1 & dat$A == 1] = y_n12
  dat$outcome[dat$Z == 1 & dat$A == -1] = y_p12
  
  ##############calculate alpha 0
  alpha_n1 = alpha_n1 #we set alpha_n1 (alpha of negative 1 treatment) 
  alpha_p1 = alpha_p1 #we set alpha_p1 (alpha of positive treatment) 
  temp_n1 = seq(-8, 8, 0.05)
  temp_p1 = temp_n1
  
  
  alpha0_n1 <-0#temp_n1[which.min(abs(1-sapply(temp_n1,w2, y = y_n1, alpha1 = alpha_n1, gamma = gamma_n1)))]
  alpha0_p1 <-0#temp_p1[which.min(abs(1-sapply(temp_p1,w2, y = y_p1, alpha1 = alpha_p1, gamma = gamma_p1)))]
  
  ####################
  #Create the density of outcome for the compliers
  
  dat$PI<-rep(-999,nrow(dat))
  #dat$PI[(S4_w_p1==1) | (S4_w_n1==1) ]<-4
  
  dat$PI[dat$Z==dat$A]<-4 #PI 4 indicates the potential compliers
  
  #generate gamma
  gamma_n1<-gamma_p1<-NULL
  
  for(i in 1:n_n1){
    temp_y_n1<-rnorm(5000, 1+2*(L1_n1[i] + L2_n1[i] ), 0.5)
    gamma_n1[i]<-mean(w(temp_y_n1,alpha0_n1,alpha_n1))
  }
  
  for(i in 1:n_p1){
    temp_y_p1<-rnorm(5000, 1-0*(L1_p1[i] + L2_p1[i] ),0.5)
    gamma_p1[i]<-mean(w(temp_y_p1,alpha0_p1,alpha_p1))
  }
  
  dat$gamma<-dat$gamma_n1<-dat$gamma_p1<-0
  dat$gamma_n1[dat$A == -1& dat$Z == -1]<-gamma_n1
  dat$gamma_p1[dat$A == 1& dat$Z == 1]<-gamma_p1
  
  dat$gamma<-dat$gamma_n1+dat$gamma_p1
  
  
  ##### sampling the outcome for compliers (with Z = 1 and Z = -1) ########
  #We are going to generate the sample of f^+_IV (complier taking treatment 1) and f^-1_IV (compliers taking treatment -1) by using the accept-reject method 
  # (https://en.wikipedia.org/wiki/Rejection_sampling)
  
  
  #generate the sample outcome (Y) for negative treatment for PI=4 (f^-1_IV) 
  M = 3.5
  #The mean of treatment each subject in the simulation depends on the 2 covariates. 
  L4_n1 = L[dat$PI ==4& dat$Z == -1,]
  L4_p1 = L[dat$PI ==4& dat$Z == 1,]
  
  temp.gamma_n1 = dat$gamma_n1[dat$PI ==4& dat$Z == -1]
  temp.gamma_p1 = dat$gamma_p1[dat$PI ==4& dat$Z == 1]
  
  mu_n1 = 1+2*(L4_n1[,1] + L4_n1[,2] )
  mu_p1 = 1-0*(L4_p1[,1] + L4_p1[,2] )
  
  #we first will generate the sample outcome (y) for the compliers that take negative treatment
  y_n1_IV = NULL
  for(i in 1:length(mu_n1)){
    set.seed(i)
    u = runif(8000)
    x = rnorm(8000,1.5,2)
    #process of sampling
    criteria = fIV(x, alpha0_n1, alpha_n1, mu_n1[i], 0.5, temp.gamma_n1[i])/(M*dnorm(x,1.5,2))
    chosen_sam = x[u < criteria]
    y_n1_IV[i] = mean(chosen_sam)#sample(chosen_sam,1)
  }
  
  #we now generate the sample outcome (y) for the compliers that take positive treatement
  y_p1_IV = NULL
  for (i in 1:length(mu_p1)) {
    set.seed(i)
    u = runif(8000)
    x = rnorm(8000,1.5,2)
    #process of sampling
    criteria = fIV(x, alpha0_p1, alpha_p1, mu_p1[i],0.5, temp.gamma_p1[i])/(M*dnorm(x,1.5,2))
    chosen_sam = x[u < criteria]
    y_p1_IV[i] = mean(chosen_sam)#sample(chosen_sam,1)
  }
  
  
  #############put the samples for compliers in the data set "dat"
  dat$outcomec<-NULL
  dat$outcomec[dat$PI ==4& dat$Z == -1] = y_n1_IV
  dat$outcomec[dat$PI ==4& dat$Z == 1] = y_p1_IV
  
  L5_n1 = L[(dat$PI != 4)& (dat$A == -1)& (dat$Z == -1),]
  
  L6_p1 = L[(dat$PI != 4)& (dat$A == 1)& (dat$Z == 1),]
  
  ###########generate the outcome for people with level of non-compliance 0.
  dat$outcome[dat$A == 0] = rnorm(length(dat$outcome[dat$A == 0]),dat$Z[dat$A==0]*5,0.1)
  
  if(dat$outcome %>% is.na() %>% sum()!= 0){
    warning("there are NAs in the simulated data")
  }
  
  #return the data set
  dat$L1 = L[,1]
  dat$L2 = L[,2]
  dat$U = U
  return(dat)
}

### Calculate the counterfactual outcome 
#Based on the generated data set from the gen.compl.data, we generate the potential outcome when the treatment Z is switchted.
#dat2 is the data that are generated from the gen.compl.data function

gen.switch.IV = function(seed.y, alpha_n1 = 0, alpha_p1 = 0, dat2){
  require(dplyr)
  dat = dat2[,c("PI","L1", "L2")] 
  L= as.matrix(dat[,c("L1","L2")])
  
  dat$Z = -1*dat2$Z #switch the treatment
  #Determine the compliance level A based on their PI strata and the treatment (Z)
  
  dat$A[dat$PI == 4] = dat$Z[dat$PI == 4]
  dat$A[dat$PI != 4] = sample(c(-1,0,1),length(dat$A[dat$PI != 4]),replace=TRUE)
  
  dat$outcome = rep(NA, nrow(dat))
  
  ####################generate the density of y
  
  n_n1 = nrow(dat[dat$Z == -1 & dat$A == -1,])
  n_p1 = nrow(dat[dat$Z == 1 & dat$A == 1,])
  
  L1_n1 = L[dat$Z == -1 & dat$A == -1,1]
  L2_n1 = L[dat$Z == -1 & dat$A == -1,2]
  
  L1_p1 = L[dat$Z == 1 & dat$A == 1, 1]
  L2_p1 = L[dat$Z == 1 & dat$A == 1, 2]
  
  #Density f_(-1) and f(1)
  set.seed(seed.y)
  y_n1 = rnorm(n_n1, 1+2*(L1_n1 + L2_n1 ), 0.5)
  y_p1 = rnorm(n_p1, 1-0*(L1_p1 + L2_p1 ), 0.5)
  
  
  n1 = nrow(dat[dat$Z == -1 & dat$A == 1,])
  L1n1 =  L[dat$Z == -1 & dat$A == 1,1]
  L2n1 =  L[dat$Z == -1 & dat$A == 1,2]
  
  #y_n12 = rnorm(n1, 2-(L1n1 + L2n1 ) , 0.1)
  y_n12 = rnorm(n1, 3 +(L1n1 + L2n1 ) , 0.5)
  #  y_n12 = rnorm(n1, -( 0.5) , 1)
  
  n2 = nrow(dat[dat$Z == 1 & dat$A == -1,])
  L1n2 = L[dat$Z == 1 & dat$A == -1,1]
  L2n2 = L[dat$Z == 1 & dat$A == -1,2]
  
  #y_p12 = rnorm(n2, 2+L1n2 + L2n2 , 0.1)
  y_p12 = rnorm(n2, -1+L1n2 + L2n2 , 0.5)
  #  y_p12 = rnorm(n2,  0.5, 1)
  
  dat$outcome[dat$Z == -1 & dat$A == -1] = y_n1
  dat$outcome[dat$Z == 1 & dat$A == 1] = y_p1
  dat$outcome[dat$Z == -1 & dat$A == 1] = y_n12
  dat$outcome[dat$Z == 1 & dat$A == -1] = y_p12
  
  ##############calculate alpha 0
  alpha_n1 = alpha_n1 #we set alpha_n1 (alpha of negative 1 treatment) 
  alpha_p1 = alpha_p1 #we set alpha_p1 (alpha of positive treatment) 
  
  
  #We set alpha0 to be 0
  alpha0_n1 <-0 
  alpha0_p1 <-0
  
  ####################
  
  #We generate the gamma+(X) and gamma-(X)
  gamma_n1<-gamma_p1<-NULL
  for(i in 1:n_n1){
    temp_y_n1<-rnorm(5000, 1+2*(L1_n1[i] + L2_n1[i] ), 0.5)
    gamma_n1[i]<-mean(w(temp_y_n1,alpha0_n1,alpha_n1))
  }
  
  for(i in 1:n_p1){
    temp_y_p1<-rnorm(5000, 1-0*(L1_p1[i] + L2_p1[i] ),0.5)
    gamma_p1[i]<-mean(w(temp_y_p1,alpha0_p1,alpha_p1))
  }
  
  dat$gamma<-dat$gamma_n1<-dat$gamma_p1<-0
  dat$gamma_n1[dat$A == -1& dat$Z == -1]<-gamma_n1
  dat$gamma_p1[dat$A == 1& dat$Z == 1]<-gamma_p1
  
  dat$gamma<-dat$gamma_n1+dat$gamma_p1
  
  
  ##### sampling the outcome for compliers
  #We are going to generate the sample of fmIV and fsIV by using the accept-reject method 
  
  
  #generate the sample outcome (y) for negative treatmet PI 4 (fsIV) 
  M = 3.5
  #The mean of treatment each subject in the simulation depends on the 2 covariates. 
  L4_n1 = L[dat$PI ==4& dat$Z == -1,]
  L4_p1 = L[dat$PI ==4& dat$Z == 1,]
  
  temp.gamma_n1 = dat$gamma_n1[dat$PI ==4& dat$Z == -1]
  temp.gamma_p1 = dat$gamma_p1[dat$PI ==4& dat$Z == 1]
  
  mu_n1 = 1+2*(L4_n1[,1] + L4_n1[,2] )
  mu_p1 = 1-0*(L4_p1[,1] + L4_p1[,2] )
  
  #we first will generate the sample outcome (y) for the compliers that take negative treatment
  y_n1_IV = NULL
  for(i in 1:length(mu_n1)){
    set.seed(i)
    u = runif(8000)
    x = rnorm(8000,1.5,2)
    #process of sampling
    criteria = fIV(x, alpha0_n1, alpha_n1, mu_n1[i], 0.5, temp.gamma_n1[i])/(M*dnorm(x,1.5,2))
    chosen_sam = x[u < criteria]
    y_n1_IV[i] = mean(chosen_sam)#sample(chosen_sam,1)
  }
  
  #we now generate the sample outcome (y) for the compliers that take positive treatement
  y_p1_IV = NULL
  for (i in 1:length(mu_p1)) {
    set.seed(i)
    u = runif(8000)
    x = rnorm(8000,1.5,2)
    #process of sampling
    criteria = fIV(x, alpha0_p1, alpha_p1, mu_p1[i],0.5, temp.gamma_p1[i])/(M*dnorm(x,1.5,2))
    chosen_sam = x[u < criteria]
    y_p1_IV[i] = mean(chosen_sam)#sample(chosen_sam,1)
  }
  
  
  #############put the samples for compliers in the data set
  dat$outcomec<-NULL
  dat$outcomec[dat$PI ==4& dat$Z == -1] = y_n1_IV
  dat$outcomec[dat$PI ==4& dat$Z == 1] = y_p1_IV
  
  L5_n1 = L[(dat$PI != 4)& (dat$A == -1)& (dat$Z == -1),]
  
  L6_p1 = L[(dat$PI != 4)& (dat$A == 1)& (dat$Z == 1),]
  
  ###########generate the outcome for people with level of non-compliance 0.
  dat$outcome[dat$A == 0] = rnorm(length(dat$outcome[dat$A == 0]),dat$Z[dat$A==0]*5,0.1)
  
  if(dat$outcome %>% is.na() %>% sum()!= 0){
    warning("there are NAs in the simulated data")
  }
  dat$L1 = L[,1]
  dat$L2 = L[,2]
  
  dat = subset(dat, select = -c(gamma_p1,gamma_n1, gamma))
  return(dat)
}