
#########################################################################
#################### Generative model ##################################
#########################################################################

#This part is similar to the Data Generator 2.R in the known alpha folder


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


#Unit test#
#y.test = rnorm(10)
#alpha1.test = 1
#alpha0.test = seq(-1,1,0.1)
#gamma.test = 0.2

#w(y.test, alpha0.test, alpha1.test)
#w2(alpha0.test[2], y.test, alpha1 =  alpha1.test, gamma.test)
#sapply(alpha0.test, w2, y = y.test, alpha1 = alpha1.test, gamma = gamma.test)
#alpha0.test + alpha1.test*y.test[1]

#length(alpha0.test)
########################################################
#We are going to generate data for the binary compliance problem. 
#We denote L be the set of covariates; A be the compliance level; Z be the treatment (IV), and Y be the outcome
#There will be two treatments: negative and positive
#There are three compliance level: negative (-1), zero (0), positive (1)

cov.gen = function(n = 500, seed){
  require(bridgedist)
  set.seed(seed)
  L = matrix(rep(NA, n*2), nrow = n)
  for(i in 1:2){
    #  L[,i] = 2*rbinom(n,1,0.5)-1
    L[,i] = runif(n,-1,1)
  }
  
  U = rbridge(n)/4
  
  dat = data.frame(L, U)
  return(dat)
}


gen.compl.data = function(seed.l, seed.y, sample.size = 500, alpha_n1 = 0, alpha_p1 = 0){
  require(dplyr)
  
  n = sample.size
  
  #generate the two covariates from uniform distribution [-1,1] and unmeasured confounders with n samples
  L = cov.gen(n = n, seed = seed.l)[,1:2] 
  U = cov.gen(n = n, seed = seed.l)[,3]
  
  
  #assigned the treatement Z. Z = -1 denotes taking treatment negative, Z = 1 denotes taking treatment positive 
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
  
  
  dat = cbind.data.frame(PI = apply(mChoices, 1, function(x) which(x==1)), Z) #dat is a data frame with 2 column A and Z. 
  dat$A = rep(NA, nrow(dat))
  
  # dat = cbind.data.frame(PI = sample(c(1:6),n,c(0.1,0.1,0.1,0.3,0.1,0.3),replace=TRUE), Z) #dat is a data frame with 2 column A and Z. 
  #  dat$A = rep(NA, nrow(dat))
  
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
  
  #Density f_(-1) and f(1)
  set.seed(seed.y)
  y_n1 = rnorm(n_n1, 1+2*(L1_n1 + L2_n1 ), 0.5)
  y_p1 = rnorm(n_p1, 1-0*(L1_p1 + L2_p1 ), 0.5)
  
  
  n1 = nrow(dat[dat$Z == -1 & dat$A == 1,])
  L1n1 =  L[dat$Z == -1 & dat$A == 1,1]
  L2n1 =  L[dat$Z == -1 & dat$A == 1,2]
  
  y_n12 = rnorm(n1, 3 +(L1n1 + L2n1 ) , 0.5)
  #  y_n12 = rnorm(n1, -( 0.5) , 1)
  
  n2 = nrow(dat[dat$Z == 1 & dat$A == -1,])
  L1n2 = L[dat$Z == 1 & dat$A == -1,1]
  L2n2 = L[dat$Z == 1 & dat$A == -1,2]
  
  y_p12 = rnorm(n2, -1+L1n2 + L2n2 , 0.5)
  #  y_p12 = rnorm(n2,  0.5, 1)
  
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
  
  # S4_w_p1<-S4_w_n1<-wpdatgen<-wndatgen<-rep(-999,nrow(dat))
  # S4_w_p1[dat$A == 1& dat$Z == 1]<-rbinom(n_p1,1,w(y_p1,alpha0_p1,alpha_p1))
  # S4_w_n1[dat$A == -1& dat$Z == -1]<-rbinom(n_n1,1,w(y_n1,alpha0_n1,alpha_n1))
  # 
  # wpdatgen[dat$A == 1& dat$Z == 1]<-w(y_p1,alpha0_p1,alpha_p1)
  # wndatgen[dat$A == -1& dat$Z == -1]<-w(y_n1,alpha0_n1,alpha_n1)
  
  
  dat$PI<-rep(-999,nrow(dat))
  #dat$PI[(S4_w_p1==1) | (S4_w_n1==1) ]<-4
  
  dat$PI[dat$Z==dat$A]<-4
  
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
  
  #sgamma_n1 = mean(dat$PI[dat$A == -1& dat$Z == -1]==4)  #PI_prop[4]/sum(PI_prop[1],PI_prop[4], PI_prop[5])
  #sgamma_p1 = mean(dat$PI[dat$A == 1& dat$Z == 1]==4)  #PI_prop[4]/sum(PI_prop[2], PI_prop[4], PI_prop[6])
  
  ##### sampling the outcome for compliers (with Z = 1 and Z = -1)
  #We are going to generate the sample of fmIV and fsIV by using the accept-reject method 
  # (https://en.wikipedia.org/wiki/Rejection_sampling)
  
  
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
  dat$L1 = L[,1]
  dat$L2 = L[,2]
  dat$U = U
  return(dat)
}


#Based on the generated data set from the gen.compl.data, we generate the potential outcome when the treatment Z is switchted.
#dat2 is the data that are generated from the gen.compl.data function
gen.switch.IV = function(seed.y, alpha_n1 = 0, alpha_p1 = 0, dat2){
  require(dplyr)
  dat = dat2[,c("PI","L1", "L2")] 
  L= as.matrix(dat[,c("L1","L2")])
  
  dat$Z = -1*dat2$Z
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



#########################################################################################################
###################### Simulate for a known alpha ###############################################
#####################################################################################################

#This code is similar to the Q mis 2.R in known alpha folder

bi.sim.Q.mis = function(seed.l, seed.y, true.al_n1, true.al_p1, alpha_n1 = 0, alpha_p1 = 0,alpha0_n1 = 0, alpha0_p1 = 0, size){
  tryCatch({
    require(locClass)
    require(dplyr)
    require(nnet)
    require(matrixStats)
    require(hal9001)
    #train data set
    dat1 = gen.compl.data(seed.l = seed.l, seed.y = seed.y, alpha_n1 = true.al_n1, 
                          alpha_p1 = true.al_p1,  sample.size = size)
    dat = dat1
    
    alpha0_n1 = alpha0_n1
    alpha0_p1 = alpha0_p1
    
    #calculate the weight w+ for Z = 1 and w- for Z = -1
    dat$wd = ifelse(dat$Z == -1, w(dat$outcome, alpha0_n1, alpha_n1), w(dat$outcome, alpha0_p1, alpha_p1))
    
    
    ############calculate propensity score
    #propensity score is correctly specified
    mod.logit  = glm(as.factor(Z) ~ L1 + L2, data = dat, family = binomial(link = logit))
    fz = predict(mod.logit, type = "response")
    fz = ifelse(dat$Z == 1, fz, 1 - fz)
    dat$fz = fz #propensity score
    
    
    #misspecified p(A|Z,X)
    mod.multi = multinom(A ~ Z + L1 + L2, data = dat)
    faz = predict(mod.multi, type = "prob")
    dat$faz = ifelse(dat$Z == 1,faz[,3],faz[,1])
    
    ############### Doubly Robust Estimator##############################
    dat.PI = dat #subset the train data to the data set where A = Z. Those are people with potential to be 
    #become compliers
    
    #calculate E[yw+ | A = 1, Z = 1, L1, L2]
    #Q+ and gamma+ are correctly specified
    mod.lm = lm(outcome ~ L1 + L2, data = dat[dat$Z == 1 & dat$A == 1, ])
    sum.modlm = summary(mod.lm)
    coef.p = sum.modlm$coefficients[,1]
    p.gammap = NA
    
    for(i in 1:nrow(dat.PI)){
      temp.yp = rnorm(5000, coef.p[1] , 0.5) #the true model for yp (y positive) has only intercept
      #Ey.wp[i] = mean(temp.yp*w(temp.yp,alpha0 = alpha0_p1, alpha1 = alpha_p1)) 
      p.gammap[i] = mean(w(temp.yp, alpha0 = alpha0_p1, alpha1 = alpha_p1)) #calculate gamma+ (E[w+| A = Z = 1, L])
    }
    
    mod.Qp.mis = lm(outcome ~ 1, data = dat[dat$Z == 1 & dat$A == 1, ])
    Ey.wp = predict(mod.Qp.mis, newdata = dat.PI) 
    
    dat.PI$Ey.wp = Ey.wp  #E[Yw+| A = Z = 1, L]
    dat.PI$p.gammap = p.gammap #E[w+| A = Z = 1, L]
    
    #calculate E[yw- | A = -1, Z = -1, L]
    #Q- and gamma- are corretly specified
    mod.lm2 = lm(outcome ~ L1 + L2, data = dat[dat$Z == -1 & dat$A == -1, ])
    sum.modlm2 = summary(mod.lm2)
    coef.n = sum.modlm2$coefficients[,1]
    p.gamman = NA
    
    for(i in 1:nrow(dat.PI)){
      temp.yn = rnorm(5000, coef.n[1] + coef.n[2]*dat.PI$L1[i] + coef.n[3]*dat.PI$L2[i] , 0.5)
      #Ey.wn[i] = mean(temp.yn*w(temp.yn,alpha0 = alpha0_n1, alpha1 = alpha_n1)) #calculate E[Yw-|A = Z = -1, L]
      p.gamman[i] = mean(w(temp.yn, alpha0 = alpha0_n1, alpha1 = alpha_n1)) #calculate gamma- (E[w-|A= Z = -1, L]) 
    }
    
    mod.Qn.mis = lm(outcome ~ 1, data = dat[dat$Z == -1 & dat$A == -1, ])
    Ey.wn = predict(mod.Qn.mis, newdata = dat.PI)
    
    dat.PI$Ey.wn = Ey.wn
    dat.PI$p.gamman = p.gamman
    
    #multiple robust estimator for E[Y1|PI = 4, L]
    delta.p = dat.PI$Ey.wp/dat.PI$p.gammap #calculate E[yw+| A=Z=1, L]/E[w+| A=Z=1, L] 
    term.p = delta.p + (dat.PI$A*(dat.PI$A+dat.PI$Z)*(dat.PI$A + 1))/(4*dat.PI$fz*dat.PI$faz*dat.PI$p.gammap)*(dat.PI$outcome*dat.PI$wd - dat.PI$Ey.wp - delta.p*(dat.PI$wd - dat.PI$p.gammap))
    
    #multiple robust estimator for E[Y-1|PI = 4, L]
    delta.n = dat.PI$Ey.wn/dat.PI$p.gamman
    term.n = delta.n + (dat.PI$A*(dat.PI$A+dat.PI$Z)*(1 - dat.PI$A))/(4*dat.PI$fz*dat.PI$p.gamman*dat.PI$faz)*(dat.PI$outcome*dat.PI$wd - dat.PI$Ey.wn - delta.n*(dat.PI$wd - dat.PI$p.gamman))
    
    dat.PI$w = dat.PI$Z*(term.p - term.n) #the weight for SVM is Z*delta
    dat.PI$lab = sign(dat.PI$w)*dat.PI$Z
    
    mod.mr = wsvm(as.factor(lab) ~ L1 + L2, data = dat.PI, case.weights = abs(dat.PI$w), kernel = "linear", 
                  cross = 10, scale =  F)
    
    #############################################################
    #############################################################
    #tchetgen and proportion and OWL method
    
    dat$w = (dat$A*(dat$Z+dat$A)*dat$outcome*dat$wd)/(2*dat$gamma*dat$fz*dat$faz+2*(dat$gamma==0)) # proposed weight
    dat$w2 = (dat$Z*dat$A*dat$outcome)/dat$fz #Tchetgen weight
    dat$w3 = dat$outcome/dat$fz #OWL weight
    
    dat$lab = sign(dat$w)*dat$Z
    dat$lab2 = sign(dat$w2)*dat$Z
    dat$lab3 = sign(dat$w3)*dat$Z
    
    
    #models building
    mod.prop = wsvm(as.factor(lab) ~ L1 + L2, data = dat[dat$A==dat$Z,], case.weights = abs(dat$w[dat$A==dat$Z]), kernel = "linear", 
                    cross = 10, scale =  F)
    
    mod.tchetgen = wsvm(as.factor(lab2) ~ L1 + L2, data = dat[dat$A!=0,], case.weights = abs(dat$w2[dat$A!=0]), kernel = "linear", 
                        cross = 10, scale =  F)
    
    mod.owl = wsvm(as.factor(lab3) ~ L1 + L2, data = dat, case.weights = abs(dat$w3), kernel = "linear", 
                   cross = 10, scale =  F)
    
    ###########################################################################################################
    ############################################## Testing ####################################################
    ###########################################################################################################
    #true.al_n1 = 0.5
    #true.al_p1 = -0.5
    dat1.test = gen.compl.data(seed.l = seed.l, seed.y = seed.y, alpha_n1 = true.al_n1, alpha_p1 = true.al_p1, sample.size = size)
    dat2.test = gen.switch.IV(seed.y = seed.y*3, alpha_n1 = true.al_n1, alpha_p1 = true.al_p1, dat2 = dat1.test)
    dat.test =  left_join(dat1.test, dat2.test, by = c("L1", "L2", "PI"), suffix = c("",".2"))
    
    fitted.prop = predict(mod.prop, newdata = dat.test[,c("L1", "L2")])
    fitted.prop = as.numeric(as.character(fitted.prop))
    
    fitted.owl = predict(mod.owl, newdata = dat.test[,c("L1","L2")])
    fitted.owl = as.numeric(as.character(fitted.owl))
    
    fitted.tchetgen = predict(mod.tchetgen, newdata = dat.test[,c("L1", "L2")])
    fitted.tchetgen = as.numeric(as.character(fitted.tchetgen))
    
    fitted.mr = predict(mod.mr, newdata = dat.test[,c("L1","L2")])
    fitted.mr = as.numeric(as.character(fitted.mr))
    
    #Calculate value function for the regime D: V[D(L)] = E[ I{D(L) = 1}*E[Y(1)|L, U] + I{D(L) = -1}*E[Y(-1)|L, U] ]
    value.tchetgen = ifelse(fitted.tchetgen[dat.test$PI == 4] == dat.test$Z[dat.test$PI == 4], dat.test$outcomec[dat.test$PI == 4]
                            , dat.test$outcomec.2[dat.test$PI == 4])
    value.tchetgen.mc = mean(value.tchetgen,na.rm = T)
    
    value.binary = ifelse(fitted.prop[dat.test$PI == 4] == dat.test$Z[dat.test$PI == 4], dat.test$outcomec[dat.test$PI == 4]
                          , dat.test$outcomec.2[dat.test$PI == 4])
    value.prop.mc = mean(value.binary,na.rm = T)
    
    value.func.owl = ifelse(fitted.owl[dat.test$PI == 4] == dat.test$Z[dat.test$PI == 4], dat.test$outcomec[dat.test$PI == 4]
                            , dat.test$outcomec.2[dat.test$PI == 4])
    value.owl.mc = mean(value.func.owl,na.rm = T)
    
    value.mr = ifelse(fitted.mr[dat.test$PI == 4] == dat.test$Z[dat.test$PI == 4], dat.test$outcomec[dat.test$PI == 4]
                      , dat.test$outcomec.2[dat.test$PI == 4])
    value.mr.mc = mean(value.mr, na.rm = T)
    
    
    
    #Correct classification rate
    opt.regime = ifelse(dat.test$outcomec > dat.test$outcomec.2, dat.test$Z, dat.test$Z.2)
    opt.regime = factor(opt.regime, levels = c(-1,1))
    
    tab1 = table(opt.regime[dat.test$PI == 4] , factor(fitted.prop[dat.test$PI == 4], levels = c(-1,1))) 
    correct.rate1 = (tab1[1,1] + tab1[2,2])/sum(tab1)
    
    tab.owl = table( opt.regime[dat.test$PI == 4], factor(fitted.owl[dat.test$PI == 4], levels = c(-1,1)))
    correct.rate.owl = (tab.owl[1,1] + tab.owl[2,2])/sum(tab.owl)
    
    tab2 = table(opt.regime[dat.test$PI == 4] , factor(fitted.tchetgen[dat.test$PI == 4], levels = c(-1,1)))
    correct.rate.tchetgen = (tab2[1,1] + tab2[2,2])/sum(tab2)
    
    tab3 = table(opt.regime[dat.test$PI == 4] , factor(fitted.mr[dat.test$PI == 4], levels = c(-1,1))) 
    correct.rate3 = (tab3[1,1] + tab3[2,2])/sum(tab3)
    
    true.value.func = mean(pmax(dat.test$outcomec, dat.test$outcomec.2),na.rm = T)
    
    
    ###########################Calculate mr value function##################################
    ########################################################################################
    
    test.df = dat.test
    #### HAL
    mod.logit3  = glm(as.factor(Z) ~ L1 + L2, data = test.df, family = binomial(link = logit))
    fz3 = predict(mod.logit3, type = "response")
    fz3 = ifelse(test.df$Z == 1, fz3, 1 - fz3)
    test.df$fz = fz3 #propensity score
    
    #correctly specified p(A|Z,X)
    mod.multi2 = multinom(A ~ Z + L1 + L2, data = test.df)
    #mod.multi2 = multinom(A ~ 1, data = test.df)
    faz2 = predict(mod.multi2, type = "prob")
    test.df$faz = ifelse(test.df$Z == 1,faz2[,3],faz2[,1])
    
    # pd1<-mean(test.df$A[test.df$Z == 1] == 1)
    # pd2<-mean(test.df$A[test.df$Z == -1] == -1)
    
    test.df$wd = ifelse(test.df$Z == -1, w(test.df$outcome, alpha0_n1, alpha_n1), w(test.df$outcome, alpha0_p1, alpha_p1))
    test.df$wd2 = test.df$wd
    test.df$wd2[test.df$A!= test.df$Z] = 0
    test.df$I = fitted.prop == test.df$Z
    test.df$Ie = fitted.tchetgen == dat.test$Z
    test.df$Ir = fitted.mr == dat.test$Z
    #calculate E[yw^z/p(A|Z,X)]
    
    hal.mod1 = fit_hal(X = matrix(c(rep(1, length(test.df$L1)),test.df$L1, test.df$L2, test.df$Z), ncol = 4),
                       Y = I(test.df$outcome*test.df$wd2/(test.df$faz*test.df$gamma + 2*(test.df$gamma == 0)) ), family = "gaussian", yolo = F)
    
    test.df$EAYw = predict(hal.mod1, new_data = matrix(c(rep(1, length(test.df$L1)),test.df$L1, test.df$L2, test.df$Z), ncol = 4)) 
    
    test.df$EAYwn = predict(hal.mod1, new_data = matrix(c(rep(1, length(test.df$L1)),test.df$L1, test.df$L2, rep(-1, length(test.df$L1))), ncol = 4))
    
    test.df$EAYwp = predict(hal.mod1, new_data = matrix(c(rep(1, length(test.df$L1)),test.df$L1, test.df$L2, rep(1, length(test.df$L1))), ncol = 4))
    
    
    #### get gamma value
    
    mod.lm = lm(outcome ~ L1 + L2, data = test.df[test.df$Z == 1 & test.df$A == 1, ])
    sum.modlm = summary(mod.lm)
    coef.p = sum.modlm$coefficients[,1]
    p.gammap = NA
    
    for(i in 1:nrow(test.df)){
      temp.yp = rnorm(5000, coef.p[1] , 0.5) #the true model for yp (y positive) has only intercept
      #Ey.wp[i] = mean(temp.yp*w(temp.yp,alpha0 = alpha0_p1, alpha1 = alpha_p1)) #calculate E[Yw+| A = Z = 1, L]
      p.gammap[i] = mean(w(temp.yp, alpha0 = alpha0_p1, alpha1 = alpha_p1)) #calculate gamma+ (E[w+| A = Z = 1, L])
    }
    
    mod.Qp.mis = lm(outcome ~ 1, data = dat[dat$Z == 1 & dat$A == 1, ])
    Ey.wp = predict(mod.Qp.mis, newdata = test.df) 
    
    test.df$p.gammap = p.gammap #E[w+| A = Z = 1, L]
    test.df$Ey.wp = Ey.wp #E[yw^+|A = Z = 1,L]
    test.df$deltap = test.df$Ey.wp/test.df$p.gammap
    
    mod.lm2 = lm(outcome ~ L1 + L2, data = test.df[test.df$Z == -1 & test.df$A == -1, ])
    sum.modlm2 = summary(mod.lm2)
    coef.n = sum.modlm2$coefficients[,1]
    p.gamman = NA
    
    for(i in 1:nrow(test.df)){
      temp.yn = rnorm(5000, coef.n[1] + coef.n[2]*test.df$L1[i] + coef.n[3]*test.df$L2[i] , 0.5)
      #Ey.wn[i] = mean(temp.yn*w(temp.yn,alpha0 = alpha0_n1, alpha1 = alpha_n1)) #calculate E[Yw-|A = Z = -1, L]
      p.gamman[i] = mean(w(temp.yn, alpha0 = alpha0_n1, alpha1 = alpha_n1)) #calculate gamma- (E[w-|A= Z = -1, L]) 
    }
    
    mod.Qn.mis = lm(outcome ~ 1, data = dat[dat$Z == -1 & dat$A == -1, ])
    Ey.wn = predict(mod.Qn.mis, newdata = test.df)
    
    test.df$p.gamman = p.gamman #(E[w-|A= Z = -1, L])
    test.df$Ey.wn = Ey.wn #E[Yw-|A = Z = -1, L]
    test.df$deltan = test.df$Ey.wn/test.df$p.gamman
    
    test.df$delta = ifelse(test.df$Z == 1, test.df$deltap, test.df$deltan)
    test.df$delta[test.df$PI == -999] <- 0
    
    test.df$p.gamma = ifelse(test.df$Z == 1, test.df$p.gammap, test.df$p.gamman)
    test.df$p.gamma[test.df$PI == -999] <- 0
    
    ############ putting everything together ############################  
    term1 = (test.df$I*test.df$outcome*test.df$wd2)/(test.df$p.gamma*test.df$fz*test.df$faz+2*(test.df$p.gamma==0) )
    
    term2.1 = (fitted.prop == test.df$Z)*test.df$EAYw/test.df$fz
    term2.2 = (fitted.prop == 1)*test.df$EAYwp + (fitted.prop == -1)*test.df$EAYwn
    
    
    term3.1 = (fitted.prop == 1)*test.df$deltap*(test.df$A == test.df$Z & test.df$Z == 1)*(test.df$wd2 - test.df$p.gamma)/(test.df$fz*test.df$faz)
    term3.2 = (fitted.prop == -1)*test.df$deltan*(test.df$A == test.df$Z & test.df$Z == -1)*(test.df$wd2 - test.df$p.gamma)/(test.df$fz*test.df$faz)
    
    
    term4.1 = (test.df$A == test.df$Z)*(fitted.prop == test.df$Z)*test.df$delta/(test.df$fz*test.df$faz)
    term4.2 = (test.df$Z ==1)*(fitted.prop == test.df$Z)*test.df$deltap/test.df$fz + 
      (test.df$Z == -1)*(fitted.prop == test.df$Z)*test.df$deltan/test.df$fz
    
    
    
    value.function.mr = mean(term1 - (term2.1 - term2.2) - (term3.1 + term3.2) - (term4.1 - term4.2))
    value.function.est = mean(term1)
    value.tchetgen.th = mean((test.df$Ie*test.df$outcome*test.df$wd2)/(test.df$p.gamma*test.df$fz*test.df$faz+2*(test.df$p.gamma==0) ))
    value.mr.th =  mean((test.df$Ir*test.df$outcome*test.df$wd2)/(test.df$p.gamma*test.df$fz*test.df$faz+2*(test.df$p.gamma==0) ))
    ######################################
    
    result = c(true.value.func, value.prop.mc, value.function.est,value.function.mr, mean(term2.1), mean(term2.2), mean(term4.1), mean(term4.2), 
               correct.rate1, value.tchetgen.mc, value.tchetgen.th, correct.rate.tchetgen, value.owl.mc, correct.rate.owl, 
               value.mr.mc, value.mr.th ,correct.rate3, nrow(dat.test[dat.test$PI == 4,]))
    
    names(result) = c("True value function", "Value Function MC (prop)", 
                      "Value Function Estimator (prop)", "Value Function MR", "Term 2.1 (prop)","Term 2.2 (prop)", "Term 4.1 (prop)", "Term 4.2 (prop)" , "Correct rate (prop)", 
                      "Value function MC (eric)", "Value Function Theory (eric)","Correct Rate (eric)", 
                      "Value Function MC(OWL)", "Correct Rate (OWL)",
                      "Value Function MC (MR)", "Value Function Theory (MR)" ,"Correct Rate (MR)", "Number of compliers")
    
    return(result)
  }, error = function(error_message){
    message(error_message)
    return(NA)
  })
}



###################################################################
#############  Sensitivity analysis for different alpha ###########
###################################################################

SLURM_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

seed.l = 1 + SLURM_ID
seed.y = 2 + SLURM_ID
true.al_n1 = 0.5
true.al_p1 = -0.5
size = 500


heat.df = expand.grid(alpha_n = seq(-0.5,1.5,length.out = 20), alpha_p = seq(-1.5,0.5, length.out = 20))
#heat.df = expand.grid(alpha_n = seq(0,1,1), alpha_p = seq(0,1,1))

temp.df = mapply(bi.sim.Q.mis, seed.l = seed.l, seed.y = seed.y, 
                 true.al_n1 = true.al_n1, true.al_p1 = true.al_p1, 
                 alpha_n1 = heat.df$alpha_n, alpha_p1 = heat.df$alpha_p, 
                 size = size)


heat.df1 = cbind(heat.df,t(temp.df)) #heat map for the case when everything is correctly specified

##### Create folders###############

#folder for the case where everything is correctly specified
correct <- paste(true.al_p1, true.al_n1, size, "Q mis" ,sep = " ")
output_correct <- paste0("output_",correct)
if (!file.exists(output_correct)){dir.create(output_correct)}

############## Save files ###############################
write.csv(heat.df1, file = paste(output_correct,"/result", SLURM_ID,".csv", sep = ""))

