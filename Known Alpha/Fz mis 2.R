###########################################################################################
########### This is the case where we misspecified f(Z|X) #################################
###########################################################################################


#Inputs: seed.l controls the randomness of the covariate L
#        seed.y controls the randomness of the outcome Y
#        alpha_n1 is the sensitivity parameter for alpha^- 
#        alpha_p1 is the sensitivity parameter for alpha^+ 
#        alpha0_n1 and alpha0_p1 is always set at 0
#        size controls the sample size of the simulation

bi.sim.fz.mis = function(seed.l, seed.y, alpha_n1 = 0, alpha_p1 = 0,alpha0_n1 = 0, alpha0_p1 = 0, size){
  tryCatch({
    require(locClass)
    require(dplyr)
    require(DynTxRegime)
    require(matrixStats)
    require(hal9001)
    require(nnet)
   
     ###################################################################################
    ######################### Training Model ##########################################
    ###################################################################################
    
    #Part I: create data set  
    #Create train data set
    
    dat1 = gen.compl.data(seed.l = seed.l, seed.y = seed.y, alpha_n1 = alpha_n1, 
                          alpha_p1 = alpha_p1,  sample.size = size)
    dat = dat1
    
    alpha0_n1 = alpha0_n1
    alpha0_p1 = alpha0_p1
    
    #calculate the weight w+ for Z = 1 and w- for Z = -1
    dat$wd = ifelse(dat$Z == -1, w(dat$outcome, alpha0_n1, alpha_n1), w(dat$outcome, alpha0_p1, alpha_p1))
    
    
    #Part II: Estimate nuisance parameters
    #f(A|Z,X) parameter is not correctly specified.
    
    ## calculate the propensity score f(Z|X)
    #propensity score is incorrectly specified
    
    mod.logit  = glm(as.factor(Z) ~ 1, data = dat, family = binomial(link = logit))
    fz = predict(mod.logit, type = "response")
    fz = ifelse(dat$Z == 1, fz, 1 - fz)
    dat$fz = fz #propensity score
    
    #specify correctly f(A|Z,X)
    mod.multi = multinom(A ~ Z + L1 + L2, data = dat)
    faz = predict(mod.multi, type = "prob")
    dat$faz = ifelse(dat$Z == 1,faz[,3],faz[,1])
    
    #Part III: Fitting model
    
    #Implement Doubly Robust Estimator (the MR method)
    dat.PI = dat #make a copy of dat
    
    #calculate E[yw+ | A = 1, Z = 1, L1, L2]
    #Q+ and gamma+ are correctly specified
    mod.lm = lm(outcome ~ L1 + L2, data = dat[dat$Z == 1 & dat$A == 1, ])
    sum.modlm = summary(mod.lm)
    coef.p = sum.modlm$coefficients[,1]
    Ey.wp = p.gammap = NA
    
    for(i in 1:nrow(dat.PI)){
      temp.yp = rnorm(5000, coef.p[1] , 0.5) #the true model for yp (y positive) has only intercept
      Ey.wp[i] = mean(temp.yp*w(temp.yp,alpha0 = alpha0_p1, alpha1 = alpha_p1)) #calculate E[Yw+| A = Z = 1, L]
      p.gammap[i] = mean(w(temp.yp, alpha0 = alpha0_p1, alpha1 = alpha_p1)) #calculate gamma+ (E[w+| A = Z = 1, L])
    }
    
    dat.PI$Ey.wp = Ey.wp  #E[Yw+| A = Z = 1, L]
    dat.PI$p.gammap = p.gammap #E[w+| A = Z = 1, L]
    
    #calculate E[yw- | A = -1, Z = -1, L]
    #Q- and gamma- are corretly specified
    mod.lm2 = lm(outcome ~ L1 + L2, data = dat[dat$Z == -1 & dat$A == -1, ])
    sum.modlm2 = summary(mod.lm2)
    coef.n = sum.modlm2$coefficients[,1]
    Ey.wn = p.gamman = NA
    
    for(i in 1:nrow(dat.PI)){
      temp.yn = rnorm(5000, coef.n[1] + coef.n[2]*dat.PI$L1[i] + coef.n[3]*dat.PI$L2[i] , 0.5)
      Ey.wn[i] = mean(temp.yn*w(temp.yn,alpha0 = alpha0_n1, alpha1 = alpha_n1)) #calculate E[Yw-|A = Z = -1, L]
      p.gamman[i] = mean(w(temp.yn, alpha0 = alpha0_n1, alpha1 = alpha_n1)) #calculate gamma- (E[w-|A= Z = -1, L]) 
    }
    
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
    
    #fitting the multiply robust model 
    mod.mr = wsvm(as.factor(lab) ~ L1 + L2, data = dat.PI, case.weights = abs(dat.PI$w), kernel = "linear", 
                  cross = 10, scale =  F)
    
    
    #Implement the IPW, IVT, and OWL method #
    dat$w = (dat$A*(dat$Z+dat$A)*dat$outcome*dat$wd)/(2*dat$gamma*dat$fz*dat$faz+2*(dat$gamma==0)) # IPW weight
    dat$w2 = (dat$Z*dat$A*dat$outcome)/dat$fz #IVT weight
    dat$w3 = dat$outcome/dat$fz #OWL weight
    
    dat$lab = sign(dat$w)*dat$Z
    dat$lab2 = sign(dat$w2)*dat$Z
    dat$lab3 = sign(dat$w3)*dat$Z
    
    #models building
    mod.prop = wsvm(as.factor(lab) ~ L1 + L2, data = dat[dat$A==dat$Z,], case.weights = abs(dat$w[dat$A==dat$Z]), kernel = "linear", 
                    cross = 10, scale =  F) #fitting IPW method
    
    mod.tchetgen = wsvm(as.factor(lab2) ~ L1 + L2, data = dat[dat$A!=0,], case.weights = abs(dat$w2[dat$A!=0]), kernel = "linear", 
                        cross = 10, scale =  F) #fitting IVT method
    
    mod.owl = wsvm(as.factor(lab3) ~ L1 + L2, data = dat, case.weights = abs(dat$w3), kernel = "linear", 
                   cross = 10, scale =  F) #fitting OWL method
    
    #####################################################################################
    ########################### Models Evaluation #######################################
    #####################################################################################
    
    
    # Part I: Obtained the test data set 
    #note that we set the seed.l differs from the one that we use to generate the train set. This will make sure that the model is not 
    #over fitting.
    
    dat1.test = gen.compl.data(seed.l = seed.l*2, seed.y = seed.y, alpha_n1 = alpha_n1, alpha_p1 = alpha_p1, sample.size = size)
    dat2.test = gen.switch.IV(seed.y = seed.y*3, alpha_n1 = alpha_n1, alpha_p1 = alpha_p1, dat2 = dat1.test)
    dat.test =  left_join(dat1.test, dat2.test, by = c("L1", "L2", "PI"), suffix = c("",".2"))
    
    #we get the predicted values for the proposed method (IPW)
    fitted.prop = predict.wsvm(mod.prop, newdata = dat.test[,c("L1", "L2")])
    fitted.prop = as.numeric(as.character(fitted.prop))
    
    #we get the predicted values for the OWL method
    fitted.owl = predict.wsvm(mod.owl, newdata = dat.test[,c("L1","L2")])
    fitted.owl = as.numeric(as.character(fitted.owl))
    
    #we get the predicted values for IVT method
    fitted.tchetgen = predict.wsvm(mod.tchetgen, newdata = dat.test[,c("L1", "L2")])
    fitted.tchetgen = as.numeric(as.character(fitted.tchetgen))
    
    #we obtain the predicted values using the multiply robust (MR) method
    fitted.mr = predict.wsvm(mod.mr, newdata = dat.test[,c("L1","L2")])
    fitted.mr = as.numeric(as.character(fitted.mr))
    
    ###############
    # Part III: Calculate the empirical value function for the regime D: V[D(L)] = E[ I{D(L) = 1}*E[Y(1)|L, U] + I{D(L) = -1}*E[Y(-1)|L, U] ]
    
    #the value function using IVT method
    value.tchetgen = ifelse(fitted.tchetgen[dat.test$PI == 4] == dat.test$Z[dat.test$PI == 4], dat.test$outcomec[dat.test$PI == 4]
                            , dat.test$outcomec.2[dat.test$PI == 4])
    value.tchetgen.mc = mean(value.tchetgen,na.rm = T)
    
    #IPW method (proposed method)
    value.binary = ifelse(fitted.prop[dat.test$PI == 4] == dat.test$Z[dat.test$PI == 4], dat.test$outcomec[dat.test$PI == 4]
                          , dat.test$outcomec.2[dat.test$PI == 4])
    value.prop.mc = mean(value.binary,na.rm = T)
    
    #the OWL method
    value.func.owl = ifelse(fitted.owl[dat.test$PI == 4] == dat.test$Z[dat.test$PI == 4], dat.test$outcomec[dat.test$PI == 4]
                            , dat.test$outcomec.2[dat.test$PI == 4])
    value.owl.mc = mean(value.func.owl,na.rm = T)
    
    #the multiply robust (MR) estimator
    value.mr = ifelse(fitted.mr[dat.test$PI == 4] == dat.test$Z[dat.test$PI == 4], dat.test$outcomec[dat.test$PI == 4]
                      , dat.test$outcomec.2[dat.test$PI == 4])
    value.mr.mc = mean(value.mr, na.rm = T)
    
    
    
    # Part IV: Estimate the correct classification rate
    
    #determine the true correct rate
    opt.regime = ifelse(dat.test$outcomec > dat.test$outcomec.2, dat.test$Z, dat.test$Z.2)
    opt.regime = factor(opt.regime, levels = c(-1,1))
    
    #determine the true correct rate for IPW method
    tab1 = table(opt.regime[dat.test$PI == 4] , factor(fitted.prop[dat.test$PI == 4], levels = c(-1,1))) 
    correct.rate1 = (tab1[1,1] + tab1[2,2])/sum(tab1)
    
    #determine the true correct rate for OWL method
    tab.owl = table( opt.regime[dat.test$PI == 4], factor(fitted.owl[dat.test$PI == 4], levels = c(-1,1)))
    correct.rate.owl = (tab.owl[1,1] + tab.owl[2,2])/sum(tab.owl)
    
    #determine the true correct rate for IVT method
    tab2 = table(opt.regime[dat.test$PI == 4] , factor(fitted.tchetgen[dat.test$PI == 4], levels = c(-1,1)))
    correct.rate.tchetgen = (tab2[1,1] + tab2[2,2])/sum(tab2)
    
    #determine the true correct rate for MR method
    tab3 = table(opt.regime[dat.test$PI == 4] , factor(fitted.mr[dat.test$PI == 4], levels = c(-1,1))) 
    correct.rate3 = (tab3[1,1] + tab3[2,2])/sum(tab3)
    
    ######### The maximal value function ###################### 
    true.value.func = mean(pmax(dat.test$outcomec, dat.test$outcomec.2),na.rm = T)
    
    # Part V: Estimate the value functions
    
    #calculate the value function using the multiple robust method
    test.df = dat.test
    
    
    #misspecified the f(Z|X)
    mod.logit3  = glm(as.factor(Z) ~ 1, data = test.df, family = binomial(link = logit))
    fz3 = predict(mod.logit3, type = "response")
    fz3 = ifelse(test.df$Z == 1, fz3, 1 - fz3)
    test.df$fz = fz3 #propensity score
    
    #Correctly specified f(A|Z,X)
    mod.multi2 = multinom(A ~ Z + L1 + L2, data = test.df)
    faz2 = predict(mod.multi2, type = "prob")
    test.df$faz = ifelse(test.df$Z == 1,faz2[,3],faz2[,1])
    
    #calculate the weight w+(Y), w-(Y) and the indicator I{pi(X) = Z}
    test.df$wd = ifelse(test.df$Z == -1, w(test.df$outcome, alpha0_n1, alpha_n1), w(test.df$outcome, alpha0_p1, alpha_p1))
    test.df$wd2 = test.df$wd
    test.df$wd2[test.df$A!= test.df$Z] = 0
    test.df$I = fitted.prop == test.df$Z
    test.df$Ie = fitted.tchetgen == dat.test$Z
    test.df$Ir = fitted.mr == dat.test$Z
    
    #calculate kappa(Z,X) = E[yw^z/p(A|Z,X)] using the HAL model
    hal.mod1 = fit_hal(X = matrix(c(rep(1, length(test.df$L1)),test.df$L1, test.df$L2, test.df$Z), ncol = 4),
                       Y = I(test.df$outcome*test.df$wd2/(test.df$faz*test.df$gamma + 2*(test.df$gamma == 0)) ), family = "gaussian", yolo = F)
    
    test.df$EAYw = predict(hal.mod1, new_data = matrix(c(rep(1, length(test.df$L1)),test.df$L1, test.df$L2, test.df$Z), ncol = 4)) 
    
    test.df$EAYwn = predict(hal.mod1, new_data = matrix(c(rep(1, length(test.df$L1)),test.df$L1, test.df$L2, rep(-1, length(test.df$L1))), ncol = 4))
    
    test.df$EAYwp = predict(hal.mod1, new_data = matrix(c(rep(1, length(test.df$L1)),test.df$L1, test.df$L2, rep(1, length(test.df$L1))), ncol = 4))
    
    #### get gamma value
    
    mod.lm = lm(outcome ~ L1 + L2, data = test.df[test.df$Z == 1 & test.df$A == 1, ])
    sum.modlm = summary(mod.lm)
    coef.p = sum.modlm$coefficients[,1]
    Ey.wp = p.gammap = NA
    
    for(i in 1:nrow(test.df)){
      temp.yp = rnorm(5000, coef.p[1] , 0.5) #the true model for yp (y positive) has only intercept
      Ey.wp[i] = mean(temp.yp*w(temp.yp,alpha0 = alpha0_p1, alpha1 = alpha_p1)) #calculate E[Yw+| A = Z = 1, L]
      p.gammap[i] = mean(w(temp.yp, alpha0 = alpha0_p1, alpha1 = alpha_p1)) #calculate gamma+ (E[w+| A = Z = 1, L])
    }
    
    test.df$p.gammap = p.gammap #E[w+| A = Z = 1, L]
    test.df$Ey.wp = Ey.wp #E[yw^+|A = Z = 1,L]
    test.df$deltap = test.df$Ey.wp/test.df$p.gammap
    
    mod.lm2 = lm(outcome ~ L1 + L2, data = test.df[test.df$Z == -1 & test.df$A == -1, ])
    sum.modlm2 = summary(mod.lm2)
    coef.n = sum.modlm2$coefficients[,1]
    Ey.wn = p.gamman = NA
    
    for(i in 1:nrow(test.df)){
      temp.yn = rnorm(5000, coef.n[1] + coef.n[2]*test.df$L1[i] + coef.n[3]*test.df$L2[i] , 0.5)
      Ey.wn[i] = mean(temp.yn*w(temp.yn,alpha0 = alpha0_n1, alpha1 = alpha_n1)) #calculate E[Yw-|A = Z = -1, L]
      p.gamman[i] = mean(w(temp.yn, alpha0 = alpha0_n1, alpha1 = alpha_n1)) #calculate gamma- (E[w-|A= Z = -1, L]) 
    }
    
    
    test.df$p.gamman = p.gamman #(E[w-|A= Z = -1, L])
    test.df$Ey.wn = Ey.wn #E[Yw-|A = Z = -1, L]
    test.df$deltan = test.df$Ey.wn/test.df$p.gamman
    
    test.df$delta = ifelse(test.df$Z == 1, test.df$deltap, test.df$deltan)
    test.df$delta[test.df$PI == -999] <- 0
    
    test.df$p.gamma = ifelse(test.df$Z == 1, test.df$p.gammap, test.df$p.gamman)
    test.df$p.gamma[test.df$PI == -999] <- 0
    
    ### putting everything together ####  
    term1 = (test.df$I*test.df$outcome*test.df$wd2)/(test.df$p.gamma*test.df$fz*test.df$faz+2*(test.df$p.gamma==0) )
    
    term2.1 = (fitted.prop == test.df$Z)*test.df$EAYw/test.df$fz
    term2.2 = (fitted.prop == 1)*test.df$EAYwp + (fitted.prop == -1)*test.df$EAYwn
    
    
    term3.1 = (fitted.prop == 1)*test.df$deltap*(test.df$A == test.df$Z & test.df$Z == 1)*(test.df$wd2 - test.df$p.gamma)/(test.df$fz*test.df$faz)
    term3.2 = (fitted.prop == -1)*test.df$deltan*(test.df$A == test.df$Z & test.df$Z == -1)*(test.df$wd2 - test.df$p.gamma)/(test.df$fz*test.df$faz)
    
    
    term4.1 = (test.df$A == test.df$Z)*(fitted.prop == test.df$Z)*test.df$delta/(test.df$fz*test.df$faz)
    term4.2 = (test.df$Z ==1)*(fitted.prop == test.df$Z)*test.df$deltap/test.df$fz + 
      (test.df$Z == -1)*(fitted.prop == test.df$Z)*test.df$deltan/test.df$fz
    
    
    #value function
    
    # MR estimator of value function
    value.function.mr = mean(term1 - (term2.1 - term2.2) - (term3.1 + term3.2) - (term4.1 - term4.2))
    
    ## IPW estimator of value function
    value.function.est = mean(term1)
    
    #calculate the value function for the method proposed in Eric's using regular value function formula (non-robust)
    value.tchetgen.th = mean((test.df$Ie*test.df$outcome*test.df$wd2)/(test.df$p.gamma*test.df$fz*test.df$faz+2*(test.df$p.gamma==0) ))
    
    #calculate the value function for the method proposed in Eric's using the robust value function formula
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

bi.sim.fz.mis(83,267,alpha_n1 = 0, alpha_p1 = 0, size = 2000)
