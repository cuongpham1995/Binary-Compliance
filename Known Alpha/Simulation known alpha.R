
#putting all functions bi.sim.all.correct, bi.sim.faz.wrong, bi.sim.fz.wrong, bi.sim.Q.wrong together to create the tables in the 
#simulation section (alpha known) in the paper
run.simulation.fix.alpha1 = function(n.cores = 8, n.sim ,  seed = seed, alpha_n1, alpha_p1, sample.size){
  set.seed(123)
  library(foreach)
  library(doParallel)
  registerDoParallel(cores=n.cores)
  
  all.correct = foreach(i = sample(1:10000,n.sim,replace = T), 
                        .combine = 'c', .inorder = F, .errorhandling = "remove") %dopar% bi.sim.correct(seed.l =  i*sample(10:100,1), seed.y = i,size = sample.size,alpha_n1 = alpha_n1, alpha_p1 = alpha_p1)
  
  all.correct1 = matrix(all.correct, ncol =18, byrow = T)
  #result.250
  result.all.correct = rbind(colMeans(all.correct1,na.rm = T), apply(all.correct1,2,sd))
  
  colnames(result.all.correct) = c("True value function", "Correct rate (prop)", 
                                   "Value Function Estimator (prop)", "Value Function MR", "Term 2.1 (prop)","Term 2.2 (prop)", "Term 4.1 (prop)", "Term 4.2 (prop)" , "Correct rate (prop)", 
                                   "Correct Rate (eric)", "Value Function Theory (eric)","Correct Rate (eric)", 
                                   "Correct Rate (OWL)", "Correct Rate (OWL)",
                                   "Correct Rate (MR)", "Value Function Theory (MR)" ,"Correct Rate (MR)", "Number of compliers")
  
  rownames(result.all.correct) = c("mean", "sd")
  print(result.all.correct)
  
  ###########################################################################
  fz.wrong = foreach(i = sample(1:10000, n.sim,replace = T), 
                     .combine = 'c', .inorder = F, .errorhandling = "remove") %dopar% bi.sim.fz.mis(seed.l =  i*sample(10:100,1), seed.y = i,size = sample.size, alpha_n1 = alpha_n1, alpha_p1 = alpha_p1)
  
  fz.wrong = matrix(fz.wrong, ncol =18, byrow = T)
  #result.250
  result.fz.wrong = rbind(colMeans(fz.wrong,na.rm = T), apply(fz.wrong,2,sd))
  
  colnames(result.fz.wrong) = c("True value function", "Correct rate (prop)", 
                                "Value Function Estimator (prop)", "Value Function MR", "Term 2.1 (prop)","Term 2.2 (prop)", "Term 4.1 (prop)", "Term 4.2 (prop)" , "Correct rate (prop)", 
                                "Correct Rate (eric)", "Value Function Theory (eric)","Correct Rate (eric)", 
                                "Correct Rate (OWL)", "Correct Rate (OWL)",
                                "Correct Rate (MR)", "Value Function Theory (MR)" ,"Correct Rate (MR)", "Number of compliers")
  rownames(result.fz.wrong) = c("mean", "sd")
  print(result.fz.wrong)
  
  #########################################################################
  faz.wrong = foreach(i = sample(1:10000,n.sim,replace = T), 
                      .combine = 'c', .inorder = F, .errorhandling = "remove") %dopar% bi.sim.faz.mis(seed.l =  i*sample(10:100,1), seed.y = i,size = sample.size,alpha_n1 = alpha_n1, alpha_p1 = alpha_p1)
  
  faz.wrong1 = matrix(faz.wrong, ncol =18, byrow = T)
  #result.250
  result.faz.wrong = rbind(colMeans(faz.wrong1,na.rm = T), apply(faz.wrong1,2,sd))
  
  colnames(result.faz.wrong) = c("True value function", "Correct rate (prop)", 
                                 "Value Function Estimator (prop)", "Value Function MR", "Term 2.1 (prop)","Term 2.2 (prop)", "Term 4.1 (prop)", "Term 4.2 (prop)" , "Correct rate (prop)", 
                                 "Correct Rate (eric)", "Value Function Theory (eric)","Correct Rate (eric)", 
                                 "Correct Rate (OWL)", "Correct Rate (OWL)",
                                 "Correct Rate (MR)", "Value Function Theory (MR)" ,"Correct Rate (MR)", "Number of compliers")
  rownames(result.faz.wrong) = c("mean", "sd")
  result.faz.wrong
  
  ########################################################################
  Q.wrong = foreach(i = sample(1:10000,n.sim,replace = T), 
                    .combine = 'c', .inorder = F, .errorhandling = "remove") %dopar% bi.sim.Q.mis(seed.l =  i*sample(10:100,1), seed.y = i,size = sample.size, alpha_n1 = alpha_n1, alpha_p1 = alpha_p1)
  
  Q.wrong1 = matrix(Q.wrong, ncol =18, byrow = T)
  #result.250
  result.Q.wrong = rbind(colMeans(Q.wrong1,na.rm = T), apply(Q.wrong1,2,sd))
  
  colnames(result.Q.wrong) = c("True value function", "Value function MC (prop)", 
                               "Value Function Estimator (prop)", "Value Function MR", "Term 2.1 (prop)","Term 2.2 (prop)", "Term 4.1 (prop)", "Term 4.2 (prop)" , "Correct rate (prop)", 
                               "Correct Rate (eric)", "Value Function Theory (eric)","Correct Rate (eric)", 
                               "Correct Rate (OWL)", "Correct Rate (OWL)",
                               "Correct Rate (MR)", "Value Function Theory (MR)" ,"Correct Rate (MR)", "Number of compliers")
  rownames(result.Q.wrong) = c("mean", "sd")
  print(result.Q.wrong)
  
  combined.result = rbind(result.all.correct, result.fz.wrong, result.faz.wrong, result.Q.wrong)
  row.names(combined.result) = c("all correct mean", "all correct sd", "fz mis mean", "fz mis sd", "faz mis mean", "faz mis sd", "Q wrong mean", "Q wrong sd")
  
  return(combined.result)
}

################################################### main simulation ################################
##### The case when (alpha_0-, alpha_0+, alpha_1-, alpha_1+) = (0,0,-0.5,0.5)

all.sim.500.n05.p05 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = -0.5, alpha_p1 = 0.5, sample.size = 500)

############# supplemtary #######################
#2000 sample size
all.sim.2000.p05.p05 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = 0.5, alpha_p1 = 0.5, sample.size = 2000)
all.sim.2000.n05.p05 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = -0.5, alpha_p1 = 0.5, sample.size = 2000)
all.sim.2000.p05.n05 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = 0.5, alpha_p1 = -0.5, sample.size = 2000)
all.sim.2000.p0.n0 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = 0, alpha_p1 = 0, sample.size = 2000)

#500 sample size
all.sim.500.p05.p05 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = 0.5, alpha_p1 = 0.5, sample.size = 500)
all.sim.500.p05.n05 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = 0.5, alpha_p1 = -0.5, sample.size = 500)
all.sim.500.p0.n0 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = 0, alpha_p1 = 0, sample.size = 500)

#200 sample size
all.sim.200.p05.p05 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = 0.5, alpha_p1 = 0.5, sample.size = 200)
all.sim.200.n05.p05 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = -0.5, alpha_p1 = 0.5, sample.size = 200)
all.sim.200.p05.n05 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = 0.5, alpha_p1 = -0.5, sample.size = 200)
all.sim.200.p0.n0 = run.simulation.fix.alpha1(n.sim = 8*64, seed = 123, alpha_n1 = 0, alpha_p1 = 0, sample.size = 200)



########################################################################################
################### Tables Maker #######################################################
########################################################################################


#alpha_n = 0.5, alpha_p = 0.5
dat1 = all.sim.2000.p05.p05

#alpha_n = -0.5, alpha_p = 0.5
dat2 = all.sim.2000.n05.p05

#alpha_n = 0.5, alpha_p = -0.5
dat3 = all.sim.2000.p05.n05

#alpha_n = 0, alpha_p = 0
dat4 = all.sim.2000.p0.n0

alpha.level = c("(0.5, 0.5)", "(-0.5, 0.5)", "(0.5, -0.5)", "(0, 0)")
options(scipen = 999)

#create tables
tab1 = data.frame(owl.vl = c(dat1[1,"Correct Rate (OWL)"], dat2[1,"Correct Rate (OWL)"], dat3[1,"Correct Rate (OWL)"], dat4[1,"Correct Rate (OWL)"]),
                  owl.vl.sd = c(dat1[2,"Correct Rate (OWL)"], dat2[2,"Correct Rate (OWL)"], dat3[2,"Correct Rate (OWL)"], dat4[2,"Correct Rate (OWL)"]),
                  eric.vl = c(dat1[1,"Correct Rate (eric)"], dat2[1,"Correct Rate (eric)"], dat3[1,"Correct Rate (eric)"], dat4[1,"Correct Rate (eric)"]),
                  eric.vl.sd = c(dat1[2,"Correct Rate (eric)"], dat2[2,"Correct Rate (eric)"], dat3[2,"Correct Rate (eric)"], dat4[2,"Correct Rate (eric)"]),
                  prop.vl = c(dat1[1,"Correct rate (prop)"], dat2[1,"Correct rate (prop)"], dat3[1,"Correct rate (prop)"], dat4[1,"Correct rate (prop)"]),
                  prop.vl.sd = c(dat1[2,"Correct rate (prop)"], dat2[2,"Correct rate (prop)"], dat3[2,"Correct rate (prop)"], dat4[2,"Correct rate (prop)"]),
                  mr.vl = c(dat1[1,"Correct Rate (MR)"], dat2[1,"Correct Rate (MR)"], dat3[1,"Correct Rate (MR)"], dat4[1,"Correct Rate (MR)"]),
                  mr.vl.sd = c(dat1[2,"Correct Rate (MR)"], dat2[2,"Correct Rate (MR)"], dat3[2,"Correct Rate (MR)"], dat4[2,"Correct Rate (MR)"]))
tab1 = round(tab1, 4)

tab1$owl.vl = paste(tab1$owl.vl, " (", tab1$owl.vl.sd, ") ", sep = "")
tab1$eric.vl = paste(tab1$eric.vl, " (", tab1$eric.vl.sd, ") ", sep = "")
tab1$prop.vl = paste(tab1$prop.vl, " (", tab1$prop.vl.sd, ") ", sep = "")
tab1$mr.vl= paste(tab1$mr.vl, " (", tab1$mr.vl.sd, ") ", sep = "")

tab1 = tab1[,c(1,3,5,7)]
tab1$alpha = alpha.level

# start table 2
tab2 = data.frame(owl.rate = c(dat1[1,"Correct Rate (OWL)"], dat2[1,"Correct Rate (OWL)"], dat3[1,"Correct Rate (OWL)"], dat4[1,"Correct Rate (OWL)"]),
                  owl.rate.sd = c(dat1[2,"Correct Rate (OWL)"], dat2[2,"Correct Rate (OWL)"], dat3[2,"Correct Rate (OWL)"], dat4[2,"Correct Rate (OWL)"]),
                  eric.rate = c(dat1[1,"Correct Rate (eric)"], dat2[1,"Correct Rate (eric)"], dat3[1,"Correct Rate (eric)"], dat4[1,"Correct Rate (eric)"]),
                  eric.rate.sd = c(dat1[2,"Correct Rate (eric)"], dat2[2,"Correct Rate (eric)"], dat3[2,"Correct Rate (eric)"], dat4[2,"Correct Rate (eric)"]),
                  prop.rate = c(dat1[1,"Correct rate (prop)"], dat2[1,"Correct rate (prop)"], dat3[1,"Correct rate (prop)"], dat4[1,"Correct rate (prop)"]),
                  prop.rate.sd = c(dat1[2,"Correct rate (prop)"], dat2[2,"Correct rate (prop)"], dat3[2,"Correct rate (prop)"], dat4[2,"Correct rate (prop)"]),
                  mr.rate = c(dat1[1,"Correct Rate (MR)"], dat2[1,"Correct Rate (MR)"], dat3[1,"Correct Rate (MR)"], dat4[1,"Correct Rate (MR)"]),
                  mr.rate.sd = c(dat1[2,"Correct Rate (MR)"], dat2[2,"Correct Rate (MR)"], dat3[2,"Correct Rate (MR)"], dat4[2,"Correct Rate (MR)"]))
tab2 = round(tab2, 4)

tab2$owl.rate = paste(tab2$owl.rate, " (", tab2$owl.rate.sd, ") ", sep = "")
tab2$eric.rate = paste(tab2$eric.rate, " (", tab2$eric.rate.sd, ") ", sep = "")
tab2$prop.rate = paste(tab2$prop.rate, " (", tab2$prop.rate.sd, ") ", sep = "")
tab2$mr.rate= paste(tab2$mr.rate, " (", tab2$mr.rate.sd, ") ", sep = "")

tab2 = tab2[,c(1,3,5,7)]
tab2$alpha = alpha.level

