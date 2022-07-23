#Generate the tables in the simulation section for known alpha 
  
  n.cores = 4 #for parallel looping
  n.sim = 500 #number of simulations
  set.seed(123) #set seed
  alpha_n1 = 0.5 #set alpha^- to be 0.5
  alpha_p1 = 0.5 #set alpha^+ to be 0.5
  sample.size = 500 #set the sample size to be 500
  library(foreach)
  library(doParallel)
  registerDoParallel(cores=n.cores)
  
  #Case I: All nuisance parameters are misspecifed
  all.correct = foreach(i = sample(1:10000,n.sim,replace = T), 
                        .combine = 'c', .inorder = F, .errorhandling = "remove") %dopar% bi.sim.correct(seed.l =  i*sample(10:100,1), seed.y = i,size = sample.size,alpha_n1 = alpha_n1, alpha_p1 = alpha_p1)
  
  all.correct1 = matrix(all.correct, ncol =18, byrow = T)
  result.all.correct = rbind(colMeans(all.correct1,na.rm = T), apply(all.correct1,2,sd))
  
  colnames(result.all.correct) = c("True value function", "Value Function MC (prop)", 
                                   "Value Function Estimator (prop)", "Value Function MR", "Term 2.1 (prop)","Term 2.2 (prop)", "Term 4.1 (prop)", "Term 4.2 (prop)" , "Correct rate (prop)", 
                                   "Value function MC (eric)", "Value Function Theory (eric)","Correct Rate (eric)", 
                                   "Value Function MC(OWL)", "Correct Rate (OWL)",
                                   "Value Function MC (MR)", "Value Function Theory (MR)" ,"Correct Rate (MR)", "Number of compliers")
  
  rownames(result.all.correct) = c("mean", "sd")
  print(result.all.correct)
  
  ########################################################################################
  #f(Z|X) and f(A|Z,X) are misspecified
  all.wrong = foreach(i = sample(1:10000,n.sim,replace = T), 
                        .combine = 'c', .inorder = F, .errorhandling = "remove") %dopar% bi.sim.wrong(seed.l =  i*sample(10:100,1), seed.y = i,size = sample.size,alpha_n1 = alpha_n1, alpha_p1 = alpha_p1)
  
  all.wrong1 = matrix(all.wrong, ncol =18, byrow = T)
  #result.250
  result.all.wrong = rbind(colMeans(all.wrong1,na.rm = T), apply(all.wrong1,2,sd))
  
  colnames(result.all.wrong) = c("True value function", "Value Function MC (prop)", 
                                   "Value Function Estimator (prop)", "Value Function MR", "Term 2.1 (prop)","Term 2.2 (prop)", "Term 4.1 (prop)", "Term 4.2 (prop)" , "Correct rate (prop)", 
                                   "Value function MC (eric)", "Value Function Theory (eric)","Correct Rate (eric)", 
                                   "Value Function MC(OWL)", "Correct Rate (OWL)",
                                   "Value Function MC (MR)", "Value Function Theory (MR)" ,"Correct Rate (MR)", "Number of compliers")
  
  rownames(result.all.wrong) = c("mean", "sd")
  print(result.all.wrong)
  
  #####################################################################################
  #Q(Z,X) is mispecified 
  Q.mis = foreach(i = sample(1:10000,n.sim,replace = T), 
                      .combine = 'c', .inorder = F, .errorhandling = "remove") %dopar% bi.sim.Q.mis(seed.l =  i*sample(10:100,1), seed.y = i,size = sample.size,alpha_n1 = alpha_n1, alpha_p1 = alpha_p1)
  
  Q.mis1 = matrix(Q.mis, ncol =18, byrow = T)
  result.Q.mis = rbind(colMeans(Q.mis1,na.rm = T), apply(Q.mis1,2,sd))
  
  colnames(result.Q.mis) = c("True value function", "Value Function MC (prop)", 
                                 "Value Function Estimator (prop)", "Value Function MR", "Term 2.1 (prop)","Term 2.2 (prop)", "Term 4.1 (prop)", "Term 4.2 (prop)" , "Correct rate (prop)", 
                                 "Value function MC (eric)", "Value Function Theory (eric)","Correct Rate (eric)", 
                                 "Value Function MC(OWL)", "Correct Rate (OWL)",
                                 "Value Function MC (MR)", "Value Function Theory (MR)" ,"Correct Rate (MR)", "Number of compliers")
  
  rownames(result.Q.mis) = c("mean", "sd")
  print(result.Q.mis)
  
  ######################################################################################
  #f(Z|X) is misspecied only
  fz.mis = foreach(i = sample(1:10000,n.sim,replace = T), 
                  .combine = 'c', .inorder = F, .errorhandling = "remove") %dopar% bi.sim.fz.mis(seed.l =  i*sample(10:100,1), seed.y = i,size = sample.size,alpha_n1 = alpha_n1, alpha_p1 = alpha_p1)
  
  fz.mis1 = matrix(fz.mis, ncol =18, byrow = T)
  result.fz.mis = rbind(colMeans(fz.mis1,na.rm = T), apply(fz.mis1,2,sd))
  
  colnames(result.fz.mis) = c("True value function", "Value Function MC (prop)", 
                             "Value Function Estimator (prop)", "Value Function MR", "Term 2.1 (prop)","Term 2.2 (prop)", "Term 4.1 (prop)", "Term 4.2 (prop)" , "Correct rate (prop)", 
                             "Value function MC (eric)", "Value Function Theory (eric)","Correct Rate (eric)", 
                             "Value Function MC(OWL)", "Correct Rate (OWL)",
                             "Value Function MC (MR)", "Value Function Theory (MR)" ,"Correct Rate (MR)", "Number of compliers")
  
  rownames(result.fz.mis) = c("mean", "sd")
  print(result.fz.mis)
  #################################################################################
  #f(A|Z,X) is misspecified
  faz.mis = foreach(i = sample(1:10000,n.sim,replace = T), 
                   .combine = 'c', .inorder = F, .errorhandling = "remove") %dopar% bi.sim.faz.mis(seed.l =  i*sample(10:100,1), seed.y = i,size = sample.size,alpha_n1 = alpha_n1, alpha_p1 = alpha_p1)
  
  faz.mis1 = matrix(faz.mis, ncol =18, byrow = T)
  result.faz.mis = rbind(colMeans(faz.mis1,na.rm = T), apply(faz.mis1,2,sd))
  
  colnames(result.faz.mis) = c("True value function", "Value Function MC (prop)", 
                              "Value Function Estimator (prop)", "Value Function MR", "Term 2.1 (prop)","Term 2.2 (prop)", "Term 4.1 (prop)", "Term 4.2 (prop)" , "Correct rate (prop)", 
                              "Value function MC (eric)", "Value Function Theory (eric)","Correct Rate (eric)", 
                              "Value Function MC(OWL)", "Correct Rate (OWL)",
                              "Value Function MC (MR)", "Value Function Theory (MR)" ,"Correct Rate (MR)", "Number of compliers")
  
  rownames(result.faz.mis) = c("mean", "sd")
  print(result.faz.mis)
  
  