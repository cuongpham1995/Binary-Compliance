#generat the plots that are used in the simulation part for unknown alpha
library(ggplot2)
library(gridExtra)
library(latex2exp)
##################################################
#file.choose()

#read the data set
get.dat = function(path, n){
  dat1 = 0
  for(i in 1:n){
    #print(i)
    tryCatch({
      dat = read.csv(paste(path,"\\result", i,".csv",sep=""))
      dat1 = dat1 + dat
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
  heat.df.dat1 = dat1/n
  
  return(heat.df.dat1)
}


#specify the path for each cases

#path to when all nuisance parameters are correctly specified 
path.all.correct = "C:\\Users\\cpham\\Documents\\Binary Compliance Addition Sim\\output_-0.5 0.5 500 all correct"
#path.all.correct = "C:\\Users\\cpham\\Documents\\Binary Compliance Addition Sim\\output_-0.5 0.5 500 sensitivity_mod_mis"
#path.all.correctPCA = "C:\\Users\\cpham\\Documents\\Binary Compliance Addition Sim\\output_-0.5 0.5 500 sensitivity_mod_misPCA"
 # "C:\\Users\\cpham\\Documents\\Binary Compliance Addition Sim\\output_-0.5 0.5 500 all correct" #"C:\\Users\\cpham\\Box\\Cuong Pham-Projects\\BlueHive\\Simulation 6 (New)\\alpha0 correct\\output_-0.5 0.5 500 all correct"

#path to the case where f(A|Z) is misspecified
path.fz.mis = "C:\\Users\\cpham\\Documents\\Binary Compliance Addition Sim\\output_-0.5 0.5 500 fz mis"  
  
#path to the case where f(A|Z,X) is misspecified
path.faz.mis = "C:\\Users\\cpham\\Documents\\Binary Compliance Addition Sim\\output_-0.5 0.5 500 faz mis" 

#path to the case where Q(Z,X) is misspecified
path.Q.mis = "C:\\Users\\cpham\\Documents\\Binary Compliance Addition Sim\\output_-0.5 0.5 500 Q mis2" 
#path.all.wrong = "C:\\Users\\cpham\\Box\\Cuong Pham-Projects\\BlueHive\\Simulation 5\\output_-0.5 0.5 500 all correct"

#create the plot 
#gen.heat.map = function(path1, path2, n = 500, type.method, plot.name){
   path1 = path.all.correct
   n = 500
   type.method = "prop"
   
   heat.df.dat1 = get.dat(path1,n)
  heat.df.dat2 = get.dat(path2,n)
 # heat.df.dat2 = get.dat(paste(path1,"2",sep = ""),n)
  
  if(type.method == "prop"){
    t.correct = "Correct.rate..prop."
    t.value = "Value.Function.MC..prop."
  }else if(type.method == "mr"){
    t.correct = "Correct.Rate..MR."
    t.value = "Value.Function.MC..MR."
  }else{
    warning("Only prop or mr is accepted in type.method")
  }
  
#  rgn = range( c(heat.df.dat1[,t.correct], heat.df.dat2[,t.correct])) #scale of correct rate
 # rgn2 = range(c(heat.df.dat1[,t.value], heat.df.dat2[,t.value])) #scale of value function
  
  
  #heatmaps
  plot11 = ggplot(data = heat.df.dat1, aes(alpha_n, alpha_p, fill = eval(parse(text = t.correct)) )) + 
    geom_tile() + scale_fill_gradient(low="white", high="dodgerblue4", limits = c(0.9,1) ) + 
    labs(title =  TeX("Correct Classification Rate", bold = T), x = TeX("$\\alpha^Y_{-1}$"), y = TeX("$\\alpha^Y_{+1}$"), fill = "Correct Rate") + 
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.key.width = unit(1.5,"cm") ) 
  

  
 plot13 = ggplot(data = heat.df.dat1, aes(alpha_n, alpha_p, fill = eval(parse(text = t.value))/1.62 )) + 
    geom_tile() + scale_fill_gradient(low="white", high="dodgerblue4", limits = c(0,1)) + 
    labs(title = TeX("Value Function", bold = T),x = TeX("$\\alpha^Y_{-1}$"), y = TeX("$\\alpha^Y_{+1}$"), fill = "Value Function") + 
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.key.width = unit(1.5, "cm")) #+ scale_size_continuous(range = c(1.0,1.60),breaks = c(1.581,1.59)) 
  

 #contour map
 plot11.contour = ggplot(data = heat.df.dat1, aes(alpha_n, alpha_p, z = eval(parse(text = t.correct)) )) + 
   geom_contour(aes(z = eval(parse(text = t.correct))), color = "black") + 
   geom_text_contour(aes(z = eval(parse(text = t.correct)), 
                         label = round(..level.., 3)), size = 4) +  # Round to 3 decimal places
   labs(title = TeX("Correct Rate", bold = T), x = TeX("$\\alpha^Y_{-1}$"), y = TeX("$\\alpha^Y_{+1}$"), fill = "Correct Rate") + 
   theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
         legend.key.width = unit(1.5, "cm")) + theme_bw()
 
 
 plot13.contour = ggplot(data = heat.df.dat1, aes(alpha_n, alpha_p, z = eval(parse(text = t.value)) )) + 
   geom_contour(aes(z = eval(parse(text = t.value))), color = "black") + 
   geom_text_contour(aes(z = eval(parse(text = t.value)), 
                         label = round(..level.., 3)), size = 4) +  # Round to 3 decimal places
   labs(title = TeX("Value Function", bold = T), x = TeX("$\\alpha^Y_{-1}$"), y = TeX("$\\alpha^Y_{+1}$"), fill = "Value Function") + 
   theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
         legend.key.width = unit(1.5, "cm")) + theme_bw()
 
### pca
  plot11pca = ggplot(data = heat.df.dat2, aes(alpha_n, alpha_p, fill = eval(parse(text = t.correct)) )) + 
    geom_tile() + scale_fill_gradient(low="white", high="dodgerblue4", limits = c(1.5,1.62) ) + 
    labs(title =  TeX("Correct Classification Rate", bold = T), x = TeX("$\\alpha^{PCA}_{-1}$"), y = TeX("$\\alpha^{PCA}_{+1}$"), fill = "Correct Rate") + 
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.key.width = unit(1.5,"cm") )

  plot13pca = ggplot(data = heat.df.dat2, aes(alpha_n, alpha_p, fill = eval(parse(text = t.value)) )) + 
    geom_tile() + scale_fill_gradient(low="white", high="dodgerblue4") + 
    labs(title = TeX("Value Function", bold = T),x = TeX("$\\alpha^{PCA}_{-1}$"), y = TeX("$\\alpha^{PCA}_{+1}$"), fill = "Value Function") + 
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.key.width = unit(1.5, "cm")) #+ scale_size_continuous(range = c(1.0,1.60),breaks = c(1.581,1.59)) 
   
 grid.arrange(plot11, plot13, plot11pca, plot13pca, ncol = 2) #for PCA 
# grid.arrange(plot11, plot13, ncol = 2)
  
#}

grid.arrange(plot11, plot13, ncol = 2) #for heatmaps  #930 by 600 
grid.arrange(plot11.contour, plot13.contour, ncol = 2) #for heatmaps

#Heat map for the proposed method
#gen.heat.map(path1 = path.all.correct, path2 = path.all.correctPCA, type.method = "mr", plot.name = "Proposed Method (All correct)")





