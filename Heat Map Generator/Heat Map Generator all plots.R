############################################################################################
############ The one that generates them all ################


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
 

#specify the path 
 path.all.correct = "C:\\Users\\cuong\\Box\\Cuong Pham-Projects\\BlueHive\\Simulation 6 (New)\\alpha0 correct\\output_-0.5 0.5 500 all correct"
 path.fz.mis = "C:\\Users\\cuong\\Box\\Cuong Pham-Projects\\BlueHive\\Simulation 6 (New)\\alpha0 correct\\output_-0.5 0.5 500 fz mis"
 path.faz.mis = "C:\\Users\\cuong\\Box\\Cuong Pham-Projects\\BlueHive\\Simulation 6 (New)\\alpha0 correct\\output_-0.5 0.5 500 faz mis"
 path.Q.mis = "C:\\Users\\cuong\\Box\\Cuong Pham-Projects\\BlueHive\\Simulation 6 (New)\\alpha0 correct\\output_-0.5 0.5 500 Q mis"


 #create the plot 
gen.heat.map = function(path1, n = 500, type.method, plot.name){
  heat.df.dat1 = get.dat(path1,n)
  heat.df.dat2 = get.dat(paste(path1,"2",sep = ""),n)
  
   if(type.method == "prop"){
     t.correct = "Correct.rate..prop."
     t.value = "Value.Function.MC..prop."
   }else if(type.method == "mr"){
     t.correct = "Correct.Rate..MR."
     t.value = "Value.Function.MC..MR."
   }else{
     warning("Only prop or mr is accepted in type.method")
  }
  
  rgn = range( c(heat.df.dat1[,t.correct], heat.df.dat2[,t.correct])) #scale of correct rate
  rgn2 = range(c(heat.df.dat1[,t.value], heat.df.dat2[,t.value])) #scale of value function

  plot11 = ggplot(data = heat.df.dat1, aes(alpha_n, alpha_p, fill = eval(parse(text = t.correct)) )) + 
    geom_tile() + scale_fill_gradient(low="white", high="dodgerblue4", limits = c(rgn[1], rgn[2]) ) + 
    labs(title =  TeX("Correct Rate ($\\alpha_0$ correctly specified)", bold =T), x = "\u03b1-", y = "\u03b1+", fill = "Correct Rate") + 
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.key.width = unit(1.5,"cm") )


  plot12 = ggplot(data = heat.df.dat2, aes(alpha_n, alpha_p, fill = eval(parse(text = t.correct)) )) + 
    geom_tile() + scale_fill_gradient(low="white", high="dodgerblue4", limits = c(rgn[1], rgn[2])) + 
    labs(title = TeX("Correct Rate ($\\alpha_0$ incorrectly specified)", bold =T) , x = "\u03b1-", y = "\u03b1+", fill = "Correct Rate") + 
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
          legend.key.width = unit(1.5, "cm")) +
    scale_x_continuous(breaks = seq(-1.5,1.5, 0.5) )
  
  plot13 = ggplot(data = heat.df.dat1, aes(alpha_n, alpha_p, fill = eval(parse(text = t.value)) )) + 
    geom_tile() + scale_fill_gradient(low="white", high="dodgerblue4",  limits = c(rgn2[1], rgn2[2])) + 
    labs(title = TeX("Value Function ($\\alpha_0$ correctly specified)", bold = T), x = "\u03b1-", y = "\u03b1+", fill = "Value Function") + 
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.key.width = unit(1.5, "cm")) #+ scale_size_continuous(range = c(1.58,1.60),breaks = c(1.581,1.59)) 
  
  
  plot14 = ggplot(data = heat.df.dat2, aes(alpha_n, alpha_p, fill = eval(parse(text = t.value)) )) + 
    geom_tile() + scale_fill_gradient(low="white", high="dodgerblue4", limits = c(rgn2[1], rgn2[2]) ) + 
    labs(title = TeX("Value Function ($\\alpha_0$ incorrectly specified)"), x = "\u03b1-", y = "\u03b1+", fill = "Value Function") + 
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.key.width = unit(1.5,"cm")) +
    scale_x_continuous(breaks = seq(-1.5,1.5, 0.5))# + scale_size_continuous(range = c(1.58,1.60),breaks = c(1.581,1.59))

  
  print(grid.arrange(plot11, plot12, plot13, plot14, ncol = 2))
  
}




#Heat map for the proposed method
gen.heat.map(path1 = path.all.correct, type.method = "prop", plot.name = "Proposed Method (All correct)")

#Heat map for the multiple robust method
gen.heat.map(path1 = path.all.correct, type.method = "mr", plot.name = "Multiply Robust (All correct)")
gen.heat.map(path1 = path.fz.mis, type.method = "mr", plot.name = "Multiply Robust (F(Z|X) mis)")
gen.heat.map(path1 = path.faz.mis, type.method = "mr", plot.name = "Multiply Robust (F(A|Z, X) mis)")
gen.heat.map(path1 = path.Q.mis, type.method = "mr", plot.name = "Multiply Robust (Q(x) mis)")
