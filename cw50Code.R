##Code for finding CW50 values associated with alleles. 
#cw50df: import the dataframe cw50df.
read.csv("cw50df")

library(mgcv)


#Data: cw50df: This is a dataframe that, for each allele, gives
#allele name, # of mutations, Prob of crossreactivity for mutated sequence and alleles, and fitclass.
#fitclass=1 is cross reactivity value is >.5 and is 0 otherwise.(Fitclass will be used for the binomial glm.)





######Method 1:  generalized additive modeling.  

#Split the data on the alleles
cw50List<-split(cw50df,cw50df$Allele)

cw50_values<-c()

for (i in 1:16) {
    x<-cw50List[[i]][,2]
    y<-cw50List[[i]][,3]
    odd<-data.frame(x,y)
    nl<-gam(y ~ s(x), data = odd)   #generalized linear model
    numberOfSteps<-max(x)
    newx <-data.frame(x=seq(0,numberOfSteps,0.01))   #make a set of predicted values
    fitline = predict(nl, newdata=newx)
    est <-data.frame(newx,fitline)             #estimate the pw50 values
    cross <-est[which.min(abs(.5-est$fitline)),]
    cw50_values<-append(cw50_values,cross$x)
  }

cw50_values

#cw50.gam.df.csv: the values for cw50, contained in the vector cw50_values, can be found in the file cw50.gam.df.csv   



##########

#Method 2: Binomial Generalized Linear Model.  
#We use the fitclass classification of cr.values, in the data cw50List, to model the data using a binomial glm. 
#fitclass classification is determined as follows- fitclass=1, if a cr.value is >.5 and 0 otherwise.

#mods-function to model fitness as a function of mutational step.
mods<-function(data) {
  glm1<-glm(fitclass~mutation,family=binomial,data=data)                #binomial glm
  if (pchisq(deviance(glm1),df.residual(glm1),lower.tail = F)<.05) {     #checking the model
    print("not a good fit") 
    break  
  }
  
  return(glm1)
}

#cw50glms- lapply the glms to the list binomial glms, indexed by alleles.
cw50glms<-lapply(cw50List,mods)



#finding the cw50 using the ld50 values of glms
ldfiveoh<-function(data) {
  return(-data$coef[1]/data$coef[2])    
}


#cw50.glm.values-returns cw50 values for all alleles.
cw50.glm.values<-lapply(cw50glms,ldfiveoh)
cw50.glm.values<-as.numeric(unlist(cw50.glm.values))

#cw50.glm.df:  dataframe of pairs (allele, cw50.glm.values)
cw50.glm.df<-data.frame(alleleVector,cw50.glm.values)

#cw50.glm.df.csv: The dataframe cw50.glm.df,containing all cw50 values, found using glm method, can be found in the file cw50.glm.df.csv   

#Compare Values of the from the gam and glm method, using the ratio (cw50.glm.values/cw50.gam.values)
plot(abs(cw50.glm.df[,2]-cw50.gam.method[,2]),main="Differences in CW50 values", ylab="Difference")

