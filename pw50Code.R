##Code for finding PW50 values-
library(mgcv)

#Data: fromTodfList: This is a list of lists of dataframes. 
#For each allele i, there is a list of 15 dataframes, each for an allele j!=i. 
#Each dataframe in fromTodfList can be called as fromTodfList[[i]][[j]].

load("fromTodfList")
head(fromTodfList[[1]][[1]])
##method 1:  generalized additive modeling.  
pw50_values<-c()
for (i in 1:16) {
  for (j in 1:15) {
    x<-fromTodfList[[i]][[j]][,3]
    y<-fromTodfList[[i]][[j]][,4]
    odd<-data.frame(x,y)
    nl<-gam(y ~ s(x), data = odd)   #generalized linear model
    numberOfSteps<-max(x)
    newx <-data.frame(x=seq(0,numberOfSteps,0.01))   #make a set of predicted values
    fitline = predict(nl, newdata=newx)
    
    est <-data.frame(newx,fitline)             #estimate the pw50 values
    cross <-est[which.min(abs(.5-est$fitline)),]
    pw50_values<-append(pw50_values,cross$x)
  }
}

#Some values in pw50_values = 0. These are changed to 1. 
pw50_values[c(16,114,143)]<-1

#the values for pw50, contained in the vector pw50_values,can be found in the file pw50.gam.df.csv    







########
#Method 2: Binomial Generalized Linear Model.  
#We use the fitclass classification of cr.values, in the data fromTodfList, to model the data using a binomial glm. 
#fitclass classification is determined as follows- fitclass=1, if a cr.value is >.5 and 0 otherwise.

#mods-function to model fitness as a function of mutational step.

mods<-function(data) {
  glm1<-glm(fitclass~mut.steps,family=binomial,data=data)                #binomial glm
  if (pchisq(deviance(glm1),df.residual(glm1),lower.tail = F)<.05) {     #checking the model
    print("not a good fit") 
    break  
  }
  
  return(glm1)
}

#applymods-function to deal with the structure of the fromTodfList.
applymods<-function(data) {
  a<-lapply(data,mods)
  return(a)
}

#FromToglms- lapply the glms to the list of lists in fromTodfList
FromToglms<-lapply(fromTodfList,applymods)



#finding the pw50 using the ld50 values 
ldfiveoh<-function(data) {
  return(-data$coef[1]/data$coef[2])    
}

applyldfiveoh<-function(data) {
  a<-lapply(data,ldfiveoh)
  return(a)
}
#FromTold50-returns pw50 value for all 
FromTold50<-lapply(FromToglms,applyldfiveoh)

FromTold50<-lapply(FromTold50,as.numeric)
pw50glmValues<-as.numeric(unlist(FromTold50))

#Fixing the data-  there are 3 pw50 values in pw50glmValues, 16, 114, and 143, which are too high. 
#After looking at the plots, it's clear that the values should be low.  For this reason, we change
# those pw50 values to 1. 
which(pw50glmValues>100)

#plot of allele 10 paired with allele 8. Note that the values are all >.5 . This indicates that 
#the pw50 shouldn't be high, but rather it should be low, close to 0. 
#For this reason we assign the pw50 value to be 1.

plot(fromTodfList[[10]][[8]][,3], fromTodfList[[10]][[8]][,4])

plot(fromTodfList[[10]][[8]][,3], fromTodfList[[10]][[8]][,5])

#set these values equal to 1
pw50glmValues[c(16,114,143)]<-1



###compare values
hist((pw50glmValues/pw50_values))

