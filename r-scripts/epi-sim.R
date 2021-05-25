#########
# Epitope and antigen crossreactivity
##########
library(stringdist)
library(GA)
library(dplyr)
library(ggplot2)

##############
# Derive ubur-allele 
##############
num.alle <- 5;
num.epi <- 10; # bit length
dist.near <- 0.5;
epi.fit <- function(x) {
  x.list <- split(x, rep(1:num.alle, each=num.epi));
  min(seq_distmatrix(x.list, method = "hamming"))/num.epi
} 
epi.ga <- ga(type = "binary", fitness = epi.fit, nBits = num.epi * num.alle, popSize = 100, maxiter = 200);

evol.pop <- epi.ga@solution[1,];
near.fit <- function(x) { # maximize total number of near's
  ref.alleles <- split(evol.pop, rep(1:num.alle, each=num.epi));
  length(which(seq_dist(x, ref.alleles, method = "hamming")/num.epi > dist.near));
} 
near.ga <- ga(type = "binary", fitness = near.fit, nBits = num.epi, popSize = 100, maxiter = 200);


##############
# Derive distribution of d 
##############
out.fit <- lapply(seq(from = 10,to = 50, by = 10), function (num.alle) {
  num.epi <- 100; # bit length
  epi.fit <- function(x) { # fitness defined by the minimum pair diff
    x.list <- split(x, rep(1:num.alle, each=num.epi));
    min(seq_distmatrix(x.list, method = "hamming"))/num.epi
  } 
  epi.ga <- ga(type = "binary", fitness = epi.fit, nBits = num.epi * num.alle);
  epi.ga;
#  c(all = num.alle, epi = num.epi, min.pair.dist = epi.ga@fitnessValue)
})

dist.pop <- lapply(out.fit, function(x) {
  num.epi <-100;
  dim <- dim(x@population);
  num.alle <- dim[2]/num.epi;
  x.list <- split(x@population[50,], rep(1:num.alle, each=num.epi));
  d.mat <- seq_distmatrix(x.list, method = "hamming")/num.epi;
  d <- c(d.mat);
  data.frame(dist=d, alle = rep(num.alle, length(d)), site = rep(num.epi, length(d)))
})  

df.dist <- bind_rows(dist.pop)
ggplot(data = df.dist, aes(x=dist, color=as.factor(alle))) + geom_density() 

# distribution comparison: sim vs obs

## a random population
num.epi <- 100; # bit length
num.alle <- 50; # ospC NE
epi.ran <- sample(0:1, size = num.epi * num.alle, replace = T)
x.list <- split(epi.ran, rep(1:num.alle, each=num.epi));
d.mat <- seq_distmatrix(x.list, method = "hamming")/num.epi;
d <- c(d.mat);
d.ran <- data.frame(dist=d, alle = rep(num.alle, length(d)), site = rep(num.epi, length(d)))

## a MAD/evolved population
  num.epi <- 100; # bit length
  num.alle <- 16; # ospC NE
  epi.fit <- function(x) { # fitness defined by the minimum pair diff
    x.list <- split(x, rep(1:num.alle, each=num.epi));
    min(seq_distmatrix(x.list, method = "hamming"))/num.epi
  } 
  epi.ga <- ga(type = "binary", fitness = epi.fit, nBits = num.epi * num.alle);
  epi.ga;
  
  x.list <- split(epi.ga@population[50,], rep(1:num.alle, each=num.epi));
  d.mat <- seq_distmatrix(x.list, method = "hamming")/num.epi;
  d <- c(d.mat);
  d.sim <- data.frame(dist=d, alle = rep(num.alle, length(d)), site = rep(num.epi, length(d)))

## ospC
  setwd("../ospc-sim//")
  d.obs <- read.table("ospC-NE.pair-diff", header = F, sep="\t")
colnames(d.obs) <- c("a1", "a2", "diff", "length", "identity")  
d.obs <- mutate(d.obs, dist = diff/100)

## plot
plot(ecdf(d.sim$dist), xlim=c(0.25,0.65), las=1, cex=1.5, main="", xlab="Distance", ylab="Cumulative Proportion")
lines(ecdf(d.obs$dist), col=2, cex=1.5)
lines(ecdf(d.ran$dist), col=4, cex=1.5)
legend("bottomright", legend = c("MAD", "OspC", "random"), col=c(1,2,4), cex=1.5, pch=16)

d.sim <- mutate(d.sim, type=rep("mad", 120))
d.ran <- mutate(d.ran, type=rep("rand", 1225))
d.oc <- mutate(d.obs, type=rep("ospc", 120))

d.all <- rbind(d.sim[,c(1,4)], d.ran[,c(1,4)], d.oc[,c(6,7)])
ggplot(data=d.all, aes(x=dist, colour=type, fill=type)) + geom_density(size=1, alpha=0.1) + theme_bw() + geom_vline(xintercept = 0.5, linetype=2)

ggplot(data=d.all, aes(x=dist, color=type)) + geom_freqpoly(binwidth=0.02, size=2) + theme_bw() + geom_vline(xintercept = 0.5, linetype=2)


##############
# Derive d(min) by length for alleles
##############
out.fit <- lapply(seq(from = 5,to = 150, by = 5), function (num.epi) {
  num.alle <- 10; # num of antigens
  epi.fit <- function(x) { # fitness defined by the minimum pair diff
    x.list <- split(x, rep(1:num.alle, each=num.epi))
    #  sum(seq_distmatrix(x.list, method = "hamming"))/(num.alle * (num.alle-1) / 2)/num.epi # average pairwise diff
    min(seq_distmatrix(x.list, method = "hamming"))/num.epi # min diff
  } 
  epi.ga <- ga(type = "binary", fitness = epi.fit, nBits = num.epi * num.alle)
  c(all = num.alle, epi = num.epi, min.pair.dist = epi.ga@fitnessValue)
})

epi.mat <- t(as.data.frame(out.fit))
rownames(epi.mat) <- 1:nrow(epi.mat)

epi.df.2 <- as.data.frame(epi.mat)
epi.df.5 <- as.data.frame(epi.mat)
epi.df.10 <- as.data.frame(epi.mat)
epi.df.20 <- as.data.frame(epi.mat)
epi.df.40 <- as.data.frame(epi.mat)
epi.df.80 <- as.data.frame(epi.mat)

epi.df <- rbind(epi.df.2, epi.df.5, epi.df.10, epi.df.20, epi.df.40, epi.df.80)
library(ggplot2)
library(dplyr)
epi.df <- mutate(epi.df, allele=as.character(paste("N=",all, sep=""))) 
epi.df <- mutate(epi.df, allele=ifelse(allele=='N=2', 'N=02', ifelse(allele=='N=5', 'N=05', allele)))
ggplot(as.data.frame(epi.df), aes(x=epi, y=min.pair.dist, color=allele)) + geom_point(size=2) + geom_line() + xlab("Num Sites") + ylab("Min Hamming Dist") +
 theme(legend.position=c(0.9,0.15)) + geom_abline(slope = 0, intercept = 0.5, linetype="dashed")
