library(GA)
library(stringdist)

p <- 0.333                                #probability of 1s

len.seq <- 100                            #length of sequences
num.alle <- 10                            #number of alleles
mypopSize <- 50                           #this is the default but can be changed
mymaxiter <- 500                          #the default is 100
myelitism <- max(1,round(mypopSize*0.1))  #default is 0.05 of pop size


#function to initialize the population with prob. p for 1 (and 1-p for 0) 
mypop <- function(object) {
     #ignore ga object
     #make a random matrix
     #with mypopSize instances (rows)
     #each with num.alle seqs and 1 centroid
     m <- matrix(rbinom(mypopSize*(len.seq*(num.alle+1)),1,p),mypopSize,len.seq*(num.alle+1))
     #initialize centroid to zero
     m[,(len.seq*num.alle+1):(len.seq*(num.alle+1))]=0
     m
}


# Fitness function, to maximize the minimal Hamming distance:
epi.fit <- function(x) { 
     x.list <- split(x, rep(1:(num.alle+1), each=len.seq))
     seqs <- x.list[1:num.alle]
     centroid <- x.list[num.alle+1]
     maxdist <- max(seq_dist(centroid, seqs, method = "hamming"))/len.seq
     mindist <- min(seq_distmatrix(seqs, method = "hamming"))/len.seq
     mindist/maxdist
}

# Run genetic algorithm
mad.ga <- ga(type = "binary", fitness = epi.fit, nBits = len.seq*(num.alle+1), maxiter=mymaxiter, pcrossover=0.8, pmutation=0.1, population=mypop, popSize=mypopSize, elitism=myelitism)

#get best solution
top.pop = mad.ga@solution[1,]
x.list <- split(top.pop, rep(1:(num.alle+1), each=len.seq))
seqs <- x.list[1:num.alle]
centroid <- x.list[num.alle+1]
maxdist <- max(seq_dist(centroid, seqs, method = "hamming"))/len.seq  #this is d
mindist <- min(seq_distmatrix(seqs, method = "hamming"))/len.seq      #this is D

