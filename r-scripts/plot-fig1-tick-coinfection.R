######################
#  vocano plot for tick data
#  with permutation test 
######################
library(tidyverse)
library(readxl)
library(ggrepel)
setwd("../Dropbox/ag-div/")
#setwd("/Users/lai/Dropbox/ag-div/")
tick.ind <- read_xlsx("data/Sup-data-master-copy.xlsx", sheet = 1)
tick.nym <- filter(tick.ind, stage=='nymph')
#tick.adt <- filter(tick.ind, stage!='nymph')
#tick.cts <- tick.nym[,1:20]
#tick.cts <- tick.adt[,1:20]

tick.ind <- tick.nym # using only nymphs

# select 17 most common alleles
tick.cts <- tick.ind %>% select(c(2,3,5,7:20)) 

# shuffle 1/0 within by columns
tick.mix <- lapply(1:1000, function(i) {apply(tick.cts, 2, function(x) sample(x)) }) 
freq <- apply(tick.cts, 2, sum)

pair.total <- 0;
for(i in 1:(ncol(tick.cts)-1)) { # for allle i
  freq.i <- freq[i];
  for(j in (i+1):ncol(tick.cts)) { # for allele j
    freq.j <- freq[j]
    pair.total <- pair.total +  length(which(tick.cts[,i]==1 & tick.cts[,j]==1));
  }
}

out <- vector("list")
item <- 1
for(i in 1:(ncol(tick.cts)-1)) { # for allle i
  freq.i <- freq[i];
  for(j in (i+1):ncol(tick.cts)) { # for allele j
    freq.j <- freq[j]
    ob.cts <- length(which(tick.cts[,i]==1 & tick.cts[,j]==1));# num of obs pairs
    ob.cts <- ifelse(ob.cts > 0, ob.cts, 1); # avoid zero counts
    exp.cts <- freq.i /sum(freq) * freq.j / sum(freq) * pair.total;
    sim.cts <- unlist(lapply(tick.mix, function(x) length(which(x[,i]==1 & x[,j]==1)))); # num of sim pairs
    p.lower <- length(which(sim.cts<=ob.cts))/1000;
    p.lower <- ifelse(p.lower > 0, p.lower, 1/1000);
    p.upper <- length(which(sim.cts>=ob.cts))/1000;
    p.upper <- ifelse(p.upper > 0, p.upper, 1/1000);
    out[[item]] <- tibble(ale.1 = colnames(tick.cts)[i], ale.2 = colnames(tick.cts)[j], odds =log2(ob.cts)-log2(exp.cts), p.value = ifelse(ob.cts > exp.cts, p.upper, p.lower));
    item <- item + 1
  }
}

out.df <- bind_rows(out)
out.df <- mutate(out.df, pair = paste(ale.1, ale.2, sep="-"))
#write.table(out.df[,3:5], "tick-pair.tsv", quote=F, col.names = F, row.names = F, sep="\t")
out.df <- mutate(out.df, class = ifelse(p.value <= 0.01, "sig", "ns"))

ggplot(data = out.df, aes(x=odds, y=p.value, color = class)) + geom_point(shape=1) + geom_vline(xintercept = 0, linetype=2) + scale_y_log10() + geom_text_repel(data = subset(out.df, class == "sig"), aes(label = pair), size = 3) + scale_color_manual(breaks=c("sig", "ns"), values = c("black", "gray")) + geom_hline(yintercept = 0.01, linetype=2) + theme_bw() + theme(legend.position = "none")
