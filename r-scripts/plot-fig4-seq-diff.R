####################################
# seq diff
####################################
library(tidyverse)
library(readxl)
library(scales)
setwd("/Users/lai/Dropbox/ag-div/")

#############
# gapless
##############

syn <- c("Root", "Consensus", "Centroid_1", "Centroid_2", "Centroid_3", "Centroid_4")
d <- read_tsv("gapless-dist.tsv", col_names = F)
d1 <- d %>% filter(!(X1 %in% syn | X2 %in% syn))

d1 %>%ggplot(aes(x=X8, fill = X9)) + geom_density(alpha=0.5) + theme_bw()+ theme(legend.position = "bottom")

d1 %>% ggplot(aes(x = X9, y = X8)) + geom_boxplot(outlier.shape = NA) + geom_jitter(shape=1) + theme_bw() 

t.test(data = d1, X8 ~ X9)

#######################
# pairwise OspC identity
#####################
x <- read_excel("data/Sup-data-master-copy.xlsx", sheet = 4)
syn <- c("Root", "Consensus", "Centroid_1", "Centroid_2", "Centroid_3", "Centroid_4")
nat <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "T", "U")

x2 <- x %>% filter(class != 'native')
out <- vector("list")
idx <- 1
for(i in 1:length(syn)) {
  ag.s <- syn[i]
  for(j in 1:length(nat)) {
    ag.n <- nat[j]
    diff <- x2 %>% filter((ag1 == ag.s & ag2 == ag.n) | (ag1 == ag.n & ag2 == ag.s)) %>% pull(frac.gapless)
    out[[idx]] <- c(syn = ag.s, nat = ag.n, diff = as.numeric(diff))
    idx <- idx + 1
  }
}
out.df <- bind_rows(out) %>% mutate(diff = as.numeric(diff))

out.mean <- out.df %>% group_by(syn) %>% summarise(diff.mean = mean(diff))

out.df %>% ggplot(aes(nat, diff, group = syn)) + geom_point(shape=1, size=1) + geom_line() + theme_bw() + facet_wrap(~syn, ncol=6) +  geom_hline(data = out.mean, aes(yintercept = diff.mean), linetype=2)

##################
# pairwise distribution
# gap diff included
################

#x3 <- x %>% filter(class != 'evol_self')
#x3 <- x3 %>% mutate(panel = if_else(class %in% c("native", "permuted"), "grp1", "grp2"))
x %>% ggplot(aes(x=frac.var, fill = class)) + geom_density(alpha=0.5) + theme_bw()  + theme(legend.position = "bottom") + geom_vline(xintercept = 0.5, linetype =2)

x.stat <- x %>% group_by(class) %>% summarise(mean=mean(frac.var), sd=sd(frac.var), median=median(frac.var), min=min(frac.var), max=max(frac.var)) %>% arrange(median)

x$class <- factor(x$class, levels = x.stat$class)
x %>% ggplot(aes(x = class, y = frac.var)) + geom_boxplot(outlier.shape = NA) + geom_jitter(shape=1, width = 0.25, color="gray") + theme_bw() 

