####################################
# seq diff
####################################
library(tidyverse)
library(readxl)
library(scales)
setwd("/Users/lai/Dropbox/ag-div/")

#######################
# pairwise OspC identity
#####################
x <- read_excel("Sup-data.xlsx", sheet = 4)
syn <- c("Root", "Consense", "Centroid_1", "Centroid_2", "Centroid_3", "Centroid_4")
nat <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "T", "U")

x2 <- x %>% filter(class != 'permuted')
out <- vector("list")
idx <- 1
for(i in 1:length(syn)) {
  ag.s <- syn[i]
  for(j in 1:length(nat)) {
    ag.n <- nat[j]
    diff <- x2 %>% filter((ag1 == ag.s & ag2 == ag.n) | (ag1 == ag.n & ag2 == ag.s)) %>% pull(diff.valid)
    out[[idx]] <- c(syn = ag.s, nat = ag.n, diff = as.numeric(diff))
    idx <- idx + 1
  }
}
out.df <- bind_rows(out) %>% mutate(diff = as.numeric(diff))

out.mean <- out.df %>% group_by(syn) %>% summarise(diff.mean = mean(diff))

out.df %>% ggplot(aes(nat, diff, group = syn)) + geom_point(shape=1, size=1) + geom_line() + theme_bw() + facet_wrap(~syn, ncol=4) +  geom_hline(data = out.mean, aes(yintercept = diff.mean), linetype=2)

##################
# pairwise distribution
################

x3 <- x %>% filter(class != 'evol_self')
x3 <- x3 %>% mutate(panel = if_else(class %in% c("native", "permuted"), "grp1", "grp2"))
x3 %>% filter(panel == 'grp1') %>% ggplot(aes(x=diff.valid, fill = class)) + geom_density(alpha=0.5) + theme_bw() + scale_fill_manual(breaks = c("native", "permuted"), values = c("black", "gray")) + theme(legend.position = "bottom")

x3 %>%ggplot(aes(x=diff.valid, fill = class)) + geom_density(alpha=0.5) + theme_bw()+ theme(legend.position = "bottom")  + facet_wrap(~panel)

x3 %>% group_by(class) %>% summarise(mean(diff.valid), sd(diff.valid))

x4 <- x3 %>% filter(panel == 'grp1')
x4 %>% ggplot(aes(x = class, y = diff.valid)) + geom_boxplot(outlier.shape = NA) + geom_jitter(shape=1) + theme_bw() 

t.test(data = x4, diff.valid ~ class)
