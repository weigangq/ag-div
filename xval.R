setwd("~/Dropbox/QiuShared/ospc-sim/")
# Bash run:
# for r in 1 2 3 4 5 6 7 8; do  # rate categorites
#  for s in 6 10 14 16 18 20 22; do # allele num for train
#    for i in {1..10}; do  # repeat
#      ./xvalidate.pl -r 2  -s $s | while read line; do echo -ne "rep$i\t$line\n">> xval.out; 
#    done; 
#   done; 
# done &
library(dplyr)
library(ggplot2)
x <- read.table("xval.out", header = F)
colnames(x) <- c("rep", "rate", "train.size", "ref", "test", "is.train", "obs", "pred", "p.near", "p.far", "p.med")
x$good <- ifelse(x$obs == x$pred, 1, 0);

x.list <- split(x, list(x$rep, x$rate, x$train.size));
x.rate <- lapply(seq_along(x.list), function(i) {
  t <- x.list[[i]];
  n.string <- names(x.list)[i];
  rate <- gsub("rep[0-9]+.(r[1-8]).[0-9]+", "\\1", n.string, perl=T);
  size <- as.numeric(gsub("rep[0-9]+.r[1-8].([0-9]+)", "\\1", n.string, perl=T));
  tb <- table(t$is.train, t$good); 
  binom.test <- binom.test(tb[1,c(2,1)]);
  binom.train <- binom.test(tb[2,c(2,1)]);
  data.frame(accuracy=c(tb[1,2]/sum(tb[1,]), tb[2,2]/sum(tb[2,])),
    lwr=c(binom.test$conf.int[1],binom.train$conf.int[1]),
    upr=c(binom.test$conf.int[2],binom.train$conf.int[2]),
    group=c("test", "train"),
    rate.cat = rep(rate, 2),
    train.size = rep(size, 2)
    ) 
  } 
)

x.df <-bind_rows(x.rate)
#library(reshape2)
#x.df <- melt(x.rate, 1)
#x.df$group <- rep(c("test", "test.lwr", "test.upr", "train", "train.lwr", "train.upr"), 50)
#x.df <- mutate(x.df, rate = paste("r", gsub("([0-9]).[0-9]+", "\\1", L1, perl=T), sep=""))
#x.df <- mutate(x.df, size = gsub("[0-9].([0-9]+)", "\\1", L1, perl=T))

#colnames(x.df) <- c("accuracy", "dummy", "group", "site.type", "train.size")
#x.df$train.size <- as.numeric(x.df$train.size)
#library(ggplot2)
#ggplot(x.df, aes(x=train.size, y=accuracy, color=group)) + geom_line(size=1) + ylim(0,1) + theme(legend.position = "bottom") + facet_wrap(~rate.cat, nrow = 2) + geom_abline(intercept = 0.5, slope = 0, linetype="dashed") + geom_errorbar(mapping=aes(ymin=lwr, ymax=upr), width=0.5, size=1) + geom_point(size=4, shape=21, fill="white") 



x.df$train.size <- as.factor(x.df$train.size)
x.df2 <- subset(x.df, rate.cat %in% c("r1", "r2", "r3", "r4"))
x.df3 <- transform(x.df2, rate.cat=factor(rate.cat,levels=c("r4","r3","r2", "r1")))
ggplot(x.df3, aes(x=train.size, y=accuracy, color=group)) +  geom_boxplot(position=position_dodge(0.8)) + geom_jitter(position=position_jitterdodge(0.8)) + facet_wrap(~rate.cat, nrow = 2) + ylim(0,1) + geom_abline(intercept = 0.333, slope = 0, linetype="dashed") + theme(legend.position = "bottom") + theme_bw()


 
