setwd("/Users/Weigang/Dropbox/ag-div")

library(tidyverse)
library(readxl)
library(pheatmap)
library(gplots)
library(ggrepel)
library(broom)

# simple linear regression
x <- read_excel("pred.xlsx")
colnames(x) = c("win", "test", "ols", "ridge", "lasso", "enet", "lars", "ag", "ab")
x <- x %>% mutate(od.obs = exp(-ridge), od.pred = exp(-test))
x %>% ggplot(aes(od.obs, od.pred, label=ag)) + geom_point(size=3, shape=1) + theme_bw() + geom_smooth(method = "lm") + geom_abline(intercept = 0, slope = 1, linetype=2) + geom_text_repel() + scale_x_log10() + scale_y_log10() + facet_wrap(~ab)

x %>% group_by(ab) %>% do(tidy(lm(ridge ~ test + 0, data = .)))
x %>% filter(ab == 'anti-CT4-1') %>% lm(od.pred ~ od.obs, data = .) %>% summary()


pos <- read_csv("paramter.csv")
pos %>% ggplot(aes(win, beta)) + geom_point(shape=1) + theme_bw() + geom_col()

# pair-diff
x <- read_tsv("pair-diff")
x.long <- x %>% pivot_longer(8:12, names_to = "class", values_to = "diff")
x.win <- x.long %>% group_by(window, class) %>% summarise(mean.diff = mean(diff))

x.range <- x.win %>% group_by(class) %>% summarise(max.diff=max(mean.diff))

x.win2 <- x.win %>% left_join(x.range, "class") %>% mutate(norm.diff = mean.diff/max.diff)

x.win2 %>% ggplot(aes(x=window, y=norm.diff, group=class)) + geom_point(shape=1) + geom_line() + theme_bw() + geom_hline(yintercept = c(0, 0.1, 0.2, 0.3), linetype = 2) + facet_wrap(~class)
