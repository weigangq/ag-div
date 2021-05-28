setwd("~/Dropbox/QiuDi/pic-ospC/ELISA/")

library(tidyverse)
library(plotly)
library(readxl)
library(ggrepel)

os = read.csv("../ospC.csv")

x =  read_xlsx("elisa-41hu-10mo.xlsx")
x = x %>% filter(OspC != "BSA" & OspC != "D17" & OspC != "V1")

# normalize z-score
x.sum =  x %>% group_by(Serum) %>% summarise(mean = mean(OD), sd = sd(OD))
y <- x %>% left_join(x.sum, "Serum") %>% mutate(odScaled = (OD-mean)/sd) %>% select(OspC,Serum,odScaled,Host) %>% rename(OD=odScaled)

y = y %>% left_join(os, "OspC")
y$allele = factor(y$allele, levels =c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","T","U","RT","CS","CT1","CT2","CT3","CT4"))

yy = y
y = yy %>% filter(Host=="Human")
y = yy %>% filter(Host!="Human")

al <-levels(y$allele)
ac.list <- vector("list")
for (i in seq_along(al)) {
  df <- y %>% filter(allele == al[i]) %>% arrange(OD) %>% mutate(cumOD = cumsum(OD))
  ac.list[[i]] <- df %>% mutate(order = 1:nrow(df))
}

ac.df <- bind_rows(ac.list)
ac.max <- ac.df %>% group_by(allele) %>% summarise(max = max(order)) %>% left_join(ac.df, c("allele", "max" = "order"))

ac.df.hu = ac.df
ac.max.hu = ac.max

ac.df.mo = ac.df
ac.max.mo = ac.max

ac.df = bind_rows(ac.df.hu, ac.df.mo)
ac.max = bind_rows(ac.max.hu, ac.max.mo)

ac.df$class = factor(ac.df$class, levels = c("Native","Root","Consense","Centroid"))
ac.df$Host = factor(ac.df$Host, levels = c("P. leucopus","Human"))
ac.max$Host = factor(ac.max$Host, levels = c("P. leucopus","Human"))

ac.df %>% ggplot(aes(x=order, y=cumOD, color=class, group=allele)) +
#  geom_hline(yintercept = 0) +
  geom_line() + theme_bw() +
#  geom_text_repel(data = ac.max, aes(x=max, y=cumOD, label=allele), colour="black", nudge_x=1) +
#  geom_point(shape=1, size=0.8) +
  geom_point(size=0.2) +
  scale_color_manual(values=c("#aaaaaa", "#E69F00", "#56B4E9", "magenta")) +
  #  theme(legend.position = "none") +
  labs(x = "Number of sera", y = "Cumulative z-score") +
  theme(legend.title = element_blank()) +
  facet_grid(Host~.)
