setwd("/Users/lai/Dropbox/ag-div/")

library(tidyverse)
library(readxl)

os = read.csv("data/ospC.csv")

x =  read_xlsx("data/Sup-data-master-copy.xlsx", sheet = 5)
x = x %>% filter(OspC != "BSA" & OspC != "D17" & OspC != "V1")

# normalize z-score

x.sum =  x %>% group_by(Serum) %>% summarise(mean = mean(OD), sd = sd(OD))
y <- x %>% left_join(x.sum, "Serum") %>% mutate(odScaled = (OD-mean)/sd) %>% select(OspC,Serum,odScaled,Host) %>% rename(OD=odScaled)

y = y %>% left_join(os, "OspC")

y$allele = factor(y$allele, levels =c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","T","U","RT","CS","CT1","CT2","CT3","CT4"))
y$class = factor(y$class, levels = c("Native","Root","Consense","Centroid"))
y$Host = factor(y$Host, levels = c("P.leucopus","Human"))

y %>% ggplot(aes(allele, OD, color=class)) +
  geom_hline(yintercept = 0, linetype=2) +
  geom_boxplot(outlier.shape = NA, fill=NA) + 
  geom_jitter(shape=1, width = 0.2, size=0.9) +
  theme_bw() + xlab("Recombinant OspCs") + ylab("z-score") +
  scale_color_manual(values=c("#aaaaaa", "#E69F00", "#56B4E9", "magenta")) +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(Host~.)
#  facet_wrap(~Host, ncol=1)

# violin
y %>% ggplot(aes(allele, OD, color=class)) +
  geom_violin() +
  geom_jitter(shape=1, width = 0.2, size=0.9) +
  theme_bw() + xlab("Recombinant OspCs") + ylab("z-score") +
  geom_point(data = y.med, aes(allele, med) , size = 2, color = "#666666") +
  scale_color_manual(values=c("#aaaaaa", "#E69F00", "#56B4E9", "magenta")) +
  geom_hline(yintercept = 0, linetype=2) +
  facet_grid(Host~.)

# t-tests
y.hu <- y %>% filter(Host == 'Human')
y.hu.root <- y.hu %>% filter(class == 'Native' | class == 'Root') 
t.test(data = y.hu.root, OD ~ class)

y.hu.consense <- y.hu %>% filter(class == 'Native' | class == 'Consense') 
t.test(data = y.hu.consense, OD ~ class)

y.hu.CT1 <- y.hu %>% filter(class == 'Native' | OspC == 'A01') 
t.test(data = y.hu.CT1, OD ~ class)

y.hu.CT2 <- y.hu %>% filter(class == 'Native' | OspC == 'A02') 
t.test(data = y.hu.CT2, OD ~ class)

y.hu.CT3 <- y.hu %>% filter(class == 'Native' | OspC == 'A03') 
t.test(data = y.hu.CT3, OD ~ class)

y.hu.CT4 <- y.hu %>% filter(class == 'Native' | OspC == 'A04') 
t.test(data = y.hu.CT4, OD ~ class)

# mouse
y.leuc <- y %>% filter(Host != 'Human')
y.leuc.root <- y.leuc %>% filter(class == 'Native' | class == 'Root') 
t.test(data = y.leuc.root, OD ~ class)

y.leuc.consense <- y.leuc %>% filter(class == 'Native' | class == 'Consense') 
t.test(data = y.leuc.consense, OD ~ class)

y.leuc.CT1 <- y.leuc %>% filter(class == 'Native' | OspC == 'A01')
t.test(data = y.leuc.CT1, OD ~ class)

y.leuc.CT2 <- y.leuc %>% filter(class == 'Native' | OspC == 'A02')
t.test(data = y.leuc.CT2, OD ~ class)

y.leuc.CT3 <- y.leuc %>% filter(class == 'Native' | OspC == 'A03')
t.test(data = y.leuc.CT3, OD ~ class)

y.leuc.CT4 <- y.leuc %>% filter(class == 'Native' | OspC == 'A04')
t.test(data = y.leuc.CT4, OD ~ class)
