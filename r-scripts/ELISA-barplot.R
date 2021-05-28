setwd("~/Dropbox/QiuDi/pic-ospC/ELISA/")

library(tidyverse)
library(plotly)

os = read.csv("../ospC.csv")

serumCDC = read_csv("serum-CDC.csv")
serumCDC$phe = paste('IgM',ifelse(serumCDC$IgM==0, '-', '+'), 'IgG',ifelse(serumCDC$IgG==0, '-', '+'), 'Lyme', ifelse(serumCDC$Stage=='Control','-',serumCDC$Stage), sep = '')
serumCDC = serumCDC %>% select(Serum, EIA, phe)

### summary
x =  read_xlsx("elisa-41hu-10mo.xlsx")
x = x %>% filter(OspC != "BSA" & OspC != "D17" & OspC != "V1")

x = x %>% left_join(os, "OspC")
x = x %>% select(1:3, 7:9)

x.hu = x %>% filter(Host=="Human")
x.mo = x %>% filter(Host!="Human")

#x.ctl <- x %>% filter(OspC=='BSA') %>% select(Serum, OD)
#x = x %>% left_join(x.ctl, "Serum") %>% mutate(OD = OD.x - OD.y)
#x = x %>% filter(OspC!="BSA")

#x = x %>% left_join(serumCDC, "Serum")
#x <- x %>% mutate(label=paste(Serum, ': EIA=', round(EIA,2), ' ', phe, sep=''))
#x = x %>% select(-5, -9,-10)

x$allele = factor(x$allele, levels =c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","T","U","RT","CS","CT1","CT2","CT3","CT4"))
x$class = factor(x$class, levels = c("Native","Root","Consense","Centroid"))

#x$allele = factor(x$allele, levels =c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","T","U","RT","CS","CT1","CT2","CT3","CT4","VlsE"))
#x$class = factor(x$class, levels = c("Native","Root","Consense","Centroid","Vls"))

x %>% ggplot(aes(x=allele, y=OD)) + geom_bar(stat = "identity", aes(fill=class))  +
  facet_wrap(~Serum, ncol=6) +
  theme_bw() + xlab("Recombinant OspCs") + ylab("OD450") +
  scale_fill_manual(values=c("#aaaaaa", "#E69F00", "#56B4E9", "magenta")) +
  #  scale_fill_manual(values=c("#aaaaaa", "#E69F00", "#56B4E9", "magenta","darkseagreen")) +
  theme(legend.key.size = unit(0.4, "cm"), legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



### individual exp

x = read_csv("210519.csv")

#dilution = x %>% select(2,4) %>% distinct_all()

#x <- x %>% mutate(label=paste(Serum, ': 1:', Dilution, sep=''))
x = x %>% left_join(os, "OspC")
x = x %>% left_join(serumCDC, "Serum")
x <- x %>% mutate(label=paste(Serum, ': EIA=', round(EIA,2), ' ', phe, ' 1:', Dilution, sep=''))
x = x %>% select(1:3,5,6,label)

x$allele = factor(x$allele, levels =c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","T","U","RT","CS","CT1","CT2","CT3","CT4","VlsE"))
x$class = factor(x$class, levels = c("Native","Root","Consense","Centroid","Vls"))
#x$label = factor(x$label, levels = paste("T13: 1:",levels(factor(x$Dilution)), sep=""))

x %>% ggplot(aes(x=allele, y=OD)) + geom_bar(stat = "identity", aes(fill=class))  +
  facet_wrap(~label, nrow = 4) +
  theme_bw() + xlab("Recombinant OspCs") + ylab("OD450") +
  scale_fill_manual(values=c("#aaaaaa", "#E69F00", "#56B4E9", "magenta","darkseagreen")) +
  theme(legend.key.size = unit(0.4, "cm"), legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
