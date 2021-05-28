setwd("~/Dropbox/QiuDi/pic-ospC/ELISA/")

library(tidyverse)
library(plotly)
library(readxl)

#library(heatmaply)
library(pheatmap)
library(gplots)
library(scales)

os = read.csv("../ospC.csv")

x =  read_xlsx("elisa-41hu-10mo.xlsx")
x = x %>% filter(OspC != "BSA" & OspC != "D17" & OspC != "V1")

# normalize z-score
x.sum =  x %>% group_by(Serum) %>% summarise(mean = mean(OD), sd = sd(OD))
y <- x %>% left_join(x.sum, "Serum") %>% mutate(odScaled = (OD-mean)/sd) %>% select(OspC,Serum,odScaled,Host) %>% rename(OD=odScaled)

y = y %>% left_join(os, "OspC")

z = y %>% select(allele,Serum,OD) %>%  spread(key=Serum, value = OD)
rownames(z) = z$allele
z = z %>% select(-allele)
z = as.data.frame(z)

pheo = y %>% select(Serum, Host) %>% distinct_all()
rownames(pheo) = pheo$Serum
pheo = pheo %>% select(-1)
pheo = as.data.frame(pheo)

alle = y %>% select(allele, class) %>% distinct_all()
rownames(alle) = alle$allele
alle = alle %>% select(class)
alle = as.data.frame(alle)

anno_colors <- list(Host = c(P.leucopus = "#a06400", Human = "#d9d9d9"), class = c(Native="#aaaaaa", Root="#E69F00", Consense="#56B4E9", Centroid="magenta"))

z %>% pheatmap(col=colorpanel(360,muted("blue"),"white","orange"), annotation_col = pheo, annotation_row = alle, annotation_colors = anno_colors)

#z %>% heatmaply(scale="column", seriate = "none", col=colorpanel(360,muted("blue"),"white","orange"), column_text_angle = 90, xlab = 'Sera', ylab = 'Recombinant OspCs', ColSideColors = pheo[colnames(z),], RowSideColors = alle[rownames(z),])
