####################################
# Immunoblot analysis
####################################
library(tidyverse)
library(readxl)
library(scales)
setwd("/Users/lai/Dropbox/ag-div/")
#################
# 10/13/2020 
# Immunoblot (Roman & Justin)
###############

x <- read_excel("Sup-data.xlsx", sheet = 3)
x <- x %>% mutate(od_adj = if_else(od_raw - bsa > 0, od_raw - bsa, 0))

# boxplot for unscaled OD
x %>% ggplot(aes(x=serum, y=od_adj)) + geom_boxplot() + theme_bw() + geom_jitter(color= "orange", shape = 1) + xlab("OspC-specific serum") + facet_wrap(~species)

# scale by serum and replicate
x.mean.sd <- x %>% group_by(serum, species) %>% summarise(meanOD = mean(od_adj), sdOD = sd(od_adj))
x2 <- x %>% left_join(x.mean.sd, c("serum", "species"))
x2 <- x2 %>% mutate(scaledOD = (od_adj - meanOD)/sdOD)

x2 %>% ggplot(aes(x=serum, y=scaledOD)) + geom_boxplot() + theme_bw() + geom_jitter(color= "orange", shape = 1) + xlab("OspC-specific serum") + facet_wrap(~species)

# plot heatmap
x2$antigen = factor(x2$antigen, levels = rev(levels(factor(x2$antigen))))
x2 %>% ggplot(aes(x=serum, y = antigen, fill = scaledOD)) +  geom_tile() + theme_bw() + scale_fill_gradient(low = "white", high = "blue") + facet_wrap(~species)

# plot bar graph
x2 %>% filter(species == 'C3H') %>% ggplot(aes(x=serum, y=scaledOD)) +  theme_bw() + geom_bar(stat = "identity", color = "black", fill = "gray") + geom_hline(yintercept = 2, color="black", linetype="dashed") + facet_wrap(~antigen, nrow = 2) + ylab("z-score") + xlab("OspC-specific sera")

x2 %>% filter(species == 'C3H') %>% ggplot(aes(x=antigen, y=scaledOD)) +  theme_bw() + geom_bar(stat = "identity", color = "black", fill = "gray") + geom_hline(yintercept = 2, color="black", linetype="dashed") + facet_wrap(~serum, nrow = 4) + ylab("z-score") + xlab("OspC-specific sera")

ag <-levels(factor(x2$antigen))
output1 <- vector("list")
for (i in seq_along(ag)) {
  output1[[i]] <- x2 %>% filter(antigen == ag[i] & species == 'C3H') %>% arrange(scaledOD) %>% mutate(od_cum = cumsum(scaledOD)) %>% mutate(order = 1:length(ag), species = 'C3H')
}

output2 <- vector("list")
for (i in seq_along(ag)) {
  output2[[i]] <- x2 %>% filter(antigen == ag[i] & species != 'C3H') %>% arrange(scaledOD) %>% mutate(od_cum = cumsum(scaledOD)) %>% mutate(order = 1:length(ag), species = 'Pleuc') 
}

cum.df <- bind_rows(output1, output2)
cum.max <- cum.df %>% group_by(antigen, species) %>% summarise(max = max(order)) %>% left_join(cum.df, c("antigen", "species", "max" = "order"))

ggplot(data=cum.df, aes(x=order, y=od_cum, label = serum, color = antigen, group = antigen)) + geom_line(size=1) + theme_bw() + geom_hline(yintercept = 0, linetype="dashed") + geom_text_repel(data = cum.max, aes(x=max, y=od_cum, label=antigen), size=5)  + xlab("Number of simulated sera") + ylab("Cumulative OD (scaled)") + facet_wrap(~species)
