#################
# Fig 2. Ivano data reanalysis
# OspC-specific serum: scaling with z-scores
###############
library(tidyverse)
library(readxl)
library(scales)
setwd("/Users/lai/Dropbox/ag-div/")
x <- read_xlsx("Sup-data.xlsx", sheet = 2)
x.long <- x %>% pivot_longer(3:18, names_to = "antigen", values_to =  "OD")

# compare between reps
x.rep <- x.long %>% pivot_wider(names_from = "rep", values_from = "OD")
colnames(x.rep) <- c("serum", "antigen", "rep1", "rep2")
x.rep %>% ggplot(aes(x = rep1, y = rep2)) + geom_point(shape=1) + theme_bw() + facet_wrap(~serum) + geom_abline(slope = 1, intercept = 0, linetype =2)

# boxplot for unscaled OD
# by serum:
x.long <- x.long %>% mutate(bind.class = if_else(str_remove(serum, "^.") == str_remove(antigen, "^."), "homol", "heter"))
x.long %>% filter(rep==1) %>% ggplot(aes(x=serum, y=OD)) + geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("OspC-specific serum") + geom_jitter(shape=1, aes(color = bind.class), size=2, width = 0.25) + theme(legend.position = "none") + scale_color_manual(breaks = c("homol", "heter"), values = c("black", "gray"))

# by antigen:
x.long.rep1 <- x.long %>% filter(rep == 1)
lm.sera <- lm(data = x.long.rep1, OD ~ serum)
anova(lm.sera)

lm.ag <- lm(data = x.long.rep1, OD ~ antigen)
anova(lm.ag)

x.med <- x.long.rep1 %>% group_by(antigen) %>% summarise(med.od = median(OD)) %>% arrange(med.od) %>% pull(antigen)

x.long.rep1$antigen <- factor(x.long.rep1$antigen, levels = x.med)  
x.long.rep1 %>% ggplot(aes(x=antigen, y=OD)) + geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("rOspC") + geom_jitter(shape=1, aes(color = bind.class), size=2, width = 0.25) + theme(legend.position = "none") + scale_color_manual(breaks = c("homol", "heter"), values = c("black", "gray"))


# scale by serum and replicate
x.mean.sd <- x.long %>% group_by(serum, rep) %>% summarise(meanOD = mean(OD), sdOD = sd(OD))
x.long2 <- x.long %>% left_join(x.mean.sd, c("serum", "rep"))
x.long2 <- x.long2 %>% mutate(scaledOD = (OD - meanOD)/sdOD, centerOD = OD - meanOD)

x.long2.rep1 <- x.long2 %>% filter(rep == 1)
lm.ag <- lm(data = x.long2.rep1, scaledOD ~ antigen)
anova(lm.ag)

lm.sera <- lm(data = x.long2.rep1, scaledOD ~ serum)
anova(lm.sera)


# obtain mean & sd for each variant
x.mean <- x.long2 %>% group_by(antigen, serum) %>% summarise(meanScaled = mean(scaledOD), sdScaled = sd(scaledOD), meanCenter = mean(centerOD), sdCenter = sd(centerOD))

# boxplots scaled OD to check
x.mean <- x.mean %>% mutate(bind.class = if_else(str_remove(serum, "^.") == str_remove(antigen, "^."), "homol", "heter"))
x.mean %>% ggplot(aes(x=serum, y=meanScaled, label = antigen)) + geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("OspC-specific serum") + geom_jitter(shape=1, aes(color = bind.class), size=2, width=0.25) + theme(legend.position = "none") + geom_hline(yintercept = 0, linetype = 2) + scale_color_manual(breaks = c("homol", "heter"), values = c("black", "gray"))

x.med <- x.mean %>% group_by(antigen) %>% summarise(med.od = median(meanScaled)) %>% arrange(med.od) %>% pull(antigen)

# an alternative to ARC:
x.mean$antigen <- factor(x.mean$antigen, levels = x.med)  
x.mean %>% ggplot(aes(x=antigen, y=meanScaled)) + geom_boxplot(outlier.shape = NA) + theme_bw() + geom_jitter(shape=1, size=2, aes(color = bind.class), width = 0.25) + xlab("rOspC") + geom_hline(yintercept = 0, linetype=2) + theme(legend.position = "none") +  scale_color_manual(breaks = c("homol", "heter"), values = c("black", "gray"))

x.mean %>% ggplot(aes(x=antigen, y=meanScaled)) + geom_violin() + theme_bw() + geom_jitter(color= "green", width = 0.25) + xlab("rOspC") 

# check high sd for sK:

x.mean %>% ggplot(aes(x=serum, y=sdScaled)) + geom_violin() + theme_bw() + geom_jitter(color= "green", width = 0.25)  

x.mean %>% ggplot(aes(x=meanScaled, y=sdScaled)) + geom_point() + theme_bw() + facet_wrap(~serum)


# plot bar graph
x.mean %>% ggplot(aes(x=serum, y=meanScaled)) +  theme_bw() + geom_bar(stat = "identity", color = "black", fill = "gray") + geom_errorbar(aes(ymin=meanScaled - sdScaled, ymax = meanScaled + sdScaled), color="red", width = 0.3)  + geom_hline(yintercept = 2, color="black", linetype="dashed") + facet_wrap(~antigen, nrow = 4) + ylab("z-score (scaled OD450)") + xlab("OspC-specific sera")

# plot heatmap

x.mean %>% ggplot(aes(x = serum, y = antigen, fill = meanScaled)) + geom_tile() + scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"))

#######
# ROC plot as specificity
##############
ag <-levels(factor(x.mean$antigen))
output <- vector("list")
for (i in seq_along(ag)) {
  output[[i]] <- x.mean %>% filter(antigen == ag[i]) %>% arrange(meanScaled) %>% mutate(od_cum = cumsum(meanScaled)) %>% mutate(order = 1:15) 
}
cum.df <- bind_rows(output)
cum.max <- cum.df %>% group_by(antigen) %>% summarise(max = max(order)) %>% left_join(cum.df, c("antigen", "max" = "order"))

ggplot(data=cum.df, aes(x=order, y=od_cum, label = serum, color=antigen, group=antigen)) + geom_line() + geom_point(size=2, shape=1) +  theme_bw() + geom_hline(yintercept = 0, linetype="dashed") + geom_text_repel(data = cum.max, aes(x=max, y=od_cum, label=antigen), size=5)  + theme(legend.position = "none") + ylim(-20,13) + xlab("Number of OspC-specific sera") + ylab("Cumulative OD (scaled)")

