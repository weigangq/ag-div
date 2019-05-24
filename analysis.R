library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)
library(MASS)
library(ape)
setwd("../Dropbox/QiuDi/ospc-sim/")

# linearize logistic growth
x <- read.table("cliff-walk.tsv", header = F)
x.long <- melt(x, measure.vars = 3:102)
colnames(x.long) <- c("allele", "pos", "id", "fit")
x.def <- x.long %>% filter(fit < 1 & fit > 0)
x.def <- mutate(x.def, f = fit/(1-fit))
x.def <- mutate(x.def, prob = pos/106)
x.list <- split(x.def, x.def$allele)
x.lm <- lapply(x.list, function(y) 
  { summary(lm(data = y, f ~ prob)) })
cliff.walk.rates <- lapply(x.lm, function(z) {z$coefficients[2,]})
cliff.rates <- do.call(rbind, cliff.walk.rates)

ggplot(data = x.def, aes(x=pos/106, y=f)) + geom_point(shape=1, color="gray") + scale_y_log10() + facet_wrap(~allele) + geom_smooth(method="lm", colour = "red", linetype="dashed") + geom_hline(yintercept = 1, linetype="dashed") + geom_vline(xintercept = 0.5, linetype = "dashed")

# plot LON
library(igraph)
links <- read.csv("network-oc-edge.csv")
nodes <- read.csv("network-oc-node.csv")
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
V(net)$size <- V(net)$cw50 * 75
V(net)$color <- rep("white",16)
E(net)$width <- 1
E(net)$weight <- E(net)$pw50 * 100
force.dir <- layout_with_fr(net)
#layout(c(1,1,2,2), byrow=T)
par(mfrow=c(1,2))
net.1 <- delete_edges(net, E(net)[pw50>=0.3])
plot(net.1, layout = force.dir, main=expression(d<0.3), edge.curved=0)
#plot(force.dir$extd_graph, main=expression(d<0.3), edge.curved=0)
net.2 <- delete_edges(net, E(net)[pw50>=0.4])
plot(net.2, layout = force.dir, main=expression(d<0.4))

# Hamming bound
perm_without_replacement <- function(n, r){
  return(factorial(n)/factorial(n - r)/factorial(r))
}

sum_error <- function(n, d) {
  den <- 0;
  error <- floor((d-1)/2);
  for(k in 0:error) {
    den <- den + perm_without_replacement(n,k)
  }
  den;
}
  

q <- 2; # binary code (Hamming code)
n <- 7; # block size
d <- 1:6; # d(min)
A <- sapply(d, function(x) { # max number
  den <- 0;
  for(k in 0:floor((x-1)/2)) {
    den <- den + perm_without_replacement(n,k) * (q-1)^k
  }
  q^n/den
})

# perfect Hamming code (n,e)
# log2(A) = n - log2(sum(M))


# map for simulated alleles6
sim <- matrix(c(0,6,6,6,6,4,6,0,6,6,6,4,6,6,0,6,6,4,6,6,6,0,6,4,6,6,6,6,0,4,4,4,4,4,4,0), nrow = 6, byrow = T)
rownames(sim) <- c("A1", "A2", "A3", "A4", "A5", "XR")
colnames(sim) <- c("A1", "A2", "A3", "A4", "A5", "XR")
sim.d <- as.dist(sim, diag = T)
sim.cmd <- cmdscale(sim.d, k = 3)
colnames(sim.cmd) <- c("d1", "d2", "d3")
library(scatterplot3d)
sim.3d <- scatterplot3d(sim.cmd, type="h", angle = 120)
sim.3d$points(sim.cmd, cex=rep(2,6), col=c(1,1,1,1,1,2), pch=rep(16,6))
sim.coord <- sim.3d$xyz.convert(sim.cmd)
text(sim.coord, labels = rownames(sim.cmd), pos=2, offset = 0.5)

library(plot3D)
text3D(sim.cmd[,1], sim.cmd[,2], sim.cmd[,3], labels  = rownames(sim.cmd), colvar = NULL, phi=20, col = c(1,1,1,1,1,2), font = 2, bty="g", adj=1)

scatter3D(sim.cmd[,1], sim.cmd[,2], sim.cmd[,3], add=T, type = "h", cex = 2, colvar = NULL, col = c(1,1,1,1,1,2), pch=19)

x <- read.table("../../LabShared/drylab-data/QiufromSulkow/pw_avg_Mat.txt", sep=" ")

x.dist <- as.dist(x)
pw <- melt(x, variable.name = "pw50")
pw <- mutate(pw, a2 = rep(rownames(x),16))
colnames(pw) <- c("a1", "pw50", "a2")

# plot distances
x <- read.table("dist.tsv5", header=T, sep="\t")
ggplot(data = x, aes(x=tree.dist, y=array.cor)) + geom_point(shape=1, size = 3) + theme_bw() + geom_smooth(method = "lm") + theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20)) + geom_hline(yintercept = 0, linetype=2)
summary(lm(data=x, array.cor ~ tree.dist))

ggplot(data = x, aes(x=seq.diff, y=array.cor)) + geom_point(shape=1, size = 3) + theme_bw() + geom_smooth(method = "lm") + theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20)) + geom_hline(yintercept = 0, linetype=2)
summary(lm(data=x, array.cor ~ seq.diff))
  
qqnorm(x$p.more)
qqline(x$p.more)
qqline(x$p.less)

# plot cor for each allele
x.ale <- lapply(oc$tip.label, function (a) { filter(x, grepl(a, pair)) %>% mutate(allele = rep(a,15)) })
x.ale.df <- do.call(rbind, x.ale)
ggplot(data=x.ale.df, aes(x=seq.diff, y=array.cor)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~allele, ncol = 4) + geom_hline(yintercept = 0, linetype=2)

library(broom)
x.lm <- x.ale.df %>% group_by(allele) %>% do(tidy(lm(array.cor ~ seq.diff, data=.)))

# CW50 analysis
x.min.seq <- x.ale.df %>% group_by(allele) %>% summarize(min(seq.diff))
cw50.brian$'min.seq.diff' <- x.min.seq$`min(seq.diff)`
ggplot(data=cw50.brian, aes(x=min.seq.diff, y=fraction, label=allele)) + geom_point() + geom_smooth(method = "lm") + geom_label() + theme_bw() + ylim(0.3, 0.7)
summary(lm(data=cw50.brian, fraction ~ min.seq.diff))

#x.min.cor <- x.ale.df %>% group_by(allele) %>% summarize(min(array.cor))
#cw50.brian$'min.cor' <- x.min.cor$`min(array.cor)`
#ggplot(data=cw50.brian, aes(x=min.cor, y=fraction, label=allele)) + geom_point() + geom_smooth(method = "lm") + geom_label() + theme_bw() + ylim(0.3,0.7)
#summary(lm(data=cw50.brian, fraction ~ min.cor))

#x.ale.closest <- lapply(x.ale, function(x) {arrange(x, seq.diff) %>% head(n=1)} )
#x.ale.cor <- lapply(x.ale, function(x) {c <- cor(x$seq.diff, x$array.cor); data.frame(al = x$allele[1], cor=c)})
#x.ale.df3 <- do.call(rbind, x.ale.cor) %>% arrange(al) # doesn't work

get.cor <- function(a,b) {cor(a,b)}
x.ale.cor <- x.ale.df %>% group_by(allele) %>% summarise(cor = get.cor(seq.diff, array.cor))
cw50.brian$'avg.cor' <- x.ale.cor$cor
ggplot(data=cw50.brian, aes(x=avg.cor, y=fraction, label=allele)) + geom_point() + geom_smooth(method = "lm") + geom_label() +  theme_bw() + ylim(0.3, 0.7)
summary(lm(data=cw50.brian, fraction ~ avg.cor))

#x.ale.mc <- x.ale.df %>% group_by(allele) %>% summarise(mean.cor = mean(array.cor))
#cw50.brian$'mean.cor' <- x.ale.mc$'mean.cor'
#ggplot(data=cw50.brian, aes(x=mean.cor, y=fraction, label=allele)) + geom_point() + geom_smooth(method = "lm") + geom_label() +  theme_bw() + ylim(0.3, 0.7)
#summary(lm(data=cw50.brian, fraction ~ mean.cor))

#x.ale.df2 <- do.call(rbind, x.ale.closest) %>% arrange(allele)
#cw50.brian$'min.seq.cor' <- x.ale.df2$array.cor
#ggplot(data=cw50.brian, aes(x=min.seq.cor, y=fraction, label=allele)) + geom_point() + geom_smooth(method = "lm") + geom_label() +  theme_bw() + ylim(0.3, 0.7)
#summary(lm(data=cw50.brian, fraction ~ min.cor))


#+  scale_colour_gradient(high="gold", low ="steelblue") 

ggplot(data = x, aes(x=array.cor, y= p.more)) + geom_point(size = 3)  + theme_bw() + geom_smooth(method = "lm")


x <- read.table("oc-tree-dist.txt", header=F, sep= " ")
colnames(x) <- c("pair", "diff", "len", "iden", "tree.dist")

ggplot(data = x, aes(x = tree.dist, y=diff/len)) + geom_point() + geom_smooth(method = "lm")


# NE genome tree
library(ape)
tr <- read.tree("bbss-cp26.us.dnd2")
tab <- read.csv("bbss-cp26.us.tab3", header=T, row.names = 2)
par.tr <- plot.phylo(tr, edge.width = 3, x.lim = 0.2, show.tip.label = F, no.margin = T)
#nodelabels(round(as.numeric(tr$node.label),2), adj = -0.2, frame = "n")
add.scale.bar()
text(rep(0.09,13), 1:13, tr$tip.label, pos = 4)
text(rep(0.12,13), 1:13, tab[tr$tip.label,2], pos = 4)
text(rep(0.15,13), 1:13, tab[tr$tip.label,3], pos = 4)
text(rep(0.16,13), 1:13, tab[tr$tip.label,4], pos = 4)

# ospC tree
library(ape)
oc <- read.tree("ospC-NE2.dnd")
plot(oc, type = "u")

# ospC pair diff
pd <- read.table("ospC-NE.pair-diff", header = F, sep="\t")
pd$all1 <- gsub("(.+)\\|.+", "\\1", pd$V1, perl=T)
pd$all2 <- gsub("(.+)\\|.+", "\\1", pd$V2, perl=T)
pd.tall <- data.frame(pd$all1, pd$all2, diff=pd$V3/pd$V5)

pd.dist <- as.dist(xtabs(pd.tall[, 3] ~ pd.tall[, 2] + pd.tall[, 1]))
pd.cmd <- cmdscale(pd.dist)
colnames(pd.cmd) <- c("x", "y")
ggplot(data = as.data.frame(pd.cmd), aes(x=x, y=y, label = rownames(pd.cmd))) + geom_label(size=5) + theme_bw()

# not working
# pd.wide <- reshape(pd.tall, direction = "wide", idvar = "pd.all1", timevar = "pd.all2")
# pd.d <- as.dist(pd.wide[,-1])
# attr(pd.d, "Labels") <- pd.d[,1]


# pair-walk
x <- read.table("pair-walk.tsv", header = F)
dim(x)
head(x)
#x.a <- filter(x, x$V1=='A' | x$V2=='A')
x.long <- melt(x, measure.vars = 5:34)
colnames(x.long) <- c("from", "to", "num.mut", "diff.site", "dummy", "fit")

x.long <- mutate(x.long, prob = num.mut/diff.site, pair = paste(from, to, sep=""))
x.def <- x.long %>% filter(fit < 1 & fit > 0)
x.def <- mutate(x.def, f = fit/(1-fit))
x.list <- split(x.def, x.def$pair)
x.lm <- lapply(x.list, function(y) 
{ summary(lm(data = y, f ~ prob)) })
pair.walk.rates <- lapply(x.lm, function(z) {z$coefficients[2,]})
pair.rates <- as.data.frame(do.call(rbind, pair.walk.rates))

pair.rates.order <- pair.rates[order(pair.rates$Estimate),]
#arrange(pair.rates, Estimate) %>% head()
#arrange(pair.rates, Estimate) %>% tail()

x.list <- split(x.long, x.long$pair)

x.mods <- lapply(x.list, function(x) glm(fit ~ prob, family=binomial(link=logit), data=x))
ld50 <- lapply(x.mods, function(m) dose.p(m,p=0.5))

out.df <-do.call("rbind", ld50)
write.table(out.df, "ld50.pair", quote = F, col.names = F)



x.a <- mutate(x, mean=apply(x, 1, function(p) mean(as.numeric(p[5:34]))), sd=apply(x, 1, function(p) sd(as.numeric(p[5:34]))))
x.a.df <- x.a[c(1:4,35:36)]
colnames(x.a.df) <- c("from", "to", "pos", "total", "fit", "sd")
ggplot(x.a.df, aes(x=pos/total, y=fit, color=to)) + geom_line() +
  geom_vline(xintercept=0.5, linetype = "dashed") +
  geom_hline(yintercept=0.5, linetype = "dashed") + theme_bw() +
  facet_wrap(~from)


x.a.df <- mutate(x.a.df, panel = ifelse(from=='A', LETTERS[to], LETTERS[from]))

x.a.df <- mutate(x.a.df, dir = ifelse(from=='A', "to", "from"))

ggplot(x.a.df, aes(x=pos/total, y=fit, color=dir)) + geom_line() +
#  geom_ribbon(aes(ymax=fit+sd, ymin=fit-sd), alpha=0.2) + 
  geom_vline(xintercept=0.5, linetype = "dashed") + theme_bw() +
 facet_wrap(~panel)



# plot cliff walk
x <- read.table("cliff-walk.tsv", header = F)
x.long <- melt(x, measure.vars = 3:102)
colnames(x.long) <- c("allele", "pos", "id", "fit")
cw50.brian <- read.table("cw50.tsv", header=T, sep="\t")
cw50.brian <- mutate(cw50.brian, x=0.5)
cw50 <- data.frame(pos=cw50.brian$fraction*106, fit=cw50.brian$x, allele=cw50.brian$allele)

x.2 <- mutate(x, mean=apply(x, 1, function(p) mean(as.numeric(p[3:102]))), sd=apply(x, 1, function(p) sd(as.numeric(p[3:102]))))
x.3 <- select(x.2, c(1,2,103,104))
colnames(x.3) <- c("allele", "pos", "mean", "sd")
 ggplot(data=x.long, aes(x=pos/106, y = fit)) + geom_point(alpha=1/20, colour="gray") + theme(axis.text.x = element_text(face="bold", size=14),  axis.text.y = element_text(face="bold", size=14)) + theme_bw() + facet_wrap(~allele) + geom_line(data=x.3, aes(x=pos/106, y = mean)) + geom_hline(yintercept=0.5, linetype = "dashed") + geom_point(data=cw50, size=3, aes(x=pos/106, y=fit))


ggplot(x.3, aes(x=pos/106, y=mean)) + geom_line() +
#geom_ribbon(aes(ymax=mean+sd, ymin=mean-sd), alpha=0.5) + 
  geom_vline(xintercept=0.5, linetype = "dashed") + facet_wrap(~allele)

# estimate ld50
x.long <- mutate(x.long, prob = pos/106)
x.list <- split(x.long, x.long$allele)
x.mods <- lapply(x.list, function(x) glm(fit ~ prob, family=binomial(link=logit), data=x))
ld50 <- lapply(x.mods, function(m) dose.p(m,p=0.5))

# plot mouse data
mouse <- read.table("baum-mouse.tsv", header=T)
m.df <- mouse[,c(3,6,7,10)]
#arr, infected, logInt, logInt_norm)
m.melt <- melt(m.df, measure.vars=3:4)
colnames(m.melt) <- c("array", "infected.allele", "group", "intensity")
m.df <- mutate(m.melt, norm=ifelse(group=="logInt", "pre", "post"))
ggplot(subset(m.df[,c(1,2,4,5)], array %in% c("A", "N", "F", "I", "T", "K")), aes(x=infected.allele, y=intensity, color=norm))  + geom_boxplot(position=position_dodge(0.8)) + geom_jitter(position=position_jitterdodge(0.8)) + facet_wrap(~array, nrow = 2) +theme(legend.position = "bottom") + theme_bw()




# plot rates for sites
x <- read.table("rate.tsv", header=F)
rate1 <-x
colnames(x) <- c("start", "end", "pos", "aa", "rate", "lwr", "upr","sd", "nseq", "is.indel")

idx.indel <- which(x$is.indel=='t')


idx.hi <- which(x$rate > 2)
x.hi <- x[idx.hi,]
df.label <- data.frame(x=x.hi$start, y=x.hi$rate, lab=paste(x.hi$start, x.hi$aa, sep="-"))

ggplot(subset(x, is.indel=='f'), aes(x=start, y=rate)) + 
  geom_pointrange(aes(ymin=lwr, ymax=upr), colour="black") +
  geom_point(size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("Alignment Position") +
  ylab("Rate") +
  geom_abline(intercept = c(0,1,2), slope = 0, linetype="dashed") + theme_bw() + annotate("text", x=x[idx.indel, 1], y=0, label = "Delta", parse =T, color="red", size=6) +
  geom_text_repel(
    data = subset(x, rate > 2 & is.indel !='t'),
    aes(label = paste(start, aa, sep="")),
    size = 4, color="blue") 








setwd("~/Desktop/ospc-simulation/")
load("lia-RData")

x = read.table("tmp", sep="|")
colnames(x) = c("id","arr","bin","intensity","spc","infected")
x$XReact = ifelse(x$bin=='t', 'Yes', 'No')
x = x[,-3]
x$logInt = log2(x$intensity)

total = x

library(dplyr)

#1. normalize sample
#x = filter(total, infected!='0', spc=='human')
x = filter(total, infected!='0')
x$label = paste(x$infected, x$id, sep='.')
#x$label = paste('h', x$id, sep='')
x$label = paste(ifelse(x$infected=='', 'h', as.character(x$infected)), x$id, sep='.')
rownames(x) = paste(x$arr, x$label, sep='|')

y = long2wide(select(x, arr, label, logInt))
y = as.matrix(y)
y.norm = normalize.quantiles(y)
rownames(y.norm) = rownames(y)
colnames(y.norm) = colnames(y)

y = matrix2long(y.norm, 'arr', 'label')
y$label = sub('A.N', 'A,N', y$label)
rownames(y) = paste(y$arr, y$label, sep='|')
x$logInt_norm = y[rownames(x),3]

sample.df = x

#2. normalize control
#ctl = filter(total, infected=='0', spc=='human')
ctl = filter(total, infected=='0')
x = ctl
#x$label = paste('mc', x$id, sep='.')
x$label = paste('c', x$id, sep='')

y = long2wide(select(x, arr, label, logInt))
y = as.matrix(y)
y.norm = normalize.quantiles(y)
rownames(y.norm) = rownames(y)
colnames(y.norm) = colnames(y)

y = matrix2long(y.norm, 'arr', 'label')

rownames(x) = paste(x$arr, x$label, sep='|')
rownames(y) = paste(y$arr, y$label, sep='|')
x$logInt_norm = y[rownames(x),3]

ctl.total = x

#3. mean
ctl.mean = tapply(x$logInt_norm, x$arr, mean)
ctl.sd = tapply(x$logInt_norm, x$arr, sd)

x = split(sample.df, sample.df$arr)
x.list = lapply(names(x), function(a) cbind(x[[a]], FC=x[[a]]$logInt_norm - ctl.mean[a]))
names(x.list) = names(x)

x.df = unsplit(x.list, f = sample.df$arr);
x.df$FC = ifelse(x.df$FC<0, 0, x.df$FC)

sample.df = x.df

#4. test by control std
#x = split(mouse.df2, mouse.df2$arr)
#x.list = lapply(names(x), function(a) cbind(x[[a]], XReact2=ifelse(x[[a]]$FC >= 3*ctl.sd[a], 1, 0))) # can't use "yes" and "no"
#names(x.list) = names(x)

#x.df = unsplit(x.list, f = mouse.df2$arr);
#x.df$XReact2 = ifelse(x.df$XReact2==1, "Yes", "No")
#mouse.df2 = x.df

#4. test by sample std
x.mean = tapply(sample.df$FC, sample.df$arr, mean)
x.sd =   tapply(sample.df$FC, sample.df$arr, sd)

x = split(sample.df, sample.df$arr)
x.list = lapply(names(x), function(a) cbind(x[[a]], XReact3=ifelse(x[[a]]$FC - x.mean[a] >= x.sd[a], 1, 0)))
names(x.list) = names(x)

x.df = unsplit(x.list, f = sample.df$arr);
x.df$XReact3 = ifelse(x.df$XReact3==1, "Yes", "No")
sample.df = x.df

#5 correlation test
x = select(sample.df, label, arr, FC)
x = long2wide(x)
x = as.matrix(x)
cor.total <- cor(x)
cor.total[which(cor.total==1)] = NA;
 # pairwise cor test (with p values)
cor.out <- lapply(1:22, function(i) { lapply((i+1):23, function(j){out <- cor.test(x[,i], x[,j]); data.frame(row.names = paste(colnames(x)[i], colnames(x)[j], sep="-"), est = out$estimate, pval = out$p.value)})})
cor.out <- unlist(cor.out, recursive = F)
cor.names <- lapply(1:22, function(i) { lapply((i+1):23, function(j){paste(colnames(x)[i], colnames(x)[j], sep="-")})})
cor.names <- unlist(cor.names)
cor.df <- unsplit(cor.out, cor.names)
cor.sig <- cor.df[which(cor.df$pval < 1e-2 & cor.df$est < 0.9),]
cor.df$allele_i <- gsub("^(.+)-(.+)$", "\\1", rownames(cor.df))
cor.df$allele_j <- gsub("^(.+)-(.+)$", "\\2", rownames(cor.df))
write.table(cor.df, "pairwise-cor.txt", quote = F, row.names = F, sep="\t")

#5 correlation with sequence diffs
# for i in {1..208}; do j=$(( i+9 )); bioaln -s "$i,$j" ospC.aln | bioaln --pair-diff | sed "s/^/$i $j /"; done > pairwise.diff
# cat pairwise.diff | tr ' ' '\t' > pairwise.diff10

# In Bash: for i in {1..198}; do j=$(( i+19 )); bioaln -s "$i,$j" ospC.aln | bioaln --pair-diff | sed "s/^/$i $j /"; done > pairwise.diff

win.diff <- read.table("pairwise.diff20", sep="\t", header=F);
colnames(win.diff) = c('i','j','a1', 'a2', 'dif', 'len')
win.list <- split(win.diff, win.diff$i);
win.cor <- lapply(win.list, function(x){
  x.df <- data.frame(row.names = paste(x$a1, x$a2, sep="-"), diff = x$dif);
  id <- unique(x$i);
  x.cor <- cor.test(10-x.df$diff, cor.df[rownames(x.df), 'est']);
  x.conf <- x.cor$conf.int;
  data.frame(win = id, cor = x.cor$estimate, min = x.conf[1], max = x.conf[2], pval = x.cor$p.value, row.names=id);
})
win.df <- unsplit(win.cor, names(win.cor));
win.df <- win.df[order(win.df$win),]

#plot(win.df$win, -1*win.df$cor, type = "l", xlab="AA position", ylab= "cor coeff", las=1, main="seq-crossreact correlation")
#abline(h=0, lty=2)

ggplot(win.df20, aes(x=win, y=cor)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_ribbon(aes(x=win, ymax=max, ymin=min), alpha=0.2) +
  geom_line() +
  labs(title="seq-crossreact correlation (win=20)", x="AA position", y="cor coeff") +
  theme_bw()

win.df10$win_size = 10
win.df20$win_size = 20
wins = rbind(win.df10, win.df20)
wins$win_size=factor(wins$win_size,levels=c(20,10))

ggplot(wins, aes(x=win, y=cor)) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_ribbon(aes(x=win, ymax=max, ymin=min), alpha=0.15) +
  geom_line() +
  labs(title="seq-crossreact correlation", x="AA position", y="cor coeff") +
  theme_bw() +
  facet_wrap(~win_size, nrow = 2, labeller = label_both) +
  scale_x_continuous(breaks = round(seq(0, max(wins$win), by=5),1))


#4. plot
library(ggplot2)
library(plotly)
x.df$XReact3=factor(x.df$XReact3,levels=c("Yes","No"))

ggplot(x.df, aes(x=infected, y=2^FC, color=XReact3)) +
  scale_y_log10(name="intensity") +
  geom_jitter(position=position_jitter(0.1), shape=1) +
  facet_wrap(~arr, nrow=4)

ggplot(x.df, aes(x=arr, y=2^FC, color=XReact3)) +
  scale_x_discrete(name ="sample") +
  scale_y_log10(name="intensity") +
  geom_point(shape=1) +
  facet_wrap(~label)

# ggplot(human.df2, aes(x=label, y=2^FC)) +
#   scale_y_log10(name="intensity") +
#   scale_x_discrete(name ="human sample") +
#   geom_boxplot(outlier.colour = NA, color="gray") +
#   geom_jitter(position=position_jitter(0.2), shape=1)

#heatmap
#x = t(x.norm)
#rownames(x) = sub('human\\|', 'h', rownames((x)))

x = select(sample.df, label, arr, FC)
x = long2wide(x)

#traditional
x = as.matrix(x)
heatmap(x, scale="none", col=colorpanel(1000,"white","blue"))

library(heatmaply)
geo = ifelse(substr(colnames((x)),2,2)=='', 'Northeast', 'other')
heatmaply(x, k_row=2, k_col=2, col=colorpanel(65,"white","blue"), margins = c(51,53), fontsize_row = 6, RowSideColors=ifelse(substr(rownames(x),1,2)=='h.', "human", "mouse"), ColSideColors = geo)

x = select(sample.df, label, arr, FC)
colnames(x) = c("name", "variable", "value")
heatmaply(long_data = x, k_row=2, k_col=2, col=colorpanel(65,"white","blue"), margins = c(51,53), fontsize_row = 6)

x = select(sample.df, label, arr, FC)
x = long2wide(x)
geo = ifelse(substr(colnames((x)),2,2)=='', 'Northeast', 'other')

x$spcies = ifelse(substr(rownames(x),1,2)=='h.', "human", "mouse")
x$infected = ifelse(substr(rownames(x),1,2)=='h.', "", gsub("(.+)\\..+$", "\\1", rownames(x)))

heatmaply(x[,-c(24,25)], seriate = "none", k_row=2, k_col=2, col=colorpanel(65,"white","blue"), margins = c(40,25), column_text_angle = 0, RowSideColors=x[,c(24,25)], ColSideColors = geo, showticklabels = c(TRUE, FALSE), xlab = 'arr', ylab = 'sample')


#heatmap of cor
library(RColorBrewer)
#display.brewer.pal(11, "BrBG")
BrBG = colorRampPalette(brewer.pal(11, "BrBG"))

geo.col = ifelse(substr(colnames((cor.total)),2,2)=='', 'Northeast', 'other')
geo.row = ifelse(substr(rownames((cor.total)),2,2)=='', 'Northeast', 'other')
heatmaply(cor.total, k_row=4, k_col=4, col=BrBG, margins = c(22,25), RowSideColors = geo.row, ColSideColors = geo.col, column_text_angle = 0)

#hc.id=hclust(as.dist(1-cor(x, use="pairwise.complete.obs")))
#hc.arr=hclust(as.dist(1-cor(t(x), use="pairwise.complete.obs")))
#heatmap(x, Colv=as.dendrogram(hc.arr), Rowv=as.dendrogram(hc.id), scale="none", col=colorpanel(1000,"white","blue"), trace="none", cexRow=0.6)


#point of cor between human and mice
x = data_frame(cor.human=c(cor.human), cor.mouse=c(cor.mouse))
ggplot(x, aes(x=cor.mouse, y=cor.human)) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_smooth(method = "lm", col="darkgray", lwd=0.5) +
  geom_point(shape=1) +
  theme_bw()


#cliff walk
x = read.table("cliff-walk.txt")
colnames(x) = c("alle", "mutNum", "fit", "std")

p = ggplot(x, aes(x=mutNum, y=fit, color=alle)) +
  geom_hline(yintercept=0.5, linetype = "dashed") +
  geom_line() +
  labs(title="Cliff walk", x="Mutation number", y="Fitness", formula=y ~ x)
ggplotly(p)
p

q = ggplot(x, aes(x=mutNum, y=fit)) +
  geom_hline(yintercept=0.5, linetype = "dashed") +
  geom_line() +
  facet_wrap(~alle, ncol=1) +
  #  geom_smooth(method = "loess", col="darkgray", lwd=0.5, position = "identity") +
  geom_ribbon(aes(x=mutNum, ymax=fit+std, ymin=fit-std), alpha=0.15) +
  labs(title="Cliff walk", x="Mutation number", y="Fitness", formula=y ~ x)
#ggplotly(q)
q

#pair fit
x=read.table("~/Dropbox/Shared_folder/ospc-sim/pair-fit.txt", header = T)
p = ggplot(x, aes(x=num.mutation, y=mean.fit, color=to)) +
  geom_hline(yintercept=0.5, linetype = "dashed") +
  geom_line() +
  labs(x="Mutation number", y="Fitness", formula=y ~ x, title="Pair Walk") +
  facet_wrap(~from, labeller = label_both)
ggplotly(p)

p

# use of dplyr
filter(x, spec=='human', arr == 'A') # filter
select(x, intensity) # select columns
summarise(x, sd(log10(intensity))) # calculate standard deviation
x %>% filter(spec=='human', arr== 'A') %>% select(intensity) %>% summarise(sd(log10(intensity))) # chain the above 3 steps


#
x = read.table("tt")
colnames(x) = c("gn", "allele", "xrect")
ggplot(x, aes(x=gn, y=xrect)) +
  geom_line() +
  labs(x="Generation", y="xRect") +
  facet_wrap(~allele, nrow=16)
