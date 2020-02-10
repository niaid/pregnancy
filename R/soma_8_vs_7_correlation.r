library(corrplot)

# Load data
load("data_generated/Pregnancy_somalogic.RData")

si = info.soma$Visit < 4
sum(si)
info.soma = info.soma[si,]
dat.soma = dat.soma[si,]
colnames(dat.soma) = ann.soma$TargetFullName

# read old 8 somamers
s8 = fread("data/pred.soma.old.txt", sep="\t", header = F) %>% pull(1)

# read new 7 somamers
s7 = fread("data/pred.soma.new.txt", sep = "\t", header = F) %>% pull(1)

ss = c(s8,s7)
i.mod = match(ss, ann.soma$TargetFullName)

dat2 = dat.soma[, i.mod]

cc.mat = cor(dat2)
colnames(cc.mat) = NULL

png("figures/soma_8_vs_7_correlation.png", w=650, h=400)
corrplot(cc.mat, col=colorRampPalette(c("blue","white","red"))(10), tl.col="black")
dev.off()
