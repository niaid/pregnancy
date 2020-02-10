
# Load data ---------------------------------------------------------------

load("data_generated/Pregnancy_flow_cytometry.RData")
load("data_generated/Pregnancy_somalogic.RData")


# load V3 vs. V1 results --------------------------------------------------

soma.comp = fread("results/soma_visits_comparison.txt")
i.soma = soma.comp$Target[soma.comp$wQ <= 0.05 & abs(log2(soma.comp$FC)) > log2(1.2)]
length(i.soma)

flow.comp = fread("results/flow_visits_comparison.txt")
i.flow = flow.comp$ID[flow.comp$wQ.3.v.1 <= 0.05 & abs(log2(flow.comp$FC.3.v.1)) > log2(1.2)]
length(i.flow)


# Soma vs. Flow -----------------------------------------------------------

fi = rownames(dat.flow) %in% rownames(dat.soma)
si = match(rownames(dat.flow)[fi],rownames(dat.soma))
identical(rownames(dat.soma)[si],rownames(dat.flow)[fi])

dat.flow = dat.flow[fi, i.flow]
info.flow = info.flow[fi,]
dat.soma = dat.soma[si, i.soma]
info.soma = info.soma[si,]

tp1 = 1
tp2 = 3

tp.lbl = paste0("V",tp2,".V",tp1)
i2 = which(info.soma$Visit == tp2)
i1 = which(info.soma$Visit == tp1)
is2 = info.soma$Subject[i2] %in% info.soma$Subject[i1]
is1 = match(info.soma$Subject[i2], info.soma$Subject[i1])
stopifnot(identical(info.soma$Subject[i2[is2]], info.soma$Subject[i1[is1]]))
dat.soma.tp = dat.soma[i2[is2],] / dat.soma[i1[is1],]
info.soma.tp = info.soma[i2,]

i2 = which(info.flow$Visit == tp2)
i1 = which(info.flow$Visit == tp1)
is2 = info.flow$Subject[i2] %in% info.flow$Subject[i1]
is1 = match(info.flow$Subject[i2], info.flow$Subject[i1])
stopifnot(identical(info.flow$Subject[i2[is2]], info.flow$Subject[i1[is1]]))
dat.flow.tp = dat.flow[i2[is2],] / dat.flow[i1[is1],]
info.flow.tp = info.flow[i2,]


ccs.sf = cor(dat.soma.tp, dat.flow.tp, use = "pairwise.complete.obs", method = "spearman")
ccs.sf.p = ccs.sf
ccs.sf.padj = ccs.sf
for(m in 1:ncol(dat.flow.tp)) {
  for(k in 1:ncol(dat.soma.tp)) {
    ccs.sf.p[k,m] = cor.test(dat.soma.tp[,k], dat.flow.tp[,m], method = "spearman",
                             alternative = "two.sided", exact=F)$p.value
  }
  ccs.sf.padj[,m] = p.adjust(ccs.sf.p[,m], method="BH") %>% matrix(dim(ccs.sf))
}
rownames(ccs.sf.padj) = rownames(ccs.sf.p)
colnames(ccs.sf.padj) = colnames(ccs.sf.p)

ccs.sf.p.log = ccs.sf
ccs.sf.p.log = -log10(ccs.sf.p) * sign(ccs.sf)
range(ccs.sf.padj)

ccs.sf.row.min = apply(ccs.sf.padj, 1, min)
i.sel = ccs.sf.row.min <= 0.05
sum(i.sel)

ccs.sf.col.min = apply(ccs.sf.padj[i.sel,], 2, min)
i.sel.fl = ccs.sf.col.min <= 0.1
sum(i.sel.fl)


DF = data.frame()
DF.ccs = data.frame()

visit_lbl = paste0("V",tp2,".V",tp1)

df.ccs.p = ccs.sf.p %>% as.data.frame() %>% 
  tibble::rownames_to_column("soma") %>% 
  gather("population", "p", -1)
df.ccs.padj = ccs.sf.padj %>% as.data.frame() %>% 
  tibble::rownames_to_column("soma") %>% 
  gather("population", "FDR", -1)
df.ccs = ccs.sf %>% as.data.frame() %>% 
  tibble::rownames_to_column("soma") %>% 
  gather("population", "rho", -1) %>% 
  mutate(p = df.ccs.p$p, FDR = df.ccs.padj$FDR, visit.label = visit_lbl) %>% 
  arrange(FDR,p)


fwrite(df.ccs, file = glue::glue("results/Soma_Flow_{visit_lbl}_spearman.txt"), sep="\t")


DF = cbind(dat.soma.tp, dat.flow.tp) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sample") %>% 
  separate(sample, c("subject","visit"), remove = F)

popt = df.ccs$population[1]
somat = df.ccs$soma[1]
pop = glue::glue("`{popt}`")
somam = glue::glue("`{somat}`")

df.text = df.ccs %>% 
  dplyr::filter(population == popt, soma == somat) %>% 
  mutate(label = glue::glue("p = {format(p, digits=2)}\nFDR = {format(FDR, digits=2)}"))
ggplot(DF, aes_string(somam, pop)) +
  geom_smooth(method = "lm", fill = NA) +
  geom_point(size=2) +
  # facet_wrap(~visit.label) +
  scale_y_log10() +
  # geom_text(aes(label=subject), size=3) +
  geom_text(data=df.text, aes(label=label), x=Inf, y=-Inf, hjust=1.1, vjust=-0.5) +
  guides(size = guide_legend(), col = F) +
  xlab(somat) +
  ylab(paste("population ", popt)) +
  theme_bw()
fn.fig = glue::glue("figures/ID{popt}_{somat}_FC_{visit_lbl}_correlation.png")
ggsave(fn.fig, w=5, h=4)

