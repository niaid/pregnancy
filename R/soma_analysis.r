library(readxl)
library(readat)

# load subject to use
fn.clin = "data/Pregnancy_clinical_data.xlsx"
df.si = read_xlsx(fn.clin) %>%
  mutate(Subject = sub("0*", "", Subject))

fn.soma = "data/Pregnancy.hybNorm.plateScale.calibrate.20170512.adat"
adat = readAdat(fn.soma)

ann.soma = getSequenceData(adat)
info.soma = getSampleData(adat) %>% 
  dplyr::rename(Subject=SampleDescription) %>% 
  mutate(Subject = sub("Pregnancy0*","",Subject))

# remove samples not related to the study
si = !info.soma$TimePoint %in% c(1,7) & info.soma$Subject %in% df.si$Subject
sum(si)

info.soma = info.soma[si,] %>% 
  group_by(Subject) %>% 
  mutate(Visit = rank(TimePoint)) %>% 
  ungroup() %>% 
  mutate(sample = paste(Subject, Visit, sep="_"))

dat.soma = readat::getIntensities(adat)[si,]
rownames(dat.soma) = info.soma$sample
colnames(dat.soma) = ann.soma$Target

# reorder by Subjects and Visits
info.soma = info.soma %>% arrange(as.numeric(Subject), Visit)
dat.soma = dat.soma[match(info.soma$sample, rownames(dat.soma)),]

dir.create("data_generated", showWarnings = F)
save(info.soma, dat.soma, ann.soma, file="data_generated/Pregnancy_somalogic.RData")

# save all data to text file
dat.soma.df = info.soma %>% 
  mutate(Subject = as.numeric(Subject)) %>% 
  cbind(as.data.frame(dat.soma)) %>% 
  arrange(Subject, Visit)

# remove columns with the same information for all samples
irm = sapply(dat.soma.df, function(x) length(unique(x))) == 1

write.table(dat.soma.df[,!irm], file = "data_generated/Pregnancy_soma_data.txt", sep="\t", quote=F, row.names = F)
write.table(ann.soma, file = "data_generated/Pregnancy_soma_somamers.txt", sep="\t", quote=F, row.names = F)



# compare Somalogic protein intensity at visit 3 vs. visit 1 --------------

fc.soma = c()
wt.soma = list()
for (ci in colnames(dat.soma)) {
  fc.soma[ci] = mean(dat.soma[info.soma$Visit==3,ci] / dat.soma[info.soma$Visit==1,ci], na.rm=T)
  wt.soma[[ci]] = wilcox.test(dat.soma[info.soma$Visit==1,ci], dat.soma[info.soma$Visit==3,ci], paired=T, exact = F, conf.int = TRUE)
}
wt.soma.p = map_dbl(wt.soma, function(x) x$p.value)
wt.soma.padj = p.adjust(wt.soma.p, method = "BH")
wt.soma.estimate = map_dbl(wt.soma, function(x) x$estimate)
i.soma = names(wt.soma)[wt.soma.padj<0.05]
length(i.soma)

S.df = data.frame(Target = colnames(dat.soma), FC = fc.soma, wP = wt.soma.p, wQ = wt.soma.padj)
fwrite(S.df, "results/soma_visits.3.v.1_comparison.txt", sep="\t")

S.df.sel = S.df %>% dplyr::filter(wQ <= 0.05) %>% arrange(wP)
fwrite(S.df.sel, "results/soma_visits.3.v.1_comparison_selected.txt", sep="\t")
