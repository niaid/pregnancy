library(readxl)

# load subject to use
fn.clin = "data_generated/Clinical_data_v123.txt"
df.si = fread(fn.clin) %>%
  mutate(Subject = sub("0*", "", Subject))

# load Somalogic data
load("data_generated/Pregnancy_somalogic.RData")

si = info.soma$Visit < 4
sum(si)
info = info.soma[si,] %>% 
  left_join(df.si %>% dplyr::select(Subject, Visit, weeks), by=c("Subject","Visit"))


# load genes from previous model
fn.m = "data/Model_genes.xlsx"
df.mod = read_xlsx(fn.m)

i.mis = which(!toupper(df.mod$`Target full name`) %in% toupper(as.character(ann.soma$TargetFullName)))
df.mod$`Target full name`[i.mis]

i.mod = toupper(ann.soma$TargetFullName) %in% toupper(df.mod$`Target full name`)
sum(i.mod)
i8 = which(df.mod$`Reduced model` == 1)
i8.mod = toupper(ann.soma$TargetFullName) %in% toupper(df.mod$`Target full name`[i8])
sum(i8.mod)

mod.list = list(ann.soma$TargetFullName[i.mod], ann.soma$TargetFullName[i8.mod])
names(mod.list) = c("70 genes", "8 genes")
saveRDS(mod.list, "data_generated/model_genes.rds")

index.list = list(i.mod, i8.mod, !i.mod, rep(TRUE, length(i.mod)))
names(index.list) = c("70genes", "8genes", "70genes.removed", "all.genes")
for (iname in names(index.list)) {
  cat("Saving data for", iname, "\n")
  dat2 = dat.soma[si, index.list[[iname]]] # removing model genes
  colnames(dat2) = ann.soma$TargetFullName[index.list[[iname]]]
  rownames(dat2) = info$sample

  DF = bind_cols(info %>% dplyr::select(sample, Subject, Visit, weeks), as.data.frame(dat2))
  fwrite(DF, glue::glue("data_generated/Soma_for_enet_{iname}.txt"), sep="\t")
}
