library(readxl)
library(ComplexHeatmap)
library(circlize)

# Flow data ---------------------------------------------------------------

fn.flow = "data/Pregnancy_flow_cytometry.xlsx"
df.flow = read_xlsx(fn.flow, sheet = "frequencies") %>% 
  dplyr::rename(Visit = TimePoint)
info.flow = df.flow[,1:2] %>% 
  mutate(sample = paste(Subject, Visit, sep="_"))
dat.flow = df.flow[,-(1:2)] %>% 
  data.matrix()
rownames(dat.flow) = info.flow$sample
ids = colnames(dat.flow)[-1]
ids = round(as.numeric(ids), digits = 1) %>% as.character()
colnames(dat.flow)[-1] = ids
colnames(df.flow)[-(1:3)] = ids

ann.flow = read_xlsx(fn.flow, sheet = "all gates") %>% 
  dplyr::select(ID = `Pop Code`, Name = `Subset name`)
ia = match(colnames(dat.flow), ann.flow$ID)
ann.flow = ann.flow[ia,]
ann.flow[is.na(ann.flow$ID),] = c("LYM","Lymphocytes")


dir.create("data_generated", showWarnings = F)
save(info.flow, dat.flow, ann.flow, file="data_generated/Pregnancy_flow_cytometry.RData")


# Compare cell population frequencies between visits ----------------------

fc31 = c()
fc41 = c()
fc43 = c()
wt31 = list()
wt41 = list()
wt43 = list()
for (ci in colnames(dat.flow)) {
  fc31[ci] = mean(dat.flow[info.flow$Visit==3,ci], na.rm=T) / mean(dat.flow[info.flow$Visit==1,ci], na.rm=T)
  fc41[ci] = mean(dat.flow[info.flow$Visit==4,ci], na.rm=T) / mean(dat.flow[info.flow$Visit==1,ci], na.rm=T)
  fc43[ci] = mean(dat.flow[info.flow$Visit==4,ci], na.rm=T) / mean(dat.flow[info.flow$Visit==3,ci], na.rm=T)
  wt31[[ci]] = wilcox.test(dat.flow[info.flow$Visit==1,ci], dat.flow[info.flow$Visit==3,ci], paired=T, exact = F)
  wt41[[ci]] = wilcox.test(dat.flow[info.flow$Visit==1 & info.flow$Subject!=8,ci], dat.flow[info.flow$Visit==4,ci], paired=T, exact = F)
  wt43[[ci]] = wilcox.test(dat.flow[info.flow$Visit==3 & info.flow$Subject!=8,ci], dat.flow[info.flow$Visit==4,ci], paired=T, exact = F)
}

wt31.p = map_dbl(wt31, function(x) x$p.value)
wt31.padj = p.adjust(wt31.p, method = "BH")
i.31 = names(wt31)[wt31.padj < 0.05]

wt41.p = map_dbl(wt41, function(x) x$p.value)
wt41.padj = p.adjust(wt41.p, method = "BH")
i.41 = names(wt41)[wt41.padj < 0.05]

wt43.p = map_dbl(wt43, function(x) x$p.value)
wt43.padj = p.adjust(wt43.p, method = "BH")
i.43 = names(wt43)[wt43.padj < 0.05]

i.flow = unique(c(i.31, i.41, i.43))
i.flow = i.flow[order(as.numeric(i.flow), na.last = T)]

# save test results to file

F.df.all = data.frame(ID = colnames(dat.flow), FC.3.v.1 = fc31, FC.4.v.3 = fc43, FC.4.v.1 = fc41,
                  wP.3.v.1 = wt31.p, wP.4.v.3 = wt43.p, wP.4.v.1 = wt41.p,
                  wQ.3.v.1 = wt31.padj, wQ.4.v.3 = wt43.padj, wQ.4.v.1 = wt41.padj)
F.df.all = left_join(F.df.all, ann.flow, by="ID")
fwrite(F.df.all, "results/flow_visits_comparison_all.txt", sep="\t")

F.df = F.df.all[match(i.flow, F.df.all$ID),]
fwrite(F.df, "results/flow_visits_comparison.txt", sep="\t")

F.df.V31 = data.frame(ID = colnames(dat.flow), FC = fc31, wP = wt31.p, wQ = wt31.padj) %>% 
  dplyr::filter(wQ <= 0.05) %>% 
  arrange(wP)
fwrite(F.df.V31, "results/flow_visits.3.v.1_comparison.txt", sep="\t")


# combine the data and filter by fold change

DF = data.frame(ID = i.flow, V31 = fc31[i.flow], V43 = fc43[i.flow], V41 = fc41[i.flow])
X = DF %>% inner_join(ann.flow, by="ID") %>% 
  mutate(Name = glue::glue("{Name} [{ID}]")) %>% 
  dplyr::select(-ID) %>% 
  tibble::column_to_rownames("Name") %>% 
  data.matrix()

X.log2 = log2(X)
X.max = apply(X.log2, 1, function(x)max(abs(x),na.rm=T))  
fc.th = 1.2
i.sel = X.max > log2(fc.th)
iord = order(X.max[i.sel], decreasing = T)
XX = X.log2[i.sel, ]#[i.ord,]
range(XX)

P.df = data.frame(ID = i.flow, V31 = wt31.padj[i.flow], V43 = wt43.padj[i.flow], V41 = wt41.padj[i.flow])
P = P.df %>% inner_join(ann.flow %>% dplyr::select(ID, Name), by="ID") %>%
  mutate(Name = glue::glue("{Name} [{ID}]")) %>% 
  tibble::column_to_rownames("Name") %>% 
  dplyr::select(-ID) %>% 
  data.matrix()
PP = P[i.sel, ]


# additional annotation of selected gates (with BH-adjusted Wilcoxon p-value < 0.05 and FC > 1.2)

ann.flow.selected = read_xlsx(fn.flow, sheet = "selected gates") %>% 
  mutate(Name2 = glue::glue("{Name} [{ID}]"))
ro = match(ann.flow.selected$ID, i.flow[i.sel])
split.vec = ann.flow.selected$group %>% fct_inorder()

# prepare data for heatmap

max.val = max(abs(X.log2))
fnt.sz = 12

colnames(XX) = c("V3 / V1", "V4 / V3", "V4 / V1")
XX = XX[ro,]
PP = PP[ro,]

rownames(XX) = ann.flow.selected$Name2
hm = Heatmap(XX, name = "Fold Change, log2", cluster_columns = F, cluster_rows = F,
             split = split.vec, 
             gap = unit(2,"mm"),
        row_names_max_width = max_text_width( rownames(XX), gp = gpar(fontsize = 12) ),
        column_names_side = "top",
        col = colorRamp2(c(-1.4,0,1.4), c("blue", "white", "red")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(PP[i, j] < 0.001) {
            grid.text("***", x, y, gp = gpar(fontsize = fnt.sz))
          } else if(PP[i, j] < 0.01) {
            grid.text("**", x, y, gp = gpar(fontsize = fnt.sz))
          } else if(PP[i, j] < 0.05) {
            grid.text("*", x, y, gp = gpar(fontsize = fnt.sz))
          }
        },
        heatmap_legend_param = list(color_bar = "continuous", 
                                    legend_width = unit(3,"cm"),
                                    legend_direction = "horizontal"))

fn.hm = glue::glue("figures/Flow_FC_heatmap_fc.th.{fc.th}")
png(paste0(fn.hm, ".png"), w=450, h=700)
draw(hm, heatmap_legend_side = "bottom")
dev.off()
pdf(paste0(fn.hm, ".pdf"), w=6, h=10)
draw(hm, heatmap_legend_side = "bottom")
dev.off()


# plot individual samples -------------------------------------------------

pop.use = i.flow[i.sel]
pop.name = paste(ann.flow.selected$subset, ann.flow.selected$ID, sep=", population ")[match(pop.use, ann.flow.selected$ID)]
DF2 = df.flow %>%
  dplyr::select(Subject, Visit, one_of(pop.use)) %>%
  gather("ID","freq", -Subject, -Visit) %>%
  arrange(Visit, Subject) %>%
  mutate(Sample = paste(Subject, Visit, sep="_") %>% fct_inorder()) %>% 
  mutate(Subject = fct_inorder(as.character(Subject))) %>%
  mutate(ID = factor(ID, levels = pop.use)) %>% 
  group_by(ID, Subject) %>% 
  mutate(FC = freq / freq[Visit==1]) %>% 
  ungroup()

dir.create("figures/flow_profiles", showWarnings = F)

for(k in seq_along(pop.use)) {
  cat(pop.use[k]," ")
  ggplot(DF2 %>% dplyr::filter(ID==pop.use[k]), aes(Visit, FC, group=Subject, col=Subject)) +
    geom_path(show.legend = F) +
    geom_hline(yintercept = 1, col="black") +
    # scale_y_continuous(trans = "log2") +
    xlab("Study visit") +
    ylab("Fold change in percent\nof parent population") +
    ggtitle(pop.name[k]) +
    theme_bw()
  ggsave(glue::glue("figures/flow_profiles/flow_profile_{pop.use[k]}.png"), w=4,h=3)
}

