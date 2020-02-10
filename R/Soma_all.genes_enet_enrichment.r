library(fgsea)

# load data and model genes using Soma_all.genes_enet_run.r

# load eNetXploer results for all somamers
dn = "results/Soma_all.genes"
fn = "enet_data.Rdata"
load(file.path(dn, fn), verbose = T)

soma.ranked = result$feature_coef_wmean[,"a0"] %>% abs()

# load genes from previous model
mod.list = readRDS("data_generated/model_genes.rds")


fres = fgsea(mod.list, soma.ranked, nperm=500)

plotEnrichment(mod.list$`70 genes`, soma.ranked) + 
  labs(title=paste0("70 proteins, p = ", format(fres[1, "pval"], digits=3)))
dev.copy(png, "figures/Soma_70genes_enrichment_plot.png", w=600, h=300)
dev.off()

plotEnrichment(mod.list$`8 genes`, soma.ranked) + 
  labs(title=paste0("8 proteins, p = ", format(fres[2, "pval"], digits=3)))
dev.copy(png, "figures/Soma_8genes_enrichment_plot.png", w=600, h=300)
dev.off()
