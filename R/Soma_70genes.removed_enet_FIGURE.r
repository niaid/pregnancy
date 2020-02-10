library(eNetXplorer)

dn.enet = "results/Soma_70genes.removed"
dn.fig = str_replace(dn.enet, "results", "figures")
dir.create(dn.fig, showWarnings = F, recursive = T)
fn.run = file.path(dn.enet, "enet_data.Rdata")

load(fn.run, verbose = T)

a = 0.2
ai = which(result$alpha == a)
albl = paste0("a",a)

feature = result$feature

freq.mean = result$feature_freq_mean[,ai]
freq.sd = result$feature_freq_sd[,ai]
freq.pval = result$feature_freq_model_vs_null_pval[,ai]
coef.mean = result$feature_coef_wmean[,ai]
coef.sd = result$feature_coef_wsd[,ai]
coef.pval = result$feature_coef_model_vs_null_pval[,ai]
DF.model = data.frame(feature, freq.mean, freq.sd, freq.pval, coef.mean, coef.sd, coef.pval, model="model")

freq.mean = result$null_feature_freq_mean[,ai]
freq.sd = result$null_feature_freq_sd[,ai]
coef.mean = result$null_feature_coef_wmean[,ai]
coef.sd = result$null_feature_coef_wsd[,ai]
DF.null = data.frame(feature, freq.mean, freq.sd, freq.pval, coef.mean, coef.sd, coef.pval, model="null")

DF.stat = rbind(DF.model, DF.null) %>% 
  mutate(model = factor(model, levels=c("null","model"))) %>% 
  mutate(feature = factor(feature, levels=sort(unique(feature), decreasing=T)))

# scatter plot of predictin performance
DF.plot = data.frame(
  sample = rownames(result$predictor),
  measured = result$response,
  pred_mean = result$predicted_values[[ai]][,1],
  pred_sd = result$predicted_values[[ai]][,2]
) %>% 
  separate(sample, c("subject","visit"), sep="_", remove=F)

df.cor = data.frame(x=Inf, y=-Inf, 
                    label=sprintf("Pearson correlation: %.2g (p = %.3g)", 
                                  result$model_QF_est[albl], result$QF_model_vs_null_pval[albl]))

p0 = ggplot(DF.plot, aes(measured, pred_mean, col=visit)) +
  geom_point() +
  geom_errorbar(aes(ymin=pred_mean-pred_sd, ymax=pred_mean+pred_sd), width=0.1) +
  geom_smooth(method="lm", alpha=0.2, col="red", lty=2) +
  # geom_text(aes(label=sub("s","",subject)), hjust=-0.3, size=3) +
  geom_text(data=df.cor, aes(x,y,label=label), hjust=1.1, vjust=-1, inherit.aes = F) +
  xlab("Actual GA") +
  ylab("Out-of-bag Predicted GA") +
  theme_bw()
p0

features.in = NULL
if(is.null(features.in)){
  fi = DF.model$freq.mean > max(DF.null$freq.mean, na.rm=T) & 
    abs(DF.model$coef.mean) > max(abs(DF.null$coef.mean[!is.na(DF.model$coef.mean)]), na.rm=T) & 
    DF.model$coef.pval < 0.2
} else {
  fi = DF.stat$feature %in% features.in
}
features.in = feature[fi]
DF.plot.2 = DF.stat %>% 
  dplyr::filter(feature %in% features.in) %>% 
  arrange(desc(model), abs(coef.mean)) %>% 
  mutate(feature = fct_inorder(as.character(feature)))
p1 = ggplot(DF.plot.2, aes(feature, freq.mean*100, fill=model)) +
  geom_col(position="dodge") +
  geom_errorbar(aes(ymin=pmax(0, (freq.mean-freq.sd)*100), ymax=pmin(100,(freq.mean+freq.sd)*100)),
                position=position_dodge(width = 0.9), width=0.25) +
  scale_fill_manual(values = c("grey60","red")) +
  scale_y_continuous(breaks=c(0, 50, 100)) +
  coord_flip() +
  xlab("") +
  ylab("% selected") +
  guides(fill=F) +
  theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
p2 = ggplot(DF.plot.2, aes(feature, coef.mean, fill=model)) +
  geom_col(position="dodge") +
  geom_errorbar(aes(ymin=coef.mean-coef.sd, ymax=coef.mean+coef.sd),
                position=position_dodge(width = 0.9), width=0.25) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("grey60","red")) +
  # coord_flip(ylim=c(-1.1,1.1)) +
  coord_flip() +
  xlab("") +
  ylab("Coefficient mean") +
  guides(fill=F) +
  theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
g0 <- ggplot_gtable(ggplot_build(p0))
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))
l = list(g0, g1, g2)

library(gridExtra)
lay  = matrix(c(1,1,1,1,2,2,2,3,3), nrow=1)
mg = grid.arrange(arrangeGrob(g0,g1,g2, layout_matrix = lay))#, width=c(0.5, 0.5))
fn.fig = file.path(dn.fig, sprintf("eNet_%s_selected", albl))
ggsave(file=paste0(fn.fig, ".png"), plot=mg, w=12, h=5)
ggsave(file=paste0(fn.fig, ".pdf"), plot=mg, w=12, h=5, useDingbats=FALSE)
