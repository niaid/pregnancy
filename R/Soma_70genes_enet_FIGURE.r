library(eNetXplorer)

dn.enet = "results/Soma_70genes/"
dn.fig = str_replace(dn.enet, "results", "figures")
dir.create(dn.fig, showWarnings = F, recursive = T)
fn.run = file.path(dn.enet, "enet_data.Rdata")

load(fn.run, verbose = T)

a = 0
ai = which.min(abs(result$alpha - a))
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

ggplot(DF.plot, aes(measured, pred_mean, col=visit)) +
  geom_point() +
  geom_errorbar(aes(ymin=pred_mean-pred_sd, ymax=pred_mean+pred_sd), width=0.1) +
  geom_smooth(method="lm", alpha=0.2, col="red", lty=2) +
  # geom_text(aes(label=sub("s","",subject)), hjust=-0.3, size=3) +
  geom_text(data=df.cor, aes(x,y,label=label), hjust=1.1, vjust=-1, inherit.aes = F) +
  xlab("Actual GA") +
  ylab("Out-of-bag Predicted GA") +
  theme_bw()
fn.fig = file.path(dn.fig, sprintf("eNet_%s_scatter", albl))
ggsave(file=paste0(fn.fig, ".png"), w=6, h=5)
ggsave(file=paste0(fn.fig, ".pdf"), w=6, h=5, useDingbats=FALSE)

