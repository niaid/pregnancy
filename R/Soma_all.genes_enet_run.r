library(eNetXplorer)

fn = "data_generated/Soma_for_enet_all.genes.txt"
df = fread(fn)

df = df %>% 
  dplyr::select(-Subject, -Visit)

y = df$weeks
x = df %>% 
  dplyr::select(-weeks) %>% 
  tibble::column_to_rownames("sample") %>% 
  data.matrix()

dn.out = "results/Soma_all.genes"
dir.create(dn.out, showWarnings = F, recursive = T)
result = eNetXplorer(x=x,y=y,family="gaussian",alpha=seq(0,0.1,by=0.1),
                     seed=123)#,nlambda.ext=1000,n_run=1000, n_perm_null=250)
save(result,file=file.path(dn.out, "enet_data.Rdata"))
summaryPDF(result,path=dn.out)
