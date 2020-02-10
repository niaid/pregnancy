library(eNetXplorer)

fn = "data_generated/Clinical_data_v123.txt"
df = fread(fn)

df = df %>% 
  dplyr::select(-Subject, -Visit)

y = df$weeks
x = df %>% 
  dplyr::select(-weeks) %>% 
  tibble::column_to_rownames("sample") %>% 
  data.matrix()

dn.out = "results/Clinical_enet"
dir.create(dn.out, showWarnings = F, recursive = T)
result = eNetXplorer(x=x,y=y,family="gaussian",alpha=seq(0,1,by=0.1),
                     seed=123)#,nlambda.ext=1000,n_run=1000, n_perm_null=250)
save(result,file=file.path(dn.out, "enet_data.Rdata"))
summaryPDF(result,path=dn.out)
