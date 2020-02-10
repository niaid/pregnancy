library(readxl)

fn = "data/Pregnancy_clinical_data.xlsx"

dir.create("data_generated", showWarnings = F)

DF = data.frame()
for (v in 1:3) {
  vn = paste0("V",v)
  cat("Reading visit", v, "\n")
  df = read_xlsx(fn, sheet = paste("Visit", v)) %>% 
    dplyr::rename(weeks = `Weeks pregnant`) %>% 
    add_column(Visit=v, .after="Subject")
  # df = df %>% 
  #   gather("measure", "value", -c(Subject, Visit, weeks))
  cat(dim(df),"\n")
  DF = rbind(DF, df)
}

DF = DF %>% 
  add_column(sample = paste(DF$Subject, DF$Visit, sep = "_"), .before="Subject")

fwrite(DF, "data_generated/Clinical_data_v123.txt", sep="\t")
