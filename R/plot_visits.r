library(readxl)

fn = "data/Pregnancy_clinical_data.xlsx"
df.subj = read_xlsx(fn, sheet = "Subjects")

df.visits = df.subj %>% 
  dplyr::select(Subject)

for (v in 1:4) {
  vn = paste0("V",v)
  cat("Reading visit", v, "\n")
  df = read_xlsx(fn, sheet = paste("Visit", v))
  df.visits = df.visits %>% 
    add_column(!!(vn) := df$`Weeks pregnant`)
}

df.visits = df.visits %>% 
  gather("visit", "week", -Subject) %>% 
  mutate(Subject = fct_rev(Subject))

ggplot(df.visits, aes(week, Subject, col=visit)) +
  geom_point(shape="circle", size=2) +
  geom_point(data=df.subj, aes(`Gestational duration`, Subject), inherit.aes = F, shape="asterisk", size=1) +
  theme_minimal()
fn.fig = "figures/pregnancy_subject_weeks"
ggsave(paste0(fn.fig, ".png"), w=6,h=4)
ggsave(paste0(fn.fig, ".pdf"), w=6,h=4)

