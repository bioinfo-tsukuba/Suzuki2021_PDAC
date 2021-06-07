################################################################################
# Initialization
################################################################################
options(warn = -1)

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, broom, tidyverse)

################################################################################
# Input and format
################################################################################

files <- dir(path= "data/ICGC/", pattern = "*.csv.gz", full.names = TRUE)
cohorts <- files %>% str_remove_all(".*survival_|.csv.gz")
df_survival <- map2_dfr(files, cohorts, ~ read_csv(.x, col_types = cols()) %>% mutate(cohort = .y))

df_lr <- read_csv("data/NATMI_LR.csv", col_names = "LR", col_types = cols()) %>%
  filter(str_detect(LR,"SHH->SMO|SEMA4B->DCBLD2"))
output <- "results/Fig1/Fig1c.pdf"

df_survival <- df_survival %>%
  mutate(status = if_else(status == "alive", 0, 1)) # alive=0, dead=1

df_lr <- df_lr %>%
  mutate(lr_pair = LR) %>%
  separate(lr_pair, c("ligand", "receptor"), sep = "->") %>%
  pivot_longer(c(ligand, receptor), names_to = "lr", values_to = "gene")

################################################################################
# Main
################################################################################

df_plot <-
  inner_join(df_survival, df_lr, by = "gene") %>%
  group_by(LR, cohort, id) %>%
  # high and low by median value
  mutate(exp_sum = sum(exp)) %>%
  ungroup(id) %>%
  mutate(exp_bin = if_else(exp_sum > median(exp_sum), "high", "low")) %>%
  select(LR, cohort, time, status, exp_bin) %>%
  distinct() %>%
  as.data.frame()

fit <- survfit(Surv(time, status) ~ exp_bin, data = df_plot)
g <- ggsurvplot_facet(fit, df_plot, scales = "free", pval = TRUE, facet.by = c("LR", "cohort"))

################################################################################
# Output
################################################################################

ggsave(output, g, dpi = 350, width = 10, height = 5)
