################################################################################
# Initialization
################################################################################
options(warn = -1)

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, broom, tidyverse)

################################################################################
# Input and format
################################################################################

files <- dir(path = "data/ICGC/", pattern = "*.csv.gz", full.names = TRUE)
cohorts <- files %>% str_remove_all(".*survival_|.csv.gz")

df_survival <-
  map2_dfr(files, cohorts, ~ read_csv(.x, col_types = cols()) %>% mutate(cohort = .y)) %>%
  mutate(status = if_else(status == "alive", 0, 1)) # alive=0, dead=1

df_tmp <-
  read_csv("results/Fig1/LR_adjPval_meanHR_screened.csv", col_types = cols()) %>%
  filter(meanHR > 1)

df_lr <-
  read_csv("data/NATMI_LR.csv", col_names = "LR", col_types = cols()) %>%
  mutate(lr_pair = LR) %>%
  inner_join(df_tmp, by = "LR") %>%
  separate(lr_pair, c("ligand", "receptor"), sep = "->") %>%
  pivot_longer(c(ligand, receptor), names_to = "lr", values_to = "gene")

################################################################################
# Main
################################################################################

# report_pval <- function(data) {
#   survdiff(Surv(time, status) ~ exp_bin, data = data, rho = 1) %>%
#     glance() %>%
#     pull(p.value)
# }

# report_hr <- function(data) {
#   cox <- coxph(Surv(time, status) ~ exp_bin, data = data)
#   1 / exp(cox$coefficients)
# }

df_exp_bin <-
  inner_join(df_survival, df_lr, by = "gene") %>%
  group_by(LR, cohort, id) %>%
  # high and low by median value
  mutate(exp_sum = sum(exp)) %>%
  ungroup(id) %>%
  mutate(exp_bin = if_else(exp_sum > median(exp_sum), "high", "low")) %>%
  select(LR, cohort, time, status, exp_bin) %>%
  distinct()

df_plot <- df_exp_bin %>%
  group_nest(cohort, LR) %>%
  mutate(surv = map(data, ~ survfit(Surv(time, status) ~ exp_bin, data = .x))) %>%
  mutate(surv = map(surv, tidy)) %>%
  select(cohort, LR, surv) %>%
  unnest(surv) %>%
  mutate(strata = str_remove_all(strata, "exp_bin="))

g <-
  ggplot(df_plot, aes(time, estimate, group = strata)) +
  geom_line(aes(color = strata)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = strata), alpha = .2) +
  labs(x = "Time (days)", y = "Survival probability %") +
  ylim(0, 1) +
  theme_bw() +
  facet_wrap(~ LR + cohort, ncol = 3, scale = "free_x")

################################################################################
# Output
################################################################################

ggsave("results/Fig1/kaplan_meier/unfavorable.pdf", g, dpi = 350, width = 10, height = 150, limitsize = FALSE)
