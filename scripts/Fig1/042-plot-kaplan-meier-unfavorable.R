################################################################################
# Initialization
################################################################################
options(warn = -1)
options(repos = "http://cran.us.r-project.org")

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, broom, tidyverse, gridExtra)

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

df_exp_bin <-
  inner_join(df_survival, df_lr, by = "gene") %>%
  group_by(LR, cohort, id) %>%
  # high and low by median value
  mutate(exp_sum = sum(exp)) %>%
  ungroup(id) %>%
  mutate(exp_bin = if_else(exp_sum > median(exp_sum), "high", "low")) %>%
  select(LR, cohort, time, status, exp_bin) %>%
  distinct()

plot_func <- function(data, LR) {
  data %>%
    mutate(LR = LR) %>%
    ggplot(aes(time, estimate, group = strata)) +
    geom_line(aes(color = strata)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = strata), alpha = .2) +
    labs(x = "Time (days)", y = "Survival probability %") +
    ylim(0, 1) +
    scale_fill_discrete(labels = c("high", "low")) +
    scale_color_discrete(labels = c("high", "low")) +
    theme_bw() +
        theme(
        axis.ticks = element_blank(),
        axis.title = element_text(size = 8, family = "Helvetica", color = "black"),
        axis.text = element_text(size = 6, family = "Helvetica", color = "black"),
        legend.text = element_text(size = 6, family = "Helvetica", color = "black"),
        legend.title = element_text(size = 6, family = "Helvetica", color = "black"),
        strip.text = element_text(size = 6, family = "Helvetica", color = "black"),
        legend.key.size = unit(0.05, "inch"),
        panel.spacing.x = unit(0, "lines"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 1.8, 0, 1.8, "cm"),
        strip.background = element_rect(fill=NA, size=0.5),
        panel.border = element_rect(size = 0.25),
        aspect.ratio = 0.5,
        legend.margin = margin(0, 0, 0,0 , "cm"),
        legend.box.margin = margin(0, 0, 0,0 , "cm"),
    ) +
    facet_grid(LR ~ cohort, scale = "free_x")
}

df_plot <-
  df_exp_bin %>%
  group_nest(cohort, LR) %>%
  mutate(surv = map(data, ~ survfit(Surv(time, status) ~ exp_bin, data = .x))) %>%
  mutate(surv = map(surv, tidy)) %>%
  select(cohort, LR, surv) %>%
  unnest(surv) %>%
  group_nest(LR) %>%
  mutate(g = map2(data, LR, plot_func))

g <- marrangeGrob(df_plot$g, nrow = 7, ncol = 1)

################################################################################
# Output
################################################################################

ggsave("results/Fig1/kaplan_meier/unfavorable.pdf", g,
  dpi = 350, width = 210, height = 297, limitsize = FALSE, units="mm"
)
