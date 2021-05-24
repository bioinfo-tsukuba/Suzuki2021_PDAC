
path_result <- "results/NATMI_each_patient/Dataframe_AtoE/Dataframe_E.csv"
path_outdir <- "results/NATMI_each_patient/Network"

# Make output directory
if (!dir.exists(path_outdir)) {
  dir.create(path_outdir, recursive = TRUE)
}

# Load library
library(tidyverse)
library(qgraph)
library(ggrepel)

# Read data
df1 <- read_csv(path_result, col_names=TRUE)

# Scatter plot
min_value <- min(c(df1$mean_NormHRL, df1$mean_NormHRH)) 
max_value <- max(c(df1$mean_NormHRL, df1$mean_NormHRH)) 

df1 %>%
    ggplot(aes(mean_NormHRL, mean_NormHRH, label=cell_type_pair)) +
    geom_point(aes(color=mean_NormHRL/mean_NormHRH > 1.5)) + 
    geom_abline(intercept=0, slope=1) +
    geom_abline(intercept=0, slope=1.5, linetype = "dashed") +
    geom_abline(intercept=0, slope=1/1.5, linetype = "dashed") +
    geom_text_repel(data=df1 %>% filter(mean_NormHRH/mean_NormHRL > 1.5 | mean_NormHRL/mean_NormHRH > 1.5)) +
    theme(legend.position="bottom") + 
    lims(x=c(min_value, max_value), y=c(min_value, max_value)) +
    labs(title="Table E, Mean of normalized average expression weight", 
    x="LR pairs with HR<1", y="LR pairs with HR>1") -> g1
g1
ggsave(file.path(path_outdir, "scatter_mean_NormHRL_vs_mean_NormHRH.pdf"), g1)


# Network
## Define edges
df1 %>% 
    separate(cell_type_pair, sep="->", into=c("from", "to")) %>%
    mutate(thickness = mean_NormHRL) %>%
    select(from, to, thickness) -> Edges_HRL

df1 %>% 
    separate(cell_type_pair, sep="->", into=c("from", "to")) %>%
    mutate(thickness = mean_NormHRH) %>%
    select(from, to, thickness) -> Edges_HRH

df1 %>% 
    separate(cell_type_pair, sep="->", into=c("from", "to")) %>%
    mutate(thickness = mean_NormHRL - mean_NormHRH) %>%
    select(from, to, thickness) -> Edges_diff

# network_cell_type_HRL
pdfname <- file.path(path_outdir, "network_cell_type_HRL.pdf")
pdf(pdfname)
qgraph(Edges_HRL, esize=5, 
    theme = 'gray', layout="circle",
    title = "LR pairs with HR<1")
dev.off()

# network_cell_type_HRH
pdfname <- file.path(path_outdir, "network_cell_type_HRH.pdf")
pdf(pdfname)
qgraph(Edges_HRH, esize=5, 
    theme = 'gray', layout="circle",
    title = "LR pairs with HR>1")
dev.off()

# network_cell_type_HRL
pdfname <- file.path(path_outdir, "network_cell_type_diff.pdf")
pdf(pdfname)
qgraph(Edges_diff, esize=5, 
    theme = 'gray', layout="circle",
    title = "Difference of HR<1 from HR>1")
dev.off()

# sessionInfo()
sessionInfo()
