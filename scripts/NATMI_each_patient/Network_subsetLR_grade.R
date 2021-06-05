
path_result <- "results/NATMI_each_patient/Dataframe_AtoE/Dataframe_D.csv"
path_metadata <- "data/WeiLin_pdac10/Patient_data.csv"
path_LR <- "results/Fig1/LR_adjPval_meanHR_screened.csv"
path_LR_all <- "results/Fig1/LR_HR_adjPval.csv"
path_outdir <- "results/NATMI_each_patient/Network_subsetLR_grade"

# Make output directory
if (!dir.exists(path_outdir)) {
  dir.create(path_outdir, recursive = TRUE)
}

# Load library
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, ggrepel, qgraph)

# Read data
df1 <- read_csv(path_result)

###COMMENT
# df1 <- read_csv(path_result) %>%
#         replace_na(list(HRL=0, HRH=0))
# df_test <- read_csv(path_result)
# all_equal(df1, df_test) # TRUE
###<<<<<<<

# Read metadata
dfmeta <- read_csv(path_metadata) %>%
    filter(Primary_or_Metasitasis=="Primary")

# Count the number of LR pairs
dflr <- read_csv(path_LR)
dflr %>%
    filter(meanHR<1) %>%
    nrow -> n_HRL
dflr %>% 
    filter(meanHR>1) %>%
    nrow -> n_HRH

# Count the number of all LR pairs
dflrall <- read_csv(path_LR_all)
dflrall %>% 
    nrow -> n_LR_all


#########################
# Network for each grade 
#########################
## Define edges
grades <- sort(unique(dfmeta$Grade))

for(grade in grades){
    dfmeta %>%
    filter(Grade == grade) %>%
    select(Patient_ID) %>%
    { filter(df1, Patient %in% .$Patient_ID)} %>%
    group_by(cell_type_pair) %>%
    summarise(
        mean_NormHRL = mean(NormHRL), 
        mean_NormHRH = mean(NormHRH),
        adjusted_mean_NormHRL = mean(NormHRL)/n_HRL*n_LR_all, 
        adjusted_mean_NormHRH = mean(NormHRH)/n_HRH*n_LR_all, 
    ) -> df1sub

   # Save data
    write_csv(df1, file.path(path_outdir, sprintf("adjusted_mean_enrichment_LRpairs_Grade%s.csv", grade)))

    min_value <- min(c(df1sub$adjusted_mean_NormHRL, df1sub$adjusted_mean_NormHRH)) 
    max_value <- max(c(df1sub$adjusted_mean_NormHRL, df1sub$adjusted_mean_NormHRH)) 

    # Define edges, HR<1
    df1sub %>%
        separate(cell_type_pair, sep="->", into=c("from", "to")) %>%
        mutate(thickness = adjusted_mean_NormHRL) %>%
        select(from, to, thickness) -> Edges

    pdfname <- file.path(
        path_outdir, 
        sprintf("network_cell_type_HRL_Grade%s.pdf", grade)
    )
    pdf(pdfname)
    qgraph(Edges, esize=5, 
        theme = 'gray', layout="circle", maximum=max_value,
        title = sprintf("LR pairs with HR<1, Grade %s", grade))
    dev.off()

    # Define edges, HR>1
    df1sub %>%
        separate(cell_type_pair, sep="->", into=c("from", "to")) %>%
        mutate(thickness = adjusted_mean_NormHRH) %>%
        select(from, to, thickness) -> Edges

    pdfname <- file.path(
        path_outdir, 
        sprintf("network_cell_type_HRH_Grade%s.pdf", grade)
    )
    pdf(pdfname)
    qgraph(Edges, esize=5, 
        theme = 'gray', layout="circle", maximum=max_value,
        title = sprintf("LR pairs with HR>1, Grade %s", grade))
    dev.off()

    # Define edges, HR<1 - HR>1
    df1sub %>%
        separate(cell_type_pair, sep="->", into=c("from", "to")) %>%
        mutate(thickness = adjusted_mean_NormHRL - adjusted_mean_NormHRH) %>%
        select(from, to, thickness) -> Edges

    pdfname <- file.path(
        path_outdir, 
        sprintf("network_cell_type_diff_Grade%s.pdf", grade)
    )
    pdf(pdfname)
    qgraph(Edges, esize=5, 
        theme = 'gray', layout="circle",
        title = sprintf("Difference of HR<1 from HR>1, Grade %s", grade))
    dev.off()

    # Define edges, HR<1 / HR>1
    df1sub %>%
        separate(cell_type_pair, sep="->", into=c("from", "to")) %>%
        mutate(thickness = log2(adjusted_mean_NormHRL+1e-5) - log2(adjusted_mean_NormHRH+1e-5) ) %>%
        select(from, to, thickness) -> Edges

    pdfname <- file.path(
        path_outdir, 
        sprintf("network_cell_type_fc_Grade%s.pdf", grade)
    )
    pdf(pdfname)
    qgraph(Edges, esize=5, 
        theme = 'gray', layout="circle",
        title = sprintf("log2 fold change of HR<1 from HR>1, Grade %s", grade))
    dev.off()


    # Scatter plot

    df1sub %>%
        ggplot(aes(adjusted_mean_NormHRL, adjusted_mean_NormHRH, label=cell_type_pair)) +
        geom_point(aes(color=adjusted_mean_NormHRL/adjusted_mean_NormHRH > 2)) + 
        geom_abline(intercept=0, slope=1) +
        geom_abline(intercept=0, slope=2, linetype = "dashed") +
        geom_abline(intercept=0, slope=1/2, linetype = "dashed") +
        geom_text_repel(
            data=df1sub %>% filter(adjusted_mean_NormHRL/adjusted_mean_NormHRH > 2 | mean_NormHRL/mean_NormHRH > 2),
            max.overlaps = Inf
            ) +
        theme(legend.position="bottom") + 
        lims(x=c(min_value, max_value), y=c(min_value, max_value)) +
        labs(title=sprintf("Adjusted mean enrichment of LR pairs, Grade%s", grade), 
        x="LR pairs with HR<1", y="LR pairs with HR>1") -> g1
    g1
    ggsave(file.path(path_outdir, sprintf("scatter_adjuested_mean_enrichment_HRL_vs_HRH_Grade%s.pdf", grade)), g1)
}

# sessionInfo()
sessionInfo()
