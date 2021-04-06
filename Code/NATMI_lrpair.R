#Preprocess the NATMI LR pair dataset
library(dplyr)
lrdbs <- read.csv("./NATMI/lrdbs/lrc2p.csv")
head(lrdbs)

NATMI_lrdbs <- lrdbs %>% as.data.frame() %>% mutate(ligandreceptor_pair = paste(!!!rlang::syms(c("Ligand.gene.symbol", "Receptor.gene.symbol")), sep="->"))

NATMI_lrpair = select(.data = NATMI_lrdbs, 3)

head(NATMI_lrpair)

write.csv(NATMI_lrpair, "./data/NATMI_lrpair")