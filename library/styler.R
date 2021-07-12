if (!require("styler")) install.packages("styler")

rscripts <- system("find . -type f | grep -i .R$", intern = TRUE)

parallel::mclapply(
  rscripts,
  function(x) styler::style_file(x),
  mc.cores = parallel::detectCores() - 1
)
