if (!require("styler")) install.packages("styler")


args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  styler::style_file(arg)
}
