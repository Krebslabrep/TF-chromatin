library(parallel)

lapply(
  list.files("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/FootprintCharter_functions/", ".R", full.names = TRUE),
  source
)