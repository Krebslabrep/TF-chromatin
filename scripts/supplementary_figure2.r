library(tidyverse)
library(ComplexHeatmap)
source("./scripts/functions/utils.r")
source("./scripts/functions/source_FootprintCharter_functions.r")

chromHMM = Load.chromHMM(GenomicTiles = TRUE)

CA_loci_df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-08-22_chromatin.influence_Sonmezer_PooledReplicates.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

# A
x = process.CA.df(x = chromatin.influence.df, cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient)

x %>%
  dplyr::select(TF.name, ChIP_annotation, CA_regulatory) %>%
  left_join(., dplyr::select(chromatin.influence.df, TF, TF.name, ChIP), by = "TF.name") %>%
  filter(!is.na(ChIP)) %>%
  group_by(TF) %>% mutate(ChIP_bin = cut(log2(ChIP+1), breaks = 50, labels = seq(2,100,2))) %>% ungroup() %>%
  mutate(ChIP_bin = factor(ChIP_bin, levels = rev(seq(2,100,2)))) %>%
  filter(ChIP_annotation == "bound") %>% 
  group_by(TF) %>% filter(n() > 30) %>% ungroup() %>%
  group_by(TF, ChIP_bin) %>%
  summarise(CA_regulatory = median(CA_regulatory, na.rm = TRUE), .groups = "drop") %>%
  mutate(TF = factor(TF, levels = c("Otx2", "Esrrb", "Yy1", "Foxd3", "Oct4", "Bach1::Mafk", "Zic3", "Sox2", "Myc", "Klf4", "Zfx", "E2f1", "Stat3", "Nrf1", "Rest", "Ctcf", "Nfya", "Banp"))) %>%
  spread(TF,  CA_regulatory, drop = TRUE) %>%
  column_to_rownames("ChIP_bin") %>%
  as.matrix() -> mat
ComplexHeatmap::Heatmap(
  matrix = mat, cluster_columns = FALSE, cluster_rows = FALSE, name = "CA frequency (%)",
  col = circlize::colorRamp2(colors = rev(colorspace::sequential_hcl(n = 100, palette = "lajolla", rev = FALSE)), breaks = seq(100)),
  row_names_side = "left", column_names_rot = 90, na_col = "white", border = TRUE, rect_gp = gpar(col = "black", lwd = 0.25), 
  width = ncol(mat)*unit(7, "mm"), height = nrow(mat)*unit(2, "mm"), row_labels = c(100, rep("", nrow(mat)-2), min(as.integer(rownames(mat)))), row_title = "ChIP-seq/-nexus percentile"
) -> pl
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs2c_", Sys.Date(), ".pdf"), width = 7, height = 4)
pl
dev.off()