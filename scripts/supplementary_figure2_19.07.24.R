library(tidyverse)
library(ComplexHeatmap)
source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")

setwd("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/")

chromatin.influence.df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-29_chromatin.influence_Sonmezer_PooledReplicates.df.qs")
chromHMM = Load.chromHMM(GenomicTiles = TRUE)
GenomicTiles = Load.GenomicTiles(tiles.width = 80)
TFBSs = Load.TFBSs()

CA_loci_df = process.CA.df(x = chromatin.influence.df, cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

# C
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

# D
TFBSs[filter(CA_loci_df, TF != "unbound")$TF.name] %>%
  SingleMoleculeFootprinting::Arrange_TFBSs_clusters(., add.single.TFs = FALSE, max_cluster_width = 300) -> TFBS.clusters
TFBS.cluster.compositions = unlist(TFBS.clusters$ClusterComposition)
TFBS.cluster.compositions$cluster.id = gsub("\\..*", "", names(TFBS.cluster.compositions))

# left_join(
#   CA_loci_df,
#   TFBS.cluster.compositions %>%
#     data.frame() %>%
#     dplyr::select(absolute.idx, cluster.id) %>%
#     dplyr::rename("TF.name" = "absolute.idx"),
#   by = "TF.name"
# ) %>%
#   na.omit() %>%
#   group_by(cluster.id, TF) %>%
#   mutate(nr.homotypic.cobinders = sum(TF == TF)) %>%
#   ungroup() %>%
#   group_by(TF, nr.homotypic.cobinders) %>% filter(n() > 20) %>% ungroup() %>%
#   ggplot(aes(as.character(nr.homotypic.cobinders), CA_regulatory)) +
#   geom_violin(draw_quantiles = c(.25,.5,.75), fill = "grey") +
#   facet_wrap(~TF) +
#   theme_bw()

left_join(
  CA_loci_df,
  TFBS.cluster.compositions %>%
    data.frame() %>%
    dplyr::select(absolute.idx, cluster.id) %>%
    dplyr::rename("TF.name" = "absolute.idx"),
  by = "TF.name"
) %>%
  na.omit() %>%
  group_by(cluster.id, TF) %>%
  mutate(nr.homotypic.cobinders = sum(TF == TF)) %>%
  ungroup() %>%
  group_by(TF, nr.homotypic.cobinders) %>% filter(n() > 20) %>% ungroup() %>%
  group_by(TF) %>% filter(all(c(1, 2) %in% unique(nr.homotypic.cobinders))) %>%
  group_by(TF, nr.homotypic.cobinders) %>% summarise(CA_regulatory = median(CA_regulatory), .groups = "drop") %>%
  spread(nr.homotypic.cobinders, CA_regulatory) %>%
  column_to_rownames("TF") %>%
  as.matrix() -> medians.mat

left_join(
  CA_loci_df,
  TFBS.cluster.compositions %>%
    data.frame() %>%
    dplyr::select(absolute.idx, cluster.id) %>%
    dplyr::rename("TF.name" = "absolute.idx"),
  by = "TF.name"
) %>%
  na.omit() %>%
  group_by(cluster.id, TF) %>%
  mutate(nr.homotypic.cobinders = sum(TF == TF)) %>%
  ungroup() %>%
  group_by(TF, nr.homotypic.cobinders) %>% filter(n() > 20) %>% ungroup() %>%
  group_by(TF) %>% filter(all(c("1", "2") %in% unique(nr.homotypic.cobinders))) %>%
  summarise(
    `1` = 1, 
    `2` = tryCatch({wilcox.test(x = ifelse(nr.homotypic.cobinders == "2", CA_regulatory, NA), y = ifelse(nr.homotypic.cobinders == "1", CA_regulatory, NA), alternative = "greater")$p.value}, error = function(e){NA}), 
    `3` = tryCatch({wilcox.test(x = ifelse(nr.homotypic.cobinders == "3", CA_regulatory, NA), y = ifelse(nr.homotypic.cobinders == "2", CA_regulatory, NA), alternative = "greater")$p.value}, error = function(e){NA}),
    `4` = tryCatch({wilcox.test(x = ifelse(nr.homotypic.cobinders == "4", CA_regulatory, NA), y = ifelse(nr.homotypic.cobinders == "3", CA_regulatory, NA), alternative = "greater")$p.value}, error = function(e){NA}),
    `5` = tryCatch({wilcox.test(x = ifelse(nr.homotypic.cobinders == "5", CA_regulatory, NA), y = ifelse(nr.homotypic.cobinders == "4", CA_regulatory, NA), alternative = "greater")$p.value}, error = function(e){NA}),
    `6` = tryCatch({wilcox.test(x = ifelse(nr.homotypic.cobinders == "6", CA_regulatory, NA), y = ifelse(nr.homotypic.cobinders == "5", CA_regulatory, NA), alternative = "greater")$p.value}, error = function(e){NA}), .groups = "drop"
  ) %>%
  column_to_rownames("TF") %>%
  as.matrix() -> pvalues.mat
pvalues.mat.stars = pvalues.mat
pvalues.mat.stars[pvalues.mat.stars<.05] = "*"
pvalues.mat.stars[pvalues.mat.stars>=.05 | is.na(pvalues.mat.stars)] = ""
Heatmap(
  matrix = medians.mat, cluster_columns = FALSE, show_row_dend = FALSE, name = "CA frequency (%)",
  rect_gp = gpar(col = "black", lwd = 0.25), col = circlize::colorRamp2(colors = colorspace::sequential_hcl(n = 100, palette = "Viridis", rev = FALSE), breaks = seq(100)),
  row_names_side = "left", column_names_rot = 0, na_col = "white", border = TRUE, column_names_centered = TRUE,
  width = ncol(medians.mat)*unit(6, "mm"), height = nrow(medians.mat)*unit(6, "mm"),
  cell_fun = function(j, i, x, y, w, h, col) {grid.text(pvalues.mat.stars[i, j], x, y, vjust = .8, gp = gpar(col = "white", fontsize = 20))}
) -> pl
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs2b_", Sys.Date(), ".pdf"), width = 4, height = 5)
pl
dev.off()
