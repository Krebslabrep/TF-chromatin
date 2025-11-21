library(tidyverse)
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/utils.r")
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/source_FootprintCharter_functions.r")
detach("package:plyranges")

chromHMM = Load.chromHMM(GenomicTiles = TRUE)
GenomicTiles = Load.GenomicTiles(tiles.width = 80)
TFBSs = Load.TFBSs()

CA_loci_df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-08-22_chromatin.influence_Sonmezer_PooledReplicates.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

CA_loci_df_dTAG = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2025-02-24_chromatin.influence_Sox2_dTAG_run2_run3_PooledReps.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

TFBSs %>%
  as.data.frame() %>%
  dplyr::select(absolute.idx, BL6.absScore) %>%
  dplyr::rename("TF.name" = "absolute.idx", "motif.score" = "BL6.absScore") -> motif.scores

# First comment
CA_loci_df %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF)) %>% 
  group_by(TF) %>% filter(n() > 30) %>% ungroup() %>%
  group_by(TF) %>% summarise(m = median(CA_regulatory), .groups = "drop") %>% 
  arrange(desc(m)) %>% dplyr::select(TF) %>% unlist(use.names = FALSE) -> tf.levels

CA_loci_df %>%
  group_by(TF) %>% filter(n() > 30) %>% ungroup() %>%
  left_join(., motif.scores, by = "TF.name") %>%
  mutate(TF = factor(TF, levels = tf.levels)) %>%
  group_by(TF) %>% mutate(scaled.motif.score = (motif.score-min(motif.score))/(max(motif.score)-min(motif.score))) %>% ungroup() %>%
  ggplot(aes(scaled.motif.score, CA_regulatory)) +
  geom_smooth(method = "loess", se = FALSE, color = "black", fill = "black") +
  facet_wrap(~TF, scales = "fixed") +
  xlab("motif score") + ylab("CA frequency (%)") +
  scale_y_continuous(limits = c(0,104), breaks = c(0,100), labels = c(0,100)) +
  scale_x_continuous(breaks = c(0,1), labels = c("min","max")) +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/motif_strenght_to_CA_freq_", Sys.Date(), ".pdf"), pl, width = 9, height = 9)

# Second comment
CA_loci_df_dTAG %>%
  filter(TF == "Sox2", ChIP_annotation == "bound") %>%
  dplyr::select(-c(chromHMM_annotation, ChIP_annotation, tot.read.count, width, CA, CA_regulatory_count, width_regulatory, -TFBS.cluster, -TF)) %>%
  left_join(., motif.scores, by = "TF.name") -> tmp

tmp %>%
  mutate(motif.score.bin = cut(motif.score, 10, include.lowest = TRUE)) %>%
  mutate(condition = ifelse(Sample == "Sox2_NT", "untreated", "2h dTAG")) %>%
  ggplot(aes(motif.score.bin, CA_regulatory, fill = condition)) +
  geom_boxplot(alpha = .5) +
  xlab("Sox2 motif score bin") + ylab("CA frequency (%)") +
  scale_fill_manual(values = c("black", "salmon"), breaks = c("untreated", "2h dTAG")) +
  scale_y_continuous(limits = c(0,100), breaks = c(0,100), labels = c(0,100)) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.title = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/sox2_motif_streght_bins_", Sys.Date(), ".pdf"), pl, width = 9, height = 6)












