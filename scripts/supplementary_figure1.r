library(tidyverse)
source("./scripts/functions/utils.r")
source("./scripts/functions/source_FootprintCharter_functions.r")
source("./scripts/functions/IGV_plotting.r")

chromHMM = Load.chromHMM(GenomicTiles = TRUE)

CA_loci_df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-21_chromatin.influence_Sonmezer_SeparateReplicates.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

# B
CA_loci_df %>%
  dplyr::select(Sample, TF.name, CA) %>%
  filter(!is.na(Sample)) %>%
  spread(Sample, CA) %>%
  na.omit() %>%
  ggplot(aes(SMF_MM_TKO_R1, SMF_MM_TKO_R2)) +
  ggpointdensity::geom_pointdensity(adjust = 3) +
  ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", size = 4.5) +
  viridis::scale_color_viridis() +
  scale_x_continuous(breaks = c(0,100), limits = c(0,100)) + 
  scale_y_continuous(breaks = c(0,100), limits = c(0,100)) + 
  xlab("replicate 1") +
  ylab("replicate 2") +
  ggtitle("frequency") +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs1c_", Sys.Date(), ".png"), pl, width = 5.5, height = 5.5, dpi = 300)

# C
CA_loci_df %>%
  dplyr::select(Sample, TF.name, width) %>%
  filter(!is.na(Sample)) %>%
  spread(Sample, width) %>%
  na.omit() %>%
  ggplot(aes(SMF_MM_TKO_R1, SMF_MM_TKO_R2)) +
  ggpointdensity::geom_pointdensity(adjust = 3) +
  ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", size = 4.5) +
  viridis::scale_color_viridis() +
  scale_x_continuous(breaks = c(0,100,500), limits = c(0,500)) +
  scale_y_continuous(breaks = c(0,100,500), limits = c(0,500)) +
  xlab("replicate 1") +
  ylab("replicate 2") +
  ggtitle("width") +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs1d_", Sys.Date(), ".png"), pl, width = 5.5, height = 5.5, dpi = 300)

# D
CA_loci_df %>% 
  filter((chromHMM_annotation %in% c("enhancer", "promoter")) | (chromHMM_annotation == "Ctcf" & TF == "Ctcf")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("intergenic", "enhancer", "promoter", "Ctcf"))) %>%
  mutate(linker = CA - CA_regulatory, accessible = CA_regulatory) %>%
  dplyr::select(TF.name, chromHMM_annotation, linker, accessible) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("Ctcf", "promoter", "enhancer"))) %>%
  ggplot(aes(linker, color = chromHMM_annotation)) +
  geom_density(linewidth = 1) +
  scale_color_manual(values = viridis::mako(6)[5:3], guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(breaks = as.integer(c(0,50,100)), limits = c(0,100)) +
  xlab("linker frequency (%)") +
  theme_bw() +
  theme(
    text = element_text(size = 18), legend.title = element_blank(), 
    legend.position = c(.77,.87), legend.background = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs1d_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5)