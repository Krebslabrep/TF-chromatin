library(tidyverse)
source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
source("/g/krebs/barzaghi/Rscripts/IGV_plotting.R")

chromHMM = Load.chromHMM(GenomicTiles = TRUE)

CA_loci_df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-21_chromatin.influence_Sonmezer_SeparateReplicates.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

rest_ko_CA = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-07-15_chromatin.influence_Rest_ko.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient)

ctrl_CA = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-07-15_chromatin.influence_AMP_Barzaghi_NRF1KD_PooledReplicates.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(Sample == "amplicon__data")

CA_loci_df_NO = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-04-29_chromatin.influence_Sonmezer_NO_PooledReplicates.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

# B
partition.collapsing.dictinary = split(1:16,1:16)[c(1,9,4,10,13,5,6,2,14,15,12,3,16,11,7,8)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "AMP_56", Sample = "amplicon_DE_data"),
  rank = 1, k = 16, pool.replicates = TRUE, resize.size = 500, deduplicate = FALSE,
  partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_amplicon", remove.TFBS.labels = TRUE
) -> pl
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs1b_", Sys.Date(), ".png"), width = 32, height = 25, units = "cm", res = 300)
pl$pl
dev.off()

CallContextMethylation(
  sampleSheet = "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/Can_amplicons_NRF1KD_QuasR_input.txt", sample = "amplicon_DE_data", 
  genome = BSgenome.Mmusculus.UCSC.mm10, RegionOfInterest = IRanges::resize(Load.Sonmezer.amplicon.GRanges()["AMP_56"], 500, "center"), 
  coverage = 20, ConvRate.thr = NULL, returnSM = TRUE, clObj = NULL
)[[2]] %>%
  PlotSM(RegionOfInterest = IRanges::resize(Load.Sonmezer.amplicon.GRanges()["AMP_56"], 500, "center"), sorting.strategy = "None", SortedReads = NULL) -> sm.plot
x.axis.breaks = as.integer(c(18990931, 18991430))
sm.plot +
  facet_null() +
  ylab(paste0("2329 molecules")) +
  scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.breaks, nsmall=1, big.mark=",")) +
  scale_discrete_manual(aesthetics = "fill", values = c("black", "grey")) + 
  theme_classic() + 
  xlab("chr7") + 
  theme(text = element_text(size = 18), axis.title.y = element_text(vjust = 0), axis.title.x = element_text(vjust = 5), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_rect(linewidth = 1, fill = "transparent"),
        legend.position = "none", ) -> pl
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs1b_extra_", Sys.Date(), ".png"), width = 9.4, height = 10, units = "cm", res = 300)
pl
dev.off()

# C
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

# D
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

# E
rbind(
  rest_ko_CA %>%
    mutate(Sample = "REST KO") %>%
    filter(TF == "Rest" & ChIP_annotation == "bound"),
  ctrl_CA %>%
    filter(Sample == "amplicon__data") %>%
    mutate(Sample = "WT") %>%
    filter(TF == "Rest" & ChIP_annotation == "bound")
) %>%
  mutate(Sample = factor(Sample, levels = c("WT", "REST KO"))) -> pl.df

pl.df %>% 
  group_by(Sample) %>%
  summarise(
    width_0.25 = quantile(width, .25), frequency_0.25 = quantile(CA, 0.25), 
    width_0.50 = quantile(width, .50), frequency_0.50 = quantile(CA, 0.50), 
    width_0.75 = quantile(width, .75), frequency_0.75 = quantile(CA, 0.75), 
    .groups = "drop"
  ) %>%
  ggplot(aes(width_0.50, frequency_0.50, color = Sample)) +
  geom_point(data = pl.df, aes(width, CA, color = Sample), alpha = 1, inherit.aes = FALSE) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = Sample), alpha = 1, linetype = 1, linewidth = 1, width = 15) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = Sample), alpha = 1, linetype = 1, linewidth = 1, height = 5) +
  scale_color_manual(values = c("black", "salmon")) +
  scale_y_continuous(breaks = as.integer(c(0,50,100)), limits = c(0,100)) +
  scale_x_continuous(breaks = as.integer(c(100,200,300)), limits = c(90,375)) +
  xlab("CA width (bp)") + ylab("CA frequency (%)") +
  theme_bw() +
  theme(text = element_text(size = 18), legend.title = element_blank(), legend.background = element_blank(), legend.position = c(.8,.17), 
        legend.text.align = 1, legend.key = element_rect(fill = "transparent")) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs2e_", Sys.Date(), ".pdf"), pl, width = 5, height = 5)

# F
remove(sampleSheet, MySamples, GenomicTiles)
patch.single.site.plots(
  interpretable.master.table = data.frame(Sample = "amplicon_DE_data", TFBS.cluster = "AMP_16"), rank = 1, pool.replicates = TRUE, resize.size = 500, k = 24,
  partition.collapsing.dict = split(seq(24), seq(24))[c(18,10,6,24,17,20,4,8,19,2,11,9,12,13,3,23,7,5,1,14,15,16,21,22)],
  data.type = "WT_amplicon", remove.TFBS.labels = TRUE, deduplicate = TRUE, reutrn.chromatin.influence.df = TRUE
) -> pl
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs2f_", Sys.Date(), ".png"), width = 25, height = 18, units = "cm", res = 300)
pl$pl
dev.off()
pl.df %>% filter(TFBS.cluster == "AMP_16", Sample == "WT")

# G
remove(sampleSheet, MySamples)
patch.single.site.plots(
  interpretable.master.table = data.frame(Sample = "Rest_ko", TFBS.cluster = "AMP_16"), rank = 1, pool.replicates = TRUE, resize.size = 500, k = 24, 
  partition.collapsing.dict = split(seq(24), seq(24))[c(12,14,22,16,4,11,1,2,3,5,6,7,8,9,10,13,15,17,18,19,20,21,23,24)],
  data.type = "Rest_ko_amplicon", remove.TFBS.labels = TRUE, deduplicate = TRUE, reutrn.chromatin.influence.df = TRUE
) -> pl
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs2g_", Sys.Date(), ".png"), width = 25, height = 18, units = "cm", res = 300)
pl$pl
dev.off()
pl.df %>% filter(TFBS.cluster == "AMP_16", Sample == "REST KO")

# H-I
CA_loci_df_NO %>%
  filter((chromHMM_annotation %in% c("intergenic", "enhancer", "promoter")) | (chromHMM_annotation == "Ctcf" & TF == "Ctcf")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("intergenic", "enhancer", "promoter", "Ctcf"))) %>%
  mutate(width = ifelse(chromHMM_annotation == "intergenic", width, width_regulatory), CA = ifelse(chromHMM_annotation == "intergenic", CA, CA_regulatory)) %>%
  mutate(Sample = factor(gsub("SMF_MM_", "", Sample), levels = c("ES", "NP", "MEL", "C2C12"))) %>%
  group_by(Sample, chromHMM_annotation) %>%
  summarise(
    width_0.25 = quantile(width, .25), frequency_0.25 = quantile(CA, 0.25), 
    width_0.50 = quantile(width, .50), frequency_0.50 = quantile(CA, 0.50), 
    width_0.75 = quantile(width, .75), frequency_0.75 = quantile(CA, 0.75), 
    .groups = "drop"
  ) %>%
  ggplot(aes(width_0.50, frequency_0.50, color = chromHMM_annotation)) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = 140, color = "black", linetype = 1, linewidth = .25) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = chromHMM_annotation), alpha = 1, linetype = 1, linewidth = 1, width = 15) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = chromHMM_annotation), alpha = 1, linetype = 1, linewidth = 1, height = 5) +
  scale_color_manual(values = viridis::mako(6)[c(1,3:5)], guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = as.integer(c(0,50,100)), limits = c(0,100)) +
  scale_x_continuous(breaks = as.integer(c(100,200,300)), limits = c(75,300)) +
  xlab("CA width (bp)") + ylab("CA frequency (%)") +
  facet_wrap(~Sample) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.title = element_blank(), legend.background = element_blank(), 
        legend.text.align = 1, legend.key = element_rect(fill = "transparent")) -> pl
