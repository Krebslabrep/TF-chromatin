library(tidyverse)
source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
detach("package:plyranges")

setwd("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/")

chromatin.influence.df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-29_chromatin.influence_Sonmezer_PooledReplicates.df.qs")
chromatin.influence.df_dTAG = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-07-25_chromatin.influence_Oct4_Sox2_KD.df.qs")
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

CA_loci_df_dTAG = chromatin.influence.df_dTAG %>%
  filter(str_detect(Sample, "Sox2")) %>%
  process.CA.df(., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
    (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription"))
  )

# A
CA_loci_df %>%
  group_by(TF) %>% filter(n() > 30) %>% ungroup() %>%
  mutate(fill.variable = ifelse(TF == "unbound", "inactive", "active")) %>%
  group_by(TF) %>% mutate(m = median(CA_regulatory)) %>% ungroup() %>% arrange(m) %>% mutate(TF = factor(TF, levels = unique(TF))) %>%
  ggplot(aes(TF, CA_regulatory, fill = fill.variable)) +
  geom_violin(draw_quantiles = c(.25, .5, .75), scale = "width") +
  ggpubr::stat_compare_means(method = "wilcox", method.args = list(alternative = "greater"), label = "p.signif", size = 5, ref.group = "unbound") +
  xlab(NULL) + ylab("CA frequency") +
  scale_fill_manual(breaks = c("inactive", "active"), values = c("transparent", "grey")) +
  scale_y_continuous(breaks = c(0,100)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.1), text = element_text(size = 18), legend.title = element_blank(), legend.position = "top") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3a_", Sys.Date(), ".pdf"), pl, width = 5, height = 3)
CA_loci_df %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
    (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
    (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF)) %>% 
  group_by(TF) %>% filter(n() > 30) %>% ungroup() %>%
  group_by(TF) %>% summarise(m = median(CA_regulatory), .groups = "drop") %>% arrange(desc(m))

# B
TFBSs[filter(CA_loci_df, ChIP_annotation == "bound")$TF.name] %>%
  SingleMoleculeFootprinting::Arrange_TFBSs_clusters(., add.single.TFs = FALSE, max_cluster_width = 300, max_cluster_size = 10) -> TFBS.clusters
TFBS.cluster.compositions = unlist(TFBS.clusters$ClusterComposition)
TFBS.cluster.compositions$cluster.id = gsub("\\..*", "", names(TFBS.cluster.compositions))

# IMPROVE TEST
# Multivariate non-parametric test: 
#   https://cran.r-project.org/web/packages/MultNonParam/index.html
#   coin::independence_test()

full_join(
  CA_loci_df %>% filter(ChIP_annotation == "bound"),
  TFBS.cluster.compositions %>%
    data.frame() %>%
    dplyr::select(absolute.idx, cluster.id) %>%
    dplyr::rename("TF.name" = "absolute.idx"),
  by = "TF.name"
) %>%
  mutate(cluster.id = ifelse(is.na(cluster.id), TF.name, cluster.id)) %>%
  na.omit() %>%
  group_by(cluster.id) %>% mutate(nr.motifs = n()) %>% ungroup() %>%
  rbind(
    ., 
    CA_loci_df %>% 
      filter(ChIP_annotation == "unbound") %>%
      mutate(cluster.id = TF.name, nr.motifs = 0)
    ) %>%
  group_by(nr.motifs) %>% filter(n() > 100) %>%
  summarise(
    frequency_0.25 = quantile(CA_regulatory, 0.25, na.rm=TRUE), frequency_0.50 = quantile(CA_regulatory, 0.50, na.rm=TRUE), frequency_0.75 = quantile(CA_regulatory, 0.75, na.rm=TRUE),
    width_0.25 = quantile(width_regulatory, 0.25, na.rm=TRUE), width_0.50 = quantile(width_regulatory, 0.50, na.rm=TRUE), width_0.75 = quantile(width_regulatory, 0.75, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(nr.motifs = factor(nr.motifs, levels = seq(0,max(nr.motifs)))) %>%
  ggplot(aes(width_0.50, frequency_0.50)) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75), alpha = .5, linetype = 1) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75), alpha = .5, linetype = 1) +
  geom_point(aes(color = nr.motifs), size = 2, inherit.aes = TRUE, show.legend = FALSE) +
  geom_text(aes(label = nr.motifs), nudge_x = -5, nudge_y = 3, show.legend = FALSE) +
  scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
  scale_x_continuous(breaks = c(100,200,300), limits = c(100,300)) +
  xlab("CA width (bp)") +
  ylab("CA frequency (%)") +
  colorspace::scale_color_discrete_sequential(palette = "Emrld", rev = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3b_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5)

# C
rbind(
  full_join(
    CA_loci_df %>% filter(ChIP_annotation == "bound"),
    TFBS.cluster.compositions %>%
      data.frame() %>%
      dplyr::select(absolute.idx, cluster.id) %>%
      dplyr::rename("TF.name" = "absolute.idx"),
    by = "TF.name"
  ) %>%
    mutate(cluster.id = ifelse(is.na(cluster.id), TF.name, cluster.id)) %>%
    na.omit() %>%
    filter(TF == "Klf4") %>%
    group_by(cluster.id) %>% mutate(nr.motifs = n()) %>% ungroup() %>%
    group_by(TF, nr.motifs) %>% filter(n() > 100) %>%
    summarise(
      frequency_0.25 = quantile(CA_regulatory, 0.25, na.rm=TRUE), frequency_0.50 = quantile(CA_regulatory, 0.50, na.rm=TRUE), frequency_0.75 = quantile(CA_regulatory, 0.75, na.rm=TRUE),
      width_0.25 = quantile(width_regulatory, 0.25, na.rm=TRUE), width_0.50 = quantile(width_regulatory, 0.50, na.rm=TRUE), width_0.75 = quantile(width_regulatory, 0.75, na.rm=TRUE),
      .groups = "drop"
    ),
  full_join(
    CA_loci_df %>% filter(ChIP_annotation == "bound"),
    TFBS.cluster.compositions %>%
      data.frame() %>%
      dplyr::select(absolute.idx, cluster.id) %>%
      dplyr::rename("TF.name" = "absolute.idx"),
    by = "TF.name"
  ) %>%
    mutate(cluster.id = ifelse(is.na(cluster.id), TF.name, cluster.id)) %>%
    na.omit() %>%
    filter(TF == "Nrf1") %>%
    group_by(cluster.id) %>% mutate(nr.motifs = n()) %>% ungroup() %>%
    group_by(TF, nr.motifs) %>% filter(n() > 100) %>%
    summarise(
      frequency_0.25 = quantile(CA_regulatory, 0.25, na.rm=TRUE), frequency_0.50 = quantile(CA_regulatory, 0.50, na.rm=TRUE), frequency_0.75 = quantile(CA_regulatory, 0.75, na.rm=TRUE),
      width_0.25 = quantile(width_regulatory, 0.25, na.rm=TRUE), width_0.50 = quantile(width_regulatory, 0.50, na.rm=TRUE), width_0.75 = quantile(width_regulatory, 0.75, na.rm=TRUE),
      .groups = "drop"
    )
) %>%
  mutate(nr.motifs = factor(nr.motifs, levels = seq(max(nr.motifs)))) %>%
  ggplot(aes(width_0.50, frequency_0.50, color = TF)) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = TF), alpha = .5, linetype = 1) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = TF), alpha = .5, linetype = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = nr.motifs), nudge_x = -5, nudge_y = 3, show.legend = FALSE) +
  scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
  scale_x_continuous(breaks = c(100,200,300), limits = c(100,300)) +
  xlab("CA width (bp)") +
  ylab("CA frequency (%)") +
  scale_color_manual(values = c("grey", "blue3"), breaks = c("Klf4", "Nrf1")) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.background = element_blank(), legend.position = c(.85,.15), legend.title = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3c_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5)

# D
partition.collapsing.dictinary = split(1:12,1:12)[c(1,2,10,12,8,6,7,4,5,9,11,3)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_39881369", Sample = "SMF_MM_TKO_DE_"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE
) -> pl
pl$chromatin.influence.df %>% mutate(ChIP = 0) %>% process.CA.df(., chromHMM, ChIP_thresholds_dictionary_lenient) # 53% | 242bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2d_", Sys.Date(), ".png"), width = 25, height = 20, units = "cm", res = 300)
pl$pl
dev.off()

# E
partition.collapsing.dictinary = split(1:10,1:10)[c(1,3,8,9,2,6,10,4,5,7)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_50232513", Sample = "SMF_MM_TKO_DE_"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE
) -> pl
pl$chromatin.influence.df %>% mutate(ChIP = 0) %>% process.CA.df(., chromHMM, ChIP_thresholds_dictionary_lenient) # 68% | 265bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2e_", Sys.Date(), ".png"), width = 25, height = 20, units = "cm", res = 300)
pl$pl
dev.off()

# F
partition.collapsing.dictinary = split(1:12,1:12)[c(1,5,4,3,7,8,12,11,10,2,6,9)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_31171557", Sample = "SMF_MM_TKO_DE_"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE
) -> pl
pl$chromatin.influence.df %>% mutate(ChIP = 0) %>% process.CA.df(., chromHMM, ChIP_thresholds_dictionary_lenient) # 68% | 117bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3f_", Sys.Date(), ".png"), width = 25, height = 20, units = "cm", res = 300)
pl$pl
dev.off()

# G
partition.collapsing.dictinary = split(1:8,1:8)[c(4,5,3,6,8,7,1,2)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_53147164", Sample = "SMF_MM_TKO_DE_"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE
) -> pl
pl$chromatin.influence.df %>% mutate(ChIP = 0) %>% process.CA.df(., chromHMM, ChIP_thresholds_dictionary_lenient) # 100% | 228bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3g_", Sys.Date(), ".png"), width = 25, height = 20, units = "cm", res = 300)
pl$pl
dev.off()

# H
dplyr::intersect(
  filter(CA_loci_df_dTAG, Sample == "Sox2NT")$TF.name, 
  filter(CA_loci_df_dTAG, Sample == "Sox22h")$TF.name
) -> motifs_intersect
CA_loci_df_dTAG %<>% filter(TF.name %in% motifs_intersect)

TFBSs[unique(CA_loci_df_dTAG$TF.name)] %>%
  SingleMoleculeFootprinting::Arrange_TFBSs_clusters(., add.single.TFs = FALSE, max_cluster_width = 300, max_cluster_size = 10) -> TFBS.clusters
TFBS.cluster.compositions = unlist(TFBS.clusters$ClusterComposition)
TFBS.cluster.compositions$cluster.id = gsub("\\..*", "", names(TFBS.cluster.compositions))

full_join(
  CA_loci_df_dTAG,
  TFBS.cluster.compositions %>%
    data.frame() %>%
    dplyr::select(absolute.idx, cluster.id) %>%
    dplyr::rename("TF.name" = "absolute.idx"),
  by = "TF.name"
) %>%
  mutate(cluster.id = ifelse(is.na(cluster.id), TF.name, cluster.id)) %>%
  na.omit() %>%
  group_by(cluster.id) %>% mutate(nr.motifs = n(), nr.sox2 = sum(TF == "Sox2")) %>% ungroup() %>%
  group_by(Sample, nr.sox2, nr.motifs) %>% filter(n() > 100) %>%
  summarise(
    frequency_0.25 = quantile(CA_regulatory, 0.25, na.rm=TRUE), frequency_0.50 = quantile(CA_regulatory, 0.50, na.rm=TRUE), frequency_0.75 = quantile(CA_regulatory, 0.75, na.rm=TRUE),
    width_0.25 = quantile(width_regulatory, 0.25, na.rm=TRUE), width_0.50 = quantile(width_regulatory, 0.50, na.rm=TRUE), width_0.75 = quantile(width_regulatory, 0.75, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  # filter(
  #   (nr.sox2 == 0) |
  #   (nr.sox2 == 2 & nr.motifs %in% c(2,4)) |
  #   (nr.sox2 == 4 & nr.motifs %in% c(4))
  # ) %>%
  mutate(
    nr.motifs = factor(nr.motifs, levels = seq(max(nr.motifs))),
    nr.sox2 = factor(paste0(nr.sox2, " SOX2"), levels = paste0(c(0,2,4), " SOX2"))
    ) %>%
  ggplot(aes(width_0.50, frequency_0.50, color = Sample)) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = Sample), alpha = .5, linetype = 1) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = Sample), alpha = .5, linetype = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = nr.motifs), nudge_x = -5, nudge_y = 3, show.legend = FALSE) +
  facet_wrap(~nr.sox2) +
  scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
  scale_x_continuous(breaks = c(150,250,350), limits = c(150,350)) +
  xlab("CA width (bp)") +
  ylab("CA frequency (%)") +
  scale_color_manual(values = c("black", "salmon"), breaks = c("Sox2NT", "Sox22h")) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.background = element_blank(), legend.position = c(.5,.87), legend.title = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3h_", Sys.Date(), ".pdf"), pl, width = 9, height = 4)

