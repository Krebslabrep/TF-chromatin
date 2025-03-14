library(tidyverse)
source("./scripts/functions/utils.r")
source("./scripts/functions/source_FootprintCharter_functions.r")
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

# A
CA_loci_df %>%
  group_by(TF) %>% filter(n() > 30) %>% ungroup() %>%
  mutate(fill.variable = ifelse(TF == "unbound", "unbound", "bound")) %>%
  group_by(TF) %>% mutate(m = median(CA_regulatory)) %>% ungroup() %>% arrange(m) %>% mutate(TF = factor(TF, levels = unique(TF))) %>%
  ggplot(aes(TF, CA_regulatory, fill = fill.variable)) +
  geom_boxplot() +
  xlab(NULL) + ylab("accessibility f (%)") +
  scale_fill_manual(breaks = c("unbound", "bound"), values = c("transparent", "grey")) +
  scale_y_continuous(breaks = c(0,50,100)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5), text = element_text(size = 18), legend.title = element_blank(), legend.position = "bottom") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2a_", Sys.Date(), ".pdf"), pl, width = 9, height = 5.5)
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
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = nr.motifs), linetype = 1, width = 15, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = nr.motifs), linetype = 1, height = 5, show.legend = FALSE) +
  geom_point(aes(color = nr.motifs), size = 2, inherit.aes = TRUE, show.legend = FALSE) +
  geom_text(aes(label = nr.motifs), nudge_x = -5, nudge_y = 3, show.legend = FALSE) +
  scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
  scale_x_continuous(breaks = c(0,100,200,300), limits = c(0,300)) +
  ggbreak::scale_x_break(breaks = c(5,100)) +
  xlab("average width (bp)") +
  ylab("accessibility f (%) \n (active molecules)") +
  colorspace::scale_color_discrete_sequential(palette = "Emrld", rev = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 18), axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2b_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5, onefile=FALSE)

# C
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
  group_by(nr.motifs) %>% filter(n() > 100) %>% ungroup() -> x

x %<>% filter(TF != "Ctcf")

CollectCompositeData(
  sampleSheet = "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/QuasR_input_AllCanWGpooled_dprm_DE_only.txt", 
  sample = "SMF_MM_TKO_DE_", 
  genome = BSgenome.Mmusculus.UCSC.mm10, 
  TFBSs = TFBSs[sample(filter(x, nr.motifs == 1)$TF.name, 3000)], 
  window = 1000, 
  coverage = 20, 
  ConvRate.thr = NULL, 
  cores = 16
) -> CompositeData_1_motif

CollectCompositeData(
  sampleSheet = "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/QuasR_input_AllCanWGpooled_dprm_DE_only.txt", 
  sample = "SMF_MM_TKO_DE_", 
  genome = BSgenome.Mmusculus.UCSC.mm10, 
  TFBSs = TFBS.clusters$ClusterCoordinates[unique(filter(x, nr.motifs == 2)$cluster.id)], 
  window = 1000, 
  coverage = 20, 
  ConvRate.thr = NULL, 
  cores = 16
) -> CompositeData_2_motifs

CollectCompositeData(
  sampleSheet = "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/QuasR_input_AllCanWGpooled_dprm_DE_only.txt", 
  sample = "SMF_MM_TKO_DE_", 
  genome = BSgenome.Mmusculus.UCSC.mm10, 
  TFBSs = TFBS.clusters$ClusterCoordinates[unique(filter(x, nr.motifs == 3)$cluster.id)], 
  window = 1000, 
  coverage = 20, 
  ConvRate.thr = NULL, 
  cores = 16
) -> CompositeData_3_motifs

CollectCompositeData(
  sampleSheet = "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/QuasR_input_AllCanWGpooled_dprm_DE_only.txt", 
  sample = "SMF_MM_TKO_DE_", 
  genome = BSgenome.Mmusculus.UCSC.mm10, 
  TFBSs = TFBS.clusters$ClusterCoordinates[unique(filter(x, nr.motifs == 4)$cluster.id)], 
  window = 1000, 
  coverage = 20, 
  ConvRate.thr = NULL, 
  cores = 16
) -> CompositeData_4_motifs

CollectCompositeData(
  sampleSheet = "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/QuasR_input_AllCanWGpooled_dprm_DE_only.txt", 
  sample = "SMF_MM_TKO_DE_", 
  genome = BSgenome.Mmusculus.UCSC.mm10, 
  TFBSs = TFBS.clusters$ClusterCoordinates[unique(filter(x, nr.motifs == 5)$cluster.id)], 
  window = 1000, 
  coverage = 20, 
  ConvRate.thr = NULL, 
  cores = 16
) -> CompositeData_5_motifs

rbind(
  mutate(CompositeData_1_motif, nr.motifs = 1),
  mutate(CompositeData_2_motifs, nr.motifs = 2),
  mutate(CompositeData_3_motifs, nr.motifs = 3),
  mutate(CompositeData_4_motifs, nr.motifs = 4),
  mutate(CompositeData_5_motifs, nr.motifs = 5)
) %>% distinct(TFBS_index)
  mutate(nr.motifs = as.factor(nr.motifs)) %>% 
  mutate(SMF = SMF*100) %>%
  ggplot(aes(RelStart, SMF, color = nr.motifs)) +
  geom_smooth(se = FALSE, method = "loess", span = .1) + 
  ylim(c(0, 1)) + 
  scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
  scale_x_continuous(breaks = c(-500,-250,0,250,500), limits = c(-500,500)) +
  xlab("Coord. relative to the motif cluster center") +
  ylab("SMF (1- meth %)") +
  colorspace::scale_color_discrete_sequential(palette = "Emrld", rev = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 18), axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank(), legend.position = "bottom") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2c_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5, onefile=FALSE)

# D
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
  CA_loci_df %>% 
    filter(ChIP_annotation == "unbound") %>%
    mutate(cluster.id = TF.name, nr.motifs = 0) %>%
    group_by(nr.motifs) %>% filter(n() > 100) %>%
    summarise(
      frequency_0.25 = quantile(CA_regulatory, 0.25, na.rm=TRUE), frequency_0.50 = quantile(CA_regulatory, 0.50, na.rm=TRUE), frequency_0.75 = quantile(CA_regulatory, 0.75, na.rm=TRUE),
      width_0.25 = quantile(width_regulatory, 0.25, na.rm=TRUE), width_0.50 = quantile(width_regulatory, 0.50, na.rm=TRUE), width_0.75 = quantile(width_regulatory, 0.75, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    mutate(TF = "Klf4")
  ) %>%
  mutate(nr.motifs = factor(nr.motifs, levels = seq(0,max(nr.motifs)))) %>%
  ggplot(aes(width_0.50, frequency_0.50)) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = nr.motifs), linetype = 1, width = 15, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = nr.motifs), linetype = 1, height = 5, show.legend = FALSE) +
  geom_point(aes(color = nr.motifs), size = 2, show.legend = FALSE) +
  geom_text(aes(label = nr.motifs), nudge_x = -5, nudge_y = 3, show.legend = FALSE) +
  scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
  scale_x_continuous(breaks = c(0,100,200,300), limits = c(0,300)) +
  ggbreak::scale_x_break(breaks = c(5,100)) +
  xlab("average width (bp) \n (active molecules)") +
  ylab("accessibility f (%) \n (active molecules)") +
  colorspace::scale_color_discrete_sequential(palette = "Emrld", rev = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 18), axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2d_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5, onefile=FALSE)

# E
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
    filter(TF == "Nrf1") %>%
    group_by(cluster.id) %>% mutate(nr.motifs = n()) %>% ungroup() %>%
    group_by(TF, nr.motifs) %>% filter(n() > 100) %>%
    summarise(
      frequency_0.25 = quantile(CA_regulatory, 0.25, na.rm=TRUE), frequency_0.50 = quantile(CA_regulatory, 0.50, na.rm=TRUE), frequency_0.75 = quantile(CA_regulatory, 0.75, na.rm=TRUE),
      width_0.25 = quantile(width_regulatory, 0.25, na.rm=TRUE), width_0.50 = quantile(width_regulatory, 0.50, na.rm=TRUE), width_0.75 = quantile(width_regulatory, 0.75, na.rm=TRUE),
      .groups = "drop"
    ),
  CA_loci_df %>% 
    filter(ChIP_annotation == "unbound") %>%
    mutate(cluster.id = TF.name, nr.motifs = 0) %>%
    group_by(nr.motifs) %>% filter(n() > 100) %>%
    summarise(
      frequency_0.25 = quantile(CA_regulatory, 0.25, na.rm=TRUE), frequency_0.50 = quantile(CA_regulatory, 0.50, na.rm=TRUE), frequency_0.75 = quantile(CA_regulatory, 0.75, na.rm=TRUE),
      width_0.25 = quantile(width_regulatory, 0.25, na.rm=TRUE), width_0.50 = quantile(width_regulatory, 0.50, na.rm=TRUE), width_0.75 = quantile(width_regulatory, 0.75, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    mutate(TF = "Nrf1")
) %>%
  mutate(nr.motifs = factor(nr.motifs, levels = seq(0,max(nr.motifs)))) %>%
  ggplot(aes(width_0.50, frequency_0.50)) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = nr.motifs), linetype = 1, width = 15, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = nr.motifs), linetype = 1, height = 5, show.legend = FALSE) +
  geom_point(aes(color = nr.motifs), size = 2, show.legend = FALSE) +
  geom_text(aes(label = nr.motifs), nudge_x = -5, nudge_y = 3, show.legend = FALSE) +
  scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
  scale_x_continuous(breaks = c(0,100,200,300), limits = c(0,300)) +
  ggbreak::scale_x_break(breaks = c(5,100)) +
  xlab("average width (bp) \n (active molecules)") +
  ylab("accessibility f (%) \n (active molecules)") +
  colorspace::scale_color_discrete_sequential(palette = "Emrld", rev = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 18), axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2e_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5, onefile=FALSE)

# F
partition.collapsing.dictinary = split(1:12,1:12)[c(12,1,10,6,7,4,3,2,11,5,8,9)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_11165737", Sample = "SMF_MM_TKO_DE_"), rank = 1, partition.collapsing.dict = partition.collapsing.dictinary,
  k = 12, pool.replicates = TRUE, resize.size = 500,
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE , plotting.TFBSs = TFBSs["TFBS_2201548"]
) -> pl
process.CA.df(x = mutate(pl$chromatin.influence.df[1,], ChIP = NA), cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) # 36%;23% | 178bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2f_", Sys.Date(), ".png"), width = 25, height = 20, units = "cm", res = 300)
pl$pl
dev.off()

# G
partition.collapsing.dictinary = split(1:8,1:8)[c(2,8,6,4,3,1,5,7)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_12826771", Sample = "SMF_MM_TKO_DE_"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE, 
  plotting.TFBSs = plyranges::filter(TFBS.cluster.compositions, cluster.id == plyranges::filter(TFBS.cluster.compositions, absolute.idx == "TFBS_2000388")$cluster.id)
) -> pl
pl$chromatin.influence.df[1,] %>% mutate(ChIP = 0) %>% process.CA.df(., chromHMM, ChIP_thresholds_dictionary_lenient) # 58%;38% | 198bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2g_", Sys.Date(), ".png"), width = 25, height = 20, units = "cm", res = 300)
pl$pl
dev.off()

# H
pl.TFBSs = plyranges::filter(TFBS.cluster.compositions, cluster.id == plyranges::filter(TFBS.cluster.compositions, absolute.idx == "TFBS_3355041")$cluster.id)
pl.TFBSs[1:3] = IRanges::resize(pl.TFBSs[1:3], 8, "center")
partition.collapsing.dictinary = split(1:9,1:9)[c(1,2,3,4,5,6,7)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_20302739", Sample = "SMF_MM_TKO_DE_"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE, 
  plotting.TFBSs = pl.TFBSs
) -> pl
pl$chromatin.influence.df[3,] %>% mutate(ChIP = 0) %>% process.CA.df(., chromHMM, ChIP_thresholds_dictionary_lenient) # 100%;0% | 253bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2h_", Sys.Date(), ".png"), width = 25, height = 20, units = "cm", res = 300)
pl$pl
dev.off()