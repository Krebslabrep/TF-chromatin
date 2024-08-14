library(tidyverse)
source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
source("/g/krebs/barzaghi/Rscripts/IGV_plotting.R")

setwd("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/")

chromatin.influence.df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-29_chromatin.influence_Sonmezer_PooledReplicates.df.qs")
chromatin.influence.df_dTAG = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-07-25_chromatin.influence_Oct4_Sox2_KD.df.qs")
chromatin.influence.df_rest.ko = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-07-15_chromatin.influence_Rest_ko.qs")
chromatin.influence.df_TKO.AMP = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-07-15_chromatin.influence_AMP_Barzaghi_NRF1KD_PooledReplicates.df.qs")

chromHMM = Load.chromHMM(GenomicTiles = TRUE)
GenomicTiles = Load.GenomicTiles(tiles.width = 80)
TFBSs = Load.TFBSs()

CA_loci_df = process.CA.df(chromatin.influence.df, cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient)
CA_loci_df_dTAG = process.CA.df(chromatin.influence.df_dTAG, cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient)
rest_ko_CA = process.CA.df(chromatin.influence.df_rest.ko, cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient)
ctrl_CA = process.CA.df(chromatin.influence.df_TKO.AMP, cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient)

# A
CA_loci_df %>%
  filter(chromHMM_annotation == "intergenic" & ChIP_annotation == "unbound") %>%
  mutate(m = median(CA), r = abs(CA-m)) %>%
  filter(tot.read.count >= 100  & CA <= m+5 & CA >= m-5) %>%
  filter(TF == "Stat3") %>%
  # filter(!TF %in% c("Bach1::Mafk", "Smad2::Smad3::Smad4", "Stat3", "Zfx")) %>%
  # filter(TFBS.cluster == "GenomicTile_23340752") %>% # This is the one I have in the figure. But it's bound
  # filter(TFBS.cluster == "GenomicTile_29037971") %>% # 10
  mutate(Sample = "SMF_MM_TKO_DE_") -> intergenic.single.sites

# GenomicTile_29037953

# TFs (left to right): Otx2, Myc, Oct4, Yy1
# 24: 3,8,9,2,10,7,1,4,5,6,11,12
# 39: 7,9,3,1,2,4,5,6,8,10,11,12
# 40: 
partition.collapsing.dictinary = split(1:12,1:12)[c(3,6,5,1,2,7)]
i=40
patch.single.site.plots(
  interpretable.master.table = intergenic.single.sites, rank = i, k = 12, 
  plotting.TFBSs = plyranges::filter_by_overlaps(TFBSs, IRanges::resize(GenomicTiles[intergenic.single.sites$TFBS.cluster[i]], 500, "center")),
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = FALSE
) -> pl
mutate(pl$chromatin.influence.df, CA = (accessible.read.count/tot.read.count)*100)$CA # 30%
pl$chromatin.influence.df$acc.width.distro # 92bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1b_", Sys.Date(), ".png"), width = 27, height = 20, units = "cm", res = 300)
pl$pl
dev.off()

RegionOfInterest = IRanges::resize(GenomicTiles["GenomicTile_8311281"], 9000, "center")

plot_genomic_track(
  sampleSheet = "/g/krebs/barzaghi/HTS/DHS/Qinput.txt",
  samples = "DNAse-seq_mESC_TKO",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 50,
  max.y.lim = 400,
  color = "black",
  y.labs = "DNAse"
) -> DHS_track

plot_genomic_track(
  sampleSheet = ChIP_data_dictionary$Otx2,
  samples = "Otx2",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 10,
  max.y.lim = 100,
  color = "black", 
  y.labs = "Oct4"
) -> Oct4_track

plot_genomic_track(
  sampleSheet = ChIP_data_dictionary$`Bach1::Mafk`,
  samples = "Bach1::Mafk",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 10,
  max.y.lim = 100,
  color = "black",
  y.labs = "Sox2"
) -> Sox2_track

plot_genomic_track(
  sampleSheet = "/g/krebs/barzaghi/HTS/MNase/Qinput.txt",
  samples = "MNase_mESC_129S6/SvEvTac_WT",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  max.y.lim = 160,
  color = "black", 
  y.labs = "MNAse",
  plot.coordinates = TRUE
) -> MNase_track

p_final <- DHS_track + Oct4_track + Sox2_track + MNase_track + 
  plot_layout(ncol = 1, heights = c(1/4, 1/4, 1/4, 1/4))

pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1b_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()

# B
CA_loci_df %>% 
  filter((chromHMM_annotation == "intergenic" & ChIP_annotation == "unbound") | (chromHMM_annotation == "Ctcf" & TF == "Ctcf" & ChIP_annotation == "bound")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("intergenic", "Ctcf"))) %>%
  group_by(chromHMM_annotation) %>% 
  summarise(median = median(CA), .groups = "drop") %>%
  readr::write_delim(., col_names = TRUE, file = paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig~2b_values_", Sys.Date(), ".tsv"), delim = "\t", append = FALSE)
CA_loci_df %>% 
  filter((chromHMM_annotation == "intergenic" & ChIP_annotation == "unbound") | (chromHMM_annotation == "Ctcf" & TF == "Ctcf" & ChIP_annotation == "bound")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("intergenic", "Ctcf"))) %>%
  ggplot(aes(chromHMM_annotation, CA)) +
  geom_violin(aes(fill = chromHMM_annotation), color = "black", linewidth = 1, alpha = 0.75, draw_quantiles = c(.25, .5, .75), scale = "width") +
  xlab(NULL) + ylab("CA frequency") +
  scale_fill_manual(values = viridis::mako(n=5)[c(1,5)]) +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig~2b_", Sys.Date(), ".pdf"), pl, width = 6, height = 5)

 # C
CA_loci_df %>%
  filter((chromHMM_annotation == "intergenic" & ChIP_annotation == "unbound") | (chromHMM_annotation == "Ctcf" & TF == "Ctcf" & ChIP_annotation == "bound")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("intergenic", "Ctcf"))) %>% 
  group_by(chromHMM_annotation) %>% 
  summarise(median = median(width), .groups = "drop") %>%
  readr::write_delim(., col_names = TRUE, file = paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig~2c_values_", Sys.Date(), ".tsv"), delim = "\t", append = FALSE)
CA_loci_df %>% 
  filter((chromHMM_annotation == "intergenic" & ChIP_annotation == "unbound") | (chromHMM_annotation == "Ctcf" & TF == "Ctcf" & ChIP_annotation == "bound")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = rev(c("intergenic", "Ctcf")))) %>% 
  ggplot(aes(chromHMM_annotation, width)) +
  geom_violin(aes(fill = chromHMM_annotation), color = "black", linewidth = 1, alpha = 0.75, draw_quantiles = c(.25, .5, .75), scale = "width") +
  xlab(NULL) + ylab("CA width") +
  scale_fill_manual(values = viridis::mako(n=5, direction = -1)[c(1,5)]) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig~2c_", Sys.Date(), ".pdf"), pl, width = 8, height = 3)

# D
bound.sox2.motifs = TFBSs[unique(filter(CA_loci_df_dTAG, TF %in% c("Sox2", "Oct4") & ChIP_annotation == "bound")$TF.name)]
CA_loci_df_dTAG %>%
  mutate(plotting_category = case_when(
    str_detect(Sample, "^Sox2") & TF == "Sox2" & chromHMM_annotation == "intergenic" & ChIP_annotation == "unbound" ~ "unbound Sox2",
    str_detect(Sample, "^Sox2") & TF != "Sox2" & chromHMM_annotation == "intergenic" & ChIP_annotation == "unbound" ~ "unbound other TFs",
    str_detect(Sample, "^Sox2") & TF == "Sox2" & chromHMM_annotation == "enhancer" & ChIP_annotation == "bound" ~ "bound Sox2",
    str_detect(Sample, "^Sox2") & !(TF %in% c("Sox2", "Ctcf")) & chromHMM_annotation == "enhancer" & ChIP_annotation == "bound" ~ "bound other TFs"
  )) %>%
  filter(!is.na(plotting_category)) %>%
  mutate(treatment = gsub("Oct4|Sox2", "", Sample), Sample = gsub("2h|NT", "", Sample)) %>%
  dplyr::select(Sample, TF.name, TF, plotting_category, chromHMM_annotation, CA, treatment) %>%
  spread(treatment, CA) %>%
  na.omit() %>%
  mutate(plotting_category = factor(plotting_category, levels = c("bound other TFs", "bound Sox2", "unbound Sox2", "unbound other TFs", "Ctcf"))) -> plotting.df
polluted.TFBSs = filter_by_overlaps(TFBSs[unique(filter(plotting.df, plotting_category != "bound Sox2")$TF.name)], IRanges::resize(bound.sox2.motifs, 300, "center"))$absolute.idx
plotting.df %>%
  filter(!TF.name %in% polluted.TFBSs) %>%
  ggplot(aes(plotting_category, `2h`-NT)) +
  geom_violin(aes(fill = chromHMM_annotation), color = "black", linewidth = 1, alpha = 0.75, draw_quantiles = c(.25, .5, .75), scale = "width") +
  geom_vline(xintercept = 2.5, linewidth = 0.25) +
  ggpubr::stat_compare_means(
    method = "wilcox", label = "p.signif", comparisons = list(1:2,3:4),
    tip.length = .01, label.y = c(55, 55), step.increase = .05, size = 5) + 
  xlab(NULL) + ylab("CA delta (dTAG)") +
  scale_fill_manual(values = viridis::mako(n=5, direction = -1)[c(4,5)]) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.title = element_blank(), legend.position = "top", legend.background = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2d_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 6.5)

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
    ) -> quantiles

quantiles %>%
  ggplot(aes(width_0.50, frequency_0.50, color = Sample)) +
  geom_point(data = pl.df, aes(width, CA, color = Sample), alpha = .5, inherit.aes = FALSE) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = Sample), alpha = .5, linetype = 1, width = 15) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = Sample), alpha = .5, linetype = 1, height = 5) +
  scale_color_manual(values = c("black", "salmon")) +
  scale_y_continuous(breaks = as.integer(c(0,50,100)), limits = c(0,100)) +
  scale_x_continuous(breaks = as.integer(c(100,200,300)), limits = c(50,375)) +
  xlab("CA width (bp)") + ylab("CA frequency (%)") +
  theme_bw() +
  theme(text = element_text(size = 18), legend.title = element_blank(), legend.position = c(.8,.15), legend.background = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2e_", Sys.Date(), ".pdf"), pl, width = 4, height = 4)

# F
remove(sampleSheet, MySamples, GenomicTiles)
patch.single.site.plots(
  interpretable.master.table = data.frame(Sample = "amplicon_DE_data", TFBS.cluster = "AMP_16"), rank = 1, pool.replicates = TRUE, resize.size = 500, k = 24,
  partition.collapsing.dict = split(seq(24), seq(24))[c(18,10,6,24,17,20,4,8,19,2,11,9,12,13,3,23,7,5,1,14,15,16,21,22)],
  data.type = "WT_amplicon", remove.TFBS.labels = TRUE, deduplicate = TRUE, reutrn.chromatin.influence.df = TRUE
) -> pl
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2f_", Sys.Date(), ".png"), width = 25, height = 18, units = "cm", res = 300)
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
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig2g_", Sys.Date(), ".png"), width = 25, height = 18, units = "cm", res = 300)
pl$pl
dev.off()
pl.df %>% filter(TFBS.cluster == "AMP_16", Sample == "REST KO")

