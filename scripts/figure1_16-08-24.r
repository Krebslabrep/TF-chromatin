library(tidyverse)
source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
source("/g/krebs/barzaghi/Rscripts/IGV_plotting.R")

setwd("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/")

chromHMM = Load.chromHMM(GenomicTiles = TRUE)
GenomicTiles = Load.GenomicTiles(tiles.width = 80)
TFBSs = Load.TFBSs()

CA_loci_df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-29_chromatin.influence_Sonmezer_PooledReplicates.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
    (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
    (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
    mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

# B
CA_loci_df %>% 
  filter((chromHMM_annotation %in% c("intergenic", "enhancer", "promoter")) | (chromHMM_annotation == "Ctcf" & TF == "Ctcf")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("intergenic", "enhancer", "promoter", "Ctcf"))) %>%
  group_by(chromHMM_annotation) %>%
  summarise(
    width_0.25 = quantile(width, .25), frequency_0.25 = quantile(CA, 0.25), 
    width_0.50 = quantile(width, .50), frequency_0.50 = quantile(CA, 0.50), 
    width_0.75 = quantile(width, .75), frequency_0.75 = quantile(CA, 0.75), 
    .groups = "drop"
  ) -> quantiles
quantiles %>% dplyr::select(chromHMM_annotation, width_0.50, frequency_0.50)
quantiles %>%
  ggplot(aes(width_0.50, frequency_0.50, color = chromHMM_annotation)) +
  geom_text(aes(300, 0, label = "Errorbars represent inter-quartile ranges, IQR"), vjust = "inward", hjust = "inward", size = 4, show.legend = FALSE, inherit.aes = FALSE) +
  geom_text(aes(300, 100, label = "n=111,226"), vjust = "inward", hjust = "inward", size = 4, show.legend = FALSE, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = chromHMM_annotation), alpha = 1, linetype = 1, linewidth = 1, width = 15, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = chromHMM_annotation), alpha = 1, linetype = 1, linewidth = 1, height = 5, show.legend = FALSE) +
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_text(aes(label = chromHMM_annotation), color = "black", nudge_x = 10, nudge_y = -5, hjust = "left", size = 6, show.legend = FALSE) +
  scale_color_manual(values = viridis::mako(6)[c(1,3:5)], guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = as.integer(c(0,50,100)), limits = c(0,100)) +
  scale_x_continuous(breaks = as.integer(c(0,100,200,300)), limits = c(0,300)) +
  xlab("average width (bp)") + ylab("accessibility f (%) \n (active + linker)") +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1b_", Sys.Date(), ".pdf"), pl, width = 5, height = 4.5)
# 
# CA_loci_df %>% 
#   filter((chromHMM_annotation %in% c("intergenic", "enhancer", "promoter")) | (chromHMM_annotation == "Ctcf" & TF == "Ctcf")) %>%
#   mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("intergenic", "enhancer", "promoter", "Ctcf"))) %>%
#   filter(chromHMM_annotation != "intergenic") %>%
#   mutate(CA = CA_regulatory, width = width_regulatory) %>%
#   group_by(chromHMM_annotation) %>%
#   summarise(
#     width_0.25 = quantile(width, .25), frequency_0.25 = quantile(CA, 0.25), 
#     width_0.50 = quantile(width, .50), frequency_0.50 = quantile(CA, 0.50), 
#     width_0.75 = quantile(width, .75), frequency_0.75 = quantile(CA, 0.75), 
#     .groups = "drop"
#   ) -> quantiles
# quantiles %>%
#   ggplot(aes(width_0.50, frequency_0.50, color = chromHMM_annotation)) +
#   geom_text(aes(300, 0, label = "Errorbars represent inter-quartile ranges, IQR"), vjust = "inward", hjust = "inward", size = 4, show.legend = FALSE, inherit.aes = FALSE) +
#   geom_text(aes(300, 100, label = "n=99,066"), vjust = "inward", hjust = "inward", size = 4, show.legend = FALSE, inherit.aes = FALSE) +
#   geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = chromHMM_annotation), alpha = 1, linetype = 5, linewidth = 1, width = 15, show.legend = FALSE) +
#   geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = chromHMM_annotation), alpha = 1, linetype = 5, linewidth = 1, height = 5, show.legend = FALSE) +
#   geom_point(size = 3.5, show.legend = FALSE) +
#   geom_text(aes(label = chromHMM_annotation), color = "black", nudge_x = 10, nudge_y = -5, hjust = "left", size = 6, show.legend = FALSE) +
#   scale_color_manual(values = viridis::mako(6)[c(3:5)], guide = guide_legend(reverse = TRUE)) +
#   scale_y_continuous(breaks = as.integer(c(0,50,100)), limits = c(0,100)) +
#   scale_x_continuous(breaks = as.integer(c(0,100,200,300)), limits = c(0,300)) +
#   xlab("average width (bp) \n (active molecules)") + ylab("accessibility f (%) \n (active molecules)") +
#   theme_bw() +
#   theme(text = element_text(size = 18)) -> pl
# ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1b_transition_", Sys.Date(), ".pdf"), pl, width = 5, height = 4.5)

CA_loci_df %>% 
  filter((chromHMM_annotation %in% c("intergenic", "enhancer", "promoter")) | (chromHMM_annotation == "Ctcf" & TF == "Ctcf")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("intergenic", "enhancer", "promoter", "Ctcf"))) %>%
  mutate(nucleosome = 100 - CA, accessible = CA) %>%
  dplyr::select(chromHMM_annotation, nucleosome, accessible) %>%
  gather(fraction, frequency, nucleosome, accessible) %>%
  ggplot(aes(chromHMM_annotation, frequency, fill = fraction)) +
  geom_boxplot() +
  scale_fill_manual(breaks = c("accessible", "nucleosome"), values = c("#8FBC8F", "#009ACD")) +
  scale_y_continuous(breaks = as.integer(c(0,50,100)), limits = c(0,100)) +
  xlab(NULL) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.title = element_blank(), legend.position = "top") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1e_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5)

CA_loci_df %>% 
  filter((chromHMM_annotation %in% c("intergenic", "enhancer", "promoter")) | (chromHMM_annotation == "Ctcf" & TF == "Ctcf")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("intergenic", "enhancer", "promoter", "Ctcf"))) %>%
  mutate(linker = CA - CA_regulatory, active = CA_regulatory) %>%
  dplyr::select(chromHMM_annotation, linker, active) %>%
  gather(fraction, frequency, linker, active) %>% # group_by(fraction, chromHMM_annotation) %>% summarise(median(frequency))
  ggplot(aes(chromHMM_annotation, frequency, fill = fraction)) +
  geom_boxplot() +
  scale_fill_manual(breaks = c("active", "linker"), values = c("#8DA28D", "#E60808")) +
  scale_y_continuous(breaks = as.integer(c(0,50,100)), limits = c(0,100)) +
  xlab(NULL) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.title = element_blank(), legend.position = "top") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1f_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5)

# C
# CA_loci_df %>%
#   filter(chromHMM_annotation == "enhancer", ChIP_annotation == "bound") %>%
#   filter(tot.read.count >= 150, TF %in% c("Oct4", "Sox2")) %>%
#   filter(TFBS.cluster == "GenomicTile_33300452") %>%
#   mutate(Sample = "SMF_MM_TKO_DE_") -> enhancer
# TFs (left to right): Sox2, Oct4
enhancer = data.frame(Sample = "SMF_MM_TKO_DE_", TFBS.cluster = "GenomicTile_33300452")
partition.collapsing.dictinary = split(1:10,1:10)[c(1,7,2,8,9,4,6,10,3,5)]
patch.single.site.plots(
  interpretable.master.table = enhancer, rank = 1, k = 16,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE, deduplicate = FALSE,
  plotting.TFBSs = TFBSs[c("TFBS_14655497", "TFBS_5774468")]
) -> pl
process.CA.df(x = mutate(pl$chromatin.influence.df, ChIP = NA), cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) # 45%;15% | 267bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1c_", Sys.Date(), ".png"), width = 32, height = 25, units = "cm", res = 300)
pl$pl
dev.off()

RegionOfInterest = IRanges::resize(GenomicTiles["GenomicTile_33300452"], 9000, "center")

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
  sampleSheet = ChIP_data_dictionary$Oct4,
  samples = "Oct4",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 10,
  max.y.lim = 400,
  color = "black",
  y.labs = "Oct4"
) -> Oct4_track

plot_genomic_track(
  sampleSheet = ChIP_data_dictionary$Sox2,
  samples = "Sox2",
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

pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1c_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()

# D
# CA_loci_df %>%
#   filter(chromHMM_annotation == "promoter", ChIP_annotation == "bound") %>%
#   mutate(Sample = "SMF_MM_TKO_DE_") %>%
#   mutate(m = 70, r = abs(CA-m)) %>%
#   arrange(r) %>%
#   filter(TFBS.cluster == "GenomicTile_11037115") %>%
#   filter(tot.read.count >= 150, TF %in% c("Oct4", "Sox2")) -> promoter
# TFs (left to right): Sox2, Oct4
promoter = data.frame(Sample = "SMF_MM_TKO_DE_", TFBS.cluster = "GenomicTile_11037115")
partition.collapsing.dictinary = split(1:8,1:8)[c(5,3,2,8,6,4,1,7)]
patch.single.site.plots(
  interpretable.master.table = promoter, rank = 1, k = 8,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE
) -> pl
process.CA.df(x = mutate(pl$chromatin.influence.df, ChIP = NA), cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) # 76% | 230bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1d_", Sys.Date(), ".png"), width = 32, height = 25, units = "cm", res = 300)
pl$pl
dev.off()

RegionOfInterest = IRanges::resize(GenomicTiles["GenomicTile_11037115"], 9000, "center")

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
  sampleSheet = ChIP_data_dictionary$Oct4,
  samples = "Oct4",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 10,
  max.y.lim = 400,
  color = "black", 
  y.labs = "Oct4"
) -> Oct4_track

plot_genomic_track(
  sampleSheet = ChIP_data_dictionary$Sox2,
  samples = "Sox2",
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

pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1d_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()

# E
# CA_loci_df %>%
#   filter(chromHMM_annotation == "intergenic" & ChIP_annotation == "unbound") %>%
#   mutate(m = median(CA), r = abs(CA-m)) %>%
#   filter(tot.read.count >= 100 & CA <= m+5 & CA >= m-5 & width < 110) %>%
#   dplyr::select(-TF) %>%
#   left_join(., dplyr::select(x, TF, TF.name), by = "TF.name") %>%
#   filter(TF %in% c("Oct4", "Sox2")) %>%
#   mutate(Sample = "SMF_MM_TKO_DE_") -> intergenic.single.sites
intergenic.single.sites = data.frame(Sample = "SMF_MM_TKO_DE_", TFBS.cluster = "GenomicTile_29023109")
partition.collapsing.dictinary = split(1:9,1:9)[c(6,8,2,1,7,3,4,5,9)]
patch.single.site.plots(
  interpretable.master.table = intergenic.single.sites, rank = 4, k = 12, 
  plotting.TFBSs = TFBSs["TFBS_5388418"],
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE
) -> pl
process.CA.df(x = mutate(pl$chromatin.influence.df, ChIP = NA), cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) # 4%;29% | 81bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1e_", Sys.Date(), ".png"), width = 27, height = 20, units = "cm", res = 300)
pl$pl
dev.off()

RegionOfInterest = IRanges::resize(GenomicTiles["GenomicTile_29023109"], 9000, "center")

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
  sampleSheet = ChIP_data_dictionary$Oct4,
  samples = "Oct4",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 10,
  max.y.lim = 400,
  color = "black", 
  y.labs = "Oct4"
) -> Oct4_track

plot_genomic_track(
  sampleSheet = ChIP_data_dictionary$Sox2,
  samples = "Sox2",
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

pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1e_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()