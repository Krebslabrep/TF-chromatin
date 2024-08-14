library(tidyverse)
source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
source("/g/krebs/barzaghi/Rscripts/IGV_plotting.R")

setwd("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/")

chromatin.influence.df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-29_chromatin.influence_Sonmezer_PooledReplicates.df.qs")
chromHMM = Load.chromHMM(GenomicTiles = TRUE)
GenomicTiles = Load.GenomicTiles(tiles.width = 80)
TFBSs = Load.TFBSs()

CA_loci_df = process.CA.df(chromatin.influence.df, cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient)

# B
CA_loci_df %>%
  filter(chromHMM_annotation == "enhancer", ChIP_annotation == "bound") %>%
  filter(tot.read.count >= 150, TF %in% c("Oct4", "Sox2")) %>%
  filter(TFBS.cluster == "GenomicTile_33300452") %>%
  mutate(Sample = "SMF_MM_TKO_DE_") -> enhancer

# TFs (left to right): Sox2, Oct4
partition.collapsing.dictinary = split(1:10,1:10)[c(1,7,2,4,6,10,8,9,3,5)]
patch.single.site.plots(
  interpretable.master.table = enhancer, rank = 1, k = 16,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE
) -> pl
mutate(pl$chromatin.influence.df, CA = (accessible.read.count/tot.read.count)*100)$CA[1] # 60%
max(pl$chromatin.influence.df$acc.width.distro[[1]]) # 363bp
sum((pl$chromatin.influence.df$acc.read.count.distro[[1]]/pl$chromatin.influence.df$tot.read.count[[1]]*100)[which(pl$chromatin.influence.df$acc.width.distro[[1]]>=100)]) # 45%
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1b_", Sys.Date(), ".png"), width = 32, height = 25, units = "cm", res = 300)
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

pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1b_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()

# C
CA_loci_df %>% 
  filter(ChIP_annotation == "bound") %>%
  filter((chromHMM_annotation %in% c("enhancer", "bivalent_promoter", "promoter")) | (chromHMM_annotation == "Ctcf" & TF == "Ctcf")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("enhancer", "bivalent_promoter", "promoter", "Ctcf"))) %>%
  group_by(chromHMM_annotation) %>% 
  summarise(median = median(CA), .groups = "drop")
CA_loci_df %>% 
  filter(ChIP_annotation == "bound") %>%
  filter((chromHMM_annotation %in% c("enhancer", "bivalent_promoter", "promoter")) | (chromHMM_annotation == "Ctcf" & TF == "Ctcf")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("enhancer", "bivalent_promoter", "promoter", "Ctcf"))) %>%
  ggplot(aes(chromHMM_annotation, CA)) +
  geom_violin(aes(fill = chromHMM_annotation), color = "black", linewidth = 1, alpha = 0.75, draw_quantiles = c(.25, .5, .75), scale = "width") +
  xlab(NULL) + ylab("CA frequency") +
  scale_fill_manual(values = viridis::mako(n=5)[2:5]) +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1c_", Sys.Date(), ".pdf"), pl, width = 7, height = 5)

# D
CA_loci_df %>% 
  filter(ChIP_annotation == "bound") %>%
  filter((chromHMM_annotation %in% c("enhancer", "bivalent_promoter", "promoter")) | (chromHMM_annotation == "Ctcf" & TF == "Ctcf")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = c("enhancer", "bivalent_promoter", "promoter", "Ctcf"))) %>% 
  group_by(chromHMM_annotation) %>% 
  summarise(median = median(width), .groups = "drop")
CA_loci_df %>% 
  filter(ChIP_annotation == "bound") %>%
  filter((chromHMM_annotation %in% c("enhancer", "bivalent_promoter", "promoter")) | (chromHMM_annotation == "Ctcf" & TF == "Ctcf")) %>%
  mutate(chromHMM_annotation = factor(chromHMM_annotation, levels = rev(c("enhancer", "bivalent_promoter", "promoter", "Ctcf")))) %>% 
  ggplot(aes(chromHMM_annotation, width)) +
  geom_violin(aes(fill = chromHMM_annotation), color = "black", linewidth = 1, alpha = 0.75, draw_quantiles = c(.25, .5, .75), scale = "width") +
  xlab(NULL) + ylab("CA width") +
  scale_fill_manual(values = viridis::mako(n=5, direction = -1)[1:4]) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1d_", Sys.Date(), ".pdf"), pl, width = 8, height = 5)

# E
CA_loci_df %>%
  filter(chromHMM_annotation == "promoter", ChIP_annotation == "bound") %>%
  mutate(Sample = "SMF_MM_TKO_DE_") %>%
  mutate(m = 70, r = abs(CA-m)) %>%
  arrange(r) %>%
  filter(TFBS.cluster == "GenomicTile_11037115") %>%
  filter(tot.read.count >= 150, TF %in% c("Oct4", "Sox2")) -> promoter

# TFs (left to right): Sox2, Oct4
partition.collapsing.dictinary = split(1:8,1:8)[c(5,3,2,8,6,4,1,7)]
patch.single.site.plots(
  interpretable.master.table = promoter, rank = 1, k = 8,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE
) -> pl
mutate(pl$chromatin.influence.df, CA = (accessible.read.count/tot.read.count)*100)$CA[1] # 76%
max(pl$chromatin.influence.df$acc.width.distro[[1]]) # 254bp
sum((pl$chromatin.influence.df$acc.read.count.distro[[1]]/pl$chromatin.influence.df$tot.read.count[[1]]*100)[which(pl$chromatin.influence.df$acc.width.distro[[1]]>=100)]) # 76%
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

# F
CA_loci_df %>%
  filter(chromHMM_annotation == "Ctcf" & TF == "Ctcf" & ChIP_annotation == "bound") %>%
  mutate(Sample = "SMF_MM_TKO_DE_") %>%
  filter(CA > 85, width > 160 & width < 200) %>%
  filter(TFBS.cluster == "GenomicTile_1014499") %>%
  filter(tot.read.count >= 150) -> insulators

partition.collapsing.dictinary = split(1:12,1:12)[c(8,3,6,2,1,5,11,10,12,7,4,9)]
patch.single.site.plots(
  interpretable.master.table = insulators, rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE, plotting.TFBSs = TFBSs[insulators$TF.name]
) -> pl

png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1f_", Sys.Date(), ".png"), width = 32, height = 25, units = "cm", res = 300)
pl$pl
dev.off()

RegionOfInterest = IRanges::resize(GenomicTiles["GenomicTile_1014499"], 9000, "center")

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
  sampleSheet = ChIP_data_dictionary$Ctcf,
  samples = "Ctcf",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 10,
  max.y.lim = 500,
  color = "black", 
  y.labs = "Ctcf"
) -> Ctcf_track

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

p_final <- DHS_track + Ctcf_track + MNase_track + 
  plot_layout(ncol = 1, heights = c(1/4, 1/4, 1/4, 1/4))

pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig1f_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()



