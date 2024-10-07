library(tidyverse)
library(magrittr)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SingleMoleculeFootprinting)
source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
source("/g/krebs/barzaghi/Rscripts/IGV_plotting.R")
detach("package:plyranges")

setwd("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/")

chromHMM = Load.chromHMM(GenomicTiles = TRUE)
GenomicTiles = Load.GenomicTiles(tiles.width = 80, GeneAssociation = FALSE, CGI = FALSE)
TFBSs = Load.TFBSs()

CA_loci_df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-08-22_chromatin.influence_Sonmezer_PooledReplicates.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

CA_loci_df_f1 = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-08-22_chromatin.influence_F1_PooledReplicates.df.qs") %>%
  process.CA.df_f1(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>% filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
    (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
    (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

CA_loci_df_dTAG = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-30_chromatin.influence_Oct4_Sox2_KD.df.qs") %>%
  filter(str_detect(Sample, "Sox2")) %>%
  process.CA.df(., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
    (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
    (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

TFBSs[unique(c(filter(CA_loci_df, ChIP_annotation == "bound")$TF.name, filter(CA_loci_df_f1, ChIP_annotation == "bound")$TF.name))] %>%
  SingleMoleculeFootprinting::Arrange_TFBSs_clusters(., add.single.TFs = FALSE, max_cluster_width = 300, max_cluster_size = 10) -> TFBS.clusters
TFBS.cluster.compositions = unlist(TFBS.clusters$ClusterComposition)
TFBS.cluster.compositions$cluster.id = gsub("\\..*", "", names(TFBS.cluster.compositions))
TFBS.cluster.compositions %>%
  data.frame() %>%
  group_by(cluster.id) %>%
  arrange(start) %>%
  mutate(
    nr.motifs = n(),
    border.motif = ifelse(absolute.idx == rev(absolute.idx)[1] | absolute.idx == absolute.idx[1], TRUE, FALSE)) %>%
  ungroup() %>%
  dplyr::select(absolute.idx, cluster.id, border.motif, nr.motifs) %>%
  dplyr::rename("TF.name" = "absolute.idx") -> TFBS.cluster.composition.df
data.frame(
  TF.name = unique(c(
    filter(filter(CA_loci_df, ChIP_annotation == "bound"), !TF.name %in% TFBS.cluster.composition.df$TF.name)$TF.name,
    filter(filter(CA_loci_df_f1, ChIP_annotation == "bound"), !TF.name %in% TFBS.cluster.composition.df$TF.name)$TF.name
  ))) %>%
  mutate(cluster.id = TF.name, border.motif = 1, nr.motifs = 1) %>%
  rbind(., TFBS.cluster.composition.df) -> TFBS.cluster.composition.df

# B left-to-right: FOXD3, STAT3, KLF4, KLF4
partition.collapsing.dictionary = split(1:12,1:12)[c(12,7,6,1,5,2,4,11,3,9,10,8)]
patch.single.site.plots(
  interpretable.master.table = data.frame(Sample = "CTKO", TFBS.cluster = "GenomicTile_17185595"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictionary, 
  data.type = "F1_bait.capture", remove.TFBS.labels = TRUE
) -> pl
process.CA.df_f1(x = mutate(pl$chromatin.influence.df[1,], ChIP = NA), cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) # 50%-30% | 301bp-272bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3b_", Sys.Date(), ".png"), width = 30, height = 30, units = "cm", res = 300)
pl$pl
dev.off()

RegionOfInterest = IRanges::resize(GenomicTiles["GenomicTile_17185595"], 9000, "center")
plot_genomic_track(
  # sampleSheet = "/g/krebs/barzaghi/analyses/nf-core_runs/130923_BL6_Spret_F1_ATAC/aln_merged_deduplicated/samplesheet_PooledReplicates.txt",
  # samples = "BL6_Spret_F1",
  sampleSheet = "/g/krebs/barzaghi/HTS/ATAC-seq/Heard_NatGen_2017/aln/Qinput_male.txt",
  samples = "male_R1",
  RegionOfInterest = RegionOfInterest,
  allelic = TRUE,
  tile.width = 200,
  tile.step = 25,
  max.y.lim = 30,
  color = c("black", "sienna"), 
  y.labs = c("Bl6 ATAC", "Cast ATAC"),
  delta = FALSE,
  normalise = FALSE
) -> ATAC_track
plot_genomic_track(
  sampleSheet = ChIP_data_dictionary$Klf4,
  samples = "Klf4",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 25,
  max.y.lim = 300,
  color = "black", 
  y.labs = "Klf4"
) -> Klf4_track
plot_genomic_track(
  sampleSheet = "/g/krebs/barzaghi/HTS/MNase/Qinput.txt",
  samples = "MNase_mESC_129S6/SvEvTac_WT",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200, 
  tile.step = 25,
  max.y.lim = 160,
  color = "black", 
  y.labs = "MNAse"
) -> MNase_track
p_final <- ATAC_track + Klf4_track + MNase_track +
  plot_layout(ncol = 1, heights = c(1/2, 1/4, 1/4))
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3e_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()

# C
inactive.median = as.integer(quantile(filter(CA_loci_df, chromHMM_annotation == "intergenic", ChIP_annotation == "unbound")$CA_regulatory, 0.50))
CA_loci_df_f1 %>%
  filter(motif.change == "l.o.f.", ChIP_annotation == "bound", TF != "Ctcf") %>%
  mutate(ref = CA_regulatory_R, alt = CA_regulatory_A, delta = CA.delta) %>%
  arrange(ref) %>%
  mutate(rank = seq(nrow(.))) -> pl.df
pl.df %>%
  dplyr::select(rank, ref, alt) %>%
  gather(allele, CA, ref, alt) %>%
  ggplot(aes(rank, CA, color = allele)) +
  geom_point(size = .25, alpha = .1) +
  geom_smooth(method = "gam") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = inactive.median, fill = viridis::mako(n=1)[1], alpha = .25) +
  xlab(NULL) + ylab(expression(paste("accessibility ", italic("f"), " (%)"))) +
  scale_color_manual(values = c("black", "salmon"), breaks = c("ref", "alt")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(0,inactive.median,50,100), limits = c(0,100)) + 
  theme_bw() +
  theme(text = element_text(size = 18), legend.position=c(.2, .8), legend.title = element_blank(), legend.background = element_blank()) -> up
pl.df %>%
  dplyr::select(rank, delta) %>%
  ggplot(aes(rank, delta, color = "delta")) +
  geom_point(size = .25, alpha = .1) +
  geom_smooth(method = "gam") +
  geom_hline(yintercept = -10) +
  xlab("motifs rank by ref activity") + ylab(expression(paste("delta (alt - ref)"))) +
  scale_color_manual(values = c("blue"), breaks = c("delta")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(-100,-50,-10,0,50), limits = c(-100,50)) +
  theme_bw() +
  theme(text = element_text(size = 18), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = c(.25, .2), legend.title = element_blank(), legend.background = element_blank()) -> down
up / down + patchwork::plot_layout(heights = c(10/25, 15/25), axes = "collect") -> pl
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3c_", Sys.Date(), ".png"), units = "cm", res = 300, width = 12, height = 18)
pl
dev.off()

# CA_loci_df_f1 %>%
#   filter(ChIP_annotation == "bound") %>%
#   mutate(class = case_when(
#     motif.change == "l.o.f." & TF == "Ctcf" ~ "diverged Ctcf",
#     motif.change == "l.o.f." & TF != "Ctcf" ~ "diverged other",
#     motif.change == "no" ~ "conserved",
#   )) %>%
#   mutate(ref = CA_regulatory_R, alt = CA_regulatory_A, delta = CA.delta) %>%
#   arrange(ref) %>%
#   group_by(class) %>% mutate(rank = seq_along(class)/length(class)) %>% ungroup() %>%
#   mutate(rank = cut(rank, breaks = seq(0,1,.1), include.lowest = TRUE, labels = seq(10,100,10))) %>%
#   dplyr::select(rank, delta, class) %>%
#   mutate(class = factor(class, levels = c("conserved", "diverged other", "diverged Ctcf"))) -> pl.df
# 
# pl.df %>%
#   ggplot(aes(rank, delta, fill = class)) +
#   geom_boxplot() +
#   xlab("TF motifs rank") + ylab("CA frequencty delta (%)") +
#   scale_fill_manual(values = viridis::cividis(3)[c(1,3,2)], breaks = c("conserved", "diverged other", "diverged Ctcf")) +
#   # scale_x_continuous(breaks = c(0,1), limits = c(0,1)) +
#   scale_y_continuous(breaks = c(-80,-10,0,40), limits = c(-80,45)) +
#   theme_bw() +
#   theme(text = element_text(size = 18), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = c(.25, .2), legend.title = element_blank(), legend.background = element_blank())
# 
# pl.df %>% ggplot(aes(rank, delta, color = class)) +
#   geom_point(data = filter(pl.df, class != "conserved"), aes(rank, delta, color = class), size = .5, alpha = .5, inherit.aes = FALSE) +
#   geom_smooth(method = "loess", span = .5, se = FALSE) +
#   geom_hline(yintercept = -10, linetype = 2, color = "grey") +
#   xlab("TF motifs rank") + ylab("CA frequencty delta (%)") +
#   scale_color_manual(values = c("black", viridis::mako(6)[c(2,5)]), breaks = c("conserved", "diverged Ctcf", "diverged other")) +
#   scale_x_continuous(breaks = c(0,1), limits = c(0,1)) +
#   scale_y_continuous(breaks = c(-80,-10,0,40), limits = c(-80,45)) +
#   theme_bw() +
#   theme(text = element_text(size = 18), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = c(.25, .2), legend.title = element_blank(), legend.background = element_blank())

# D
CA_loci_df_f1 %>%
  filter(ChIP_annotation == "unbound", motif.change == "l.o.f.") %>%
  dplyr::select(Sample, TF, TF.name, CA_regulatory_R, CA_regulatory_A, width_regulatory_R, width_regulatory_A) %>%
  mutate(nr.motifs = 0) -> x.0

CA_loci_df_f1 %>% 
  filter(ChIP_annotation == "bound") %>%
  full_join(., TFBS.cluster.composition.df, by = "TF.name") %>% #filter(!is.na(cluster.id)) %>% 
  mutate(cluster.id = ifelse(is.na(cluster.id), TF.name, cluster.id), nr.motifs = ifelse(is.na(nr.motifs), 1, nr.motifs), border.motif = ifelse(is.na(border.motif), TRUE, border.motif)) %>%
  group_by(Sample, cluster.id) %>% 
  mutate(nr.lofs = sum(motif.change == "l.o.f.")) %>% 
  filter(any(motif.change == "l.o.f." & border.motif)) %>%
  ungroup() %>%
  filter(motif.change == "l.o.f." & nr.lofs == 1) %>% dplyr::select(-nr.lofs, -border.motif) %>%
  group_by(nr.motifs) %>% filter(n() > 100) %>% ungroup() %>%
  dplyr::select(Sample, TF, TF.name, nr.motifs, CA_regulatory_R, CA_regulatory_A, width_regulatory_R, width_regulatory_A) %>%
  rbind(., x.0) %>%
  pivot_longer(cols = c(CA_regulatory_R, CA_regulatory_A, width_regulatory_R, width_regulatory_A)) %>%
  mutate(name = gsub("_regulatory", "", gsub("CA", "frequency", name))) %>%
  separate("name", into = c("name", "allele"), sep = "_") %>%
  spread(name, value) %>%
  mutate(allele = factor(ifelse(allele == "R", "ref", "alt"), levels = c("ref", "alt"))) %>%
  group_by(allele, nr.motifs) %>%
  summarise(
    frequency_0.25 = quantile(frequency, 0.25, na.rm=TRUE), frequency_0.50 = quantile(frequency, 0.50, na.rm=TRUE), frequency_0.75 = quantile(frequency, 0.75, na.rm=TRUE), 
    width_0.25 = quantile(width, 0.25, na.rm=TRUE), width_0.50 = quantile(width, 0.50, na.rm=TRUE), width_0.75 = quantile(width, 0.75, na.rm=TRUE), 
    .groups = "drop"
  ) %>%
  ggplot(aes(width_0.50, frequency_0.50, color = allele)) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75), alpha = .5, linetype = 1, width = 15, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75), alpha = .5, linetype = 1, height = 5, show.legend = FALSE) +
  geom_point(size = 2) +
  geom_text(aes(label = nr.motifs), nudge_x = -5, nudge_y = 3, show.legend = FALSE) +
  scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
  scale_x_continuous(breaks = c(0,100,200,300), limits = c(0,300)) +
  ggbreak::scale_x_break(breaks = c(5,100)) +
  xlab("average width (bp)") + ylab(expression(paste("accessibility ", italic("f"), " (%)"))) +
  scale_color_manual(values = c("black","salmon"), breaks = c("ref", "alt")) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.title = element_blank(), legend.background = element_blank(), axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3d_", Sys.Date(), ".pdf"), pl, width = 5.5, height = 4.5, onefile = FALSE)

# CA_loci_df_f1 %>% 
#   filter(ChIP_annotation == "bound") %>%
#   full_join(., TFBS.cluster.composition.df, by = "TF.name") %>% #filter(!is.na(cluster.id)) %>% 
#   mutate(cluster.id = ifelse(is.na(cluster.id), TF.name, cluster.id), nr.motifs = ifelse(is.na(nr.motifs), 1, nr.motifs), border.motif = ifelse(is.na(border.motif), TRUE, border.motif)) %>%
#   group_by(Sample, cluster.id) %>% 
#   mutate(nr.lofs = sum(motif.change == "l.o.f.")) %>% 
#   filter(any(motif.change == "l.o.f." & border.motif)) %>%
#   ungroup() %>%
#   filter(motif.change == "l.o.f." & nr.lofs == 1) %>% dplyr::select(-nr.lofs, -border.motif) %>%
#   group_by(nr.motifs) %>% filter(n() > 100) %>% ungroup() %>%
#   dplyr::select(Sample, TF, TF.name, nr.motifs, CA_regulatory_R, CA_regulatory_A, width_regulatory_R, width_regulatory_A) %>%
#   rbind(., x.0) %>%
#   pivot_longer(cols = c(CA_regulatory_R, CA_regulatory_A, width_regulatory_R, width_regulatory_A)) %>%
#   mutate(name = gsub("_regulatory", "", gsub("CA", "frequency", name))) %>%
#   separate("name", into = c("name", "allele"), sep = "_") %>%
#   spread(name, value) %>%
#   mutate(allele = factor(ifelse(allele == "R", "ref", "alt"), levels = c("ref", "alt"))) %>%
#   mutate(nr.motifs = factor(nr.motifs, levels = seq(0,3))) %>%
#   ggplot(aes(nr.motifs, frequency, fill = allele)) +
#   geom_boxplot(color = "black") +
#   ggpubr::stat_compare_means(method = "wilcox", label = "p.signif", size = 5) +
#   scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
#   xlab("Nb motifs per CRE") + ylab("frequency (%)") +
#   scale_fill_manual(values = c("grey30","salmon"), breaks = c("ref", "alt")) +
#   theme_bw() +
#   theme(text = element_text(size = 18)) -> pl
# ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3d_", Sys.Date(), ".pdf"), pl, width = 5.5, height = 5)

# E
inactive.median = as.integer(quantile(filter(CA_loci_df, chromHMM_annotation == "intergenic", ChIP_annotation == "unbound")$CA_regulatory, 0.50))
CA_loci_df_dTAG %>%
  filter(TF == "Sox2", ChIP_annotation == "bound") %>%
  dplyr::select(-c(chromHMM_annotation, ChIP_annotation, tot.read.count, width, CA, CA_regulatory_count, width_regulatory)) %>%
  spread(Sample, CA_regulatory) %>%
  dplyr::rename("NT" = "Sox2NT", "2h" = "Sox22h") %>%
  na.omit() %>%
  mutate(delta = `2h` - NT) %>%
  arrange(NT) %>%
  mutate(rank = seq(nrow(.))) -> pl.df
pl.df %>%
  dplyr::select(rank, NT, `2h`) %>%
  gather(treatment, CA, NT, `2h`) %>%
  ggplot(aes(rank, CA, color = treatment)) +
  geom_point(size = .25, alpha = .1) +
  geom_smooth(method = "gam") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = inactive.median, fill = viridis::mako(n=1)[1], alpha = .25) +
  xlab(NULL) + ylab(expression(paste("accessibility ", italic("f"), " (%)"))) +
  scale_color_manual(values = c("black", "salmon"), breaks = c("NT", "2h")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(0,inactive.median,50,100), limits = c(0,100)) + 
  theme_bw() +
  theme(text = element_text(size = 18), legend.position=c(.2, .8), legend.title = element_blank(), legend.background = element_blank()) -> up
pl.df %>%
  dplyr::select(rank, delta) %>%
  ggplot(aes(rank, delta, color = "delta")) +
  geom_point(size = .25, alpha = .1) +
  geom_smooth(method = "gam") +
  geom_hline(yintercept = -10) +
  xlab("motifs rank by NT activity") + ylab(expression(paste("delta (2h - NT)"))) +
  scale_color_manual(values = c("blue"), breaks = c("delta")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(-100,-50,-10,0,50), limits = c(-100,50)) +
  theme_bw() +
  theme(text = element_text(size = 18), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = c(.25, .2), legend.title = element_blank(), legend.background = element_blank()) -> down
up / down + patchwork::plot_layout(heights = c(10/25, 15/25), axes = "collect") -> pl
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3e_", Sys.Date(), ".png"), units = "cm", res = 300, width = 12, height = 18)
pl
dev.off()

# F
CA_loci_df_dTAG %>%
  filter(ChIP_annotation == "unbound") %>%
  dplyr::select(Sample, TF, TF.name, CA_regulatory, width_regulatory) %>%
  pivot_wider(id_cols = c(TF.name, TF), names_from = Sample, values_from = c(CA_regulatory, width_regulatory)) %>% na.omit() %>%
  pivot_longer(cols = contains("_Sox2")) %>%
  separate(name, c("name", "Sample"), "_Sox2") %>%
  mutate(Sample = paste0("Sox2", Sample)) %>%
  spread(name, value) %>%
  mutate(nr.motifs = 0, nr.sox2 = 0) -> x.0

CA_loci_df_dTAG %>%
  filter(ChIP_annotation == "bound") %>%
  full_join(., TFBS.cluster.composition.df, by = "TF.name") %>%
  filter(!is.na(cluster.id) & !is.na(Sample)) %>%
  group_by(Sample, cluster.id) %>% mutate(nr.sox2 = sum(TF == "Sox2")) %>% ungroup() %>%
  filter(nr.sox2 == 0 | TF == "Sox2") %>%
  pivot_wider(
    id_cols = c(TFBS.cluster, TF.name, TF, chromHMM_annotation, ChIP_annotation, cluster.id, border.motif, nr.motifs, nr.sox2), 
    names_from = Sample, values_from = c(tot.read.count, width, CA, CA_regulatory_count, CA_regulatory, width_regulatory)) %>% na.omit() %>%
  pivot_longer(cols = contains("_Sox2")) %>% separate(name, c("name", "Sample"), "_Sox2") %>% mutate(Sample = paste0("Sox2", Sample)) %>% spread(name, value) %>%
  group_by(nr.sox2, nr.motifs, Sample) %>% filter(n() > 100) %>% ungroup() %>%
  dplyr::select(Sample, TF, TF.name, nr.motifs, nr.sox2, CA_regulatory, width_regulatory) %>%
  dplyr::rename(frequency = CA_regulatory, width = width_regulatory) %>%
  mutate(Sample = factor(ifelse(Sample == "Sox22h", "2h", "NT"), levels = c("NT", "2h"))) %>%
  mutate(nr.motifs = as.factor(nr.motifs), nr.sox2 = as.factor(nr.sox2)) %>%
  pivot_longer(cols = c(frequency, width)) %>%
  spread(Sample, value) %>%
  mutate(delta = `2h` - NT) %>%
  dplyr::select(-NT, -`2h`) %>%
  spread(name, delta) %>%
  group_by(nr.sox2) %>%
  summarise(
    frequency_0.25 = quantile(frequency, 0.25, na.rm=TRUE), frequency_0.50 = quantile(frequency, 0.50, na.rm=TRUE), frequency_0.75 = quantile(frequency, 0.75, na.rm=TRUE), 
    width_0.25 = quantile(width, 0.25, na.rm=TRUE), width_0.50 = quantile(width, 0.50, na.rm=TRUE), width_0.75 = quantile(width, 0.75, na.rm=TRUE), 
    .groups = "drop"
  ) -> pl.df

pl.df %>%
  ggplot(aes(width_0.50, frequency_0.50)) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75, color = nr.sox2), linetype = 1, linewidth = 1, width = 3, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75, color = nr.sox2), linetype = 1, linewidth = 1, height = 1, show.legend = FALSE) +
  geom_point(aes(color = nr.sox2), size = 3.5, show.legend = FALSE) +
  geom_text(aes(label = nr.sox2), nudge_x = -1.5, nudge_y = .75, show.legend = FALSE) +
  xlab("average width delta (2h - NT)") + ylab("accessibility delta (2h - NT)") +
  scale_color_manual(values = colorspace::sequential_hcl(n = 7, palette = "Blues3", rev = TRUE)[c(3,5,7)], breaks = 0:2) +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3f_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5)



