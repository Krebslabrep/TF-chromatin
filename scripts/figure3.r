library(tidyverse)
library(magrittr)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SingleMoleculeFootprinting)
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/utils.r")
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/source_FootprintCharter_functions.r")
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/IGV_plotting.r")
detach("package:plyranges")

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

# CA_loci_df_dTAG = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2025-02-24_chromatin.influence_Sox2_dTAG_run2_run3_PooledReps.df.qs") %>%
CA_loci_df_dTAG = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2025-09-05_chromatin.influence_Sox2_dTAG_2h_24h_PooledReps.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
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

# B
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
  xlab("motifs rank by ref CA frequency") + ylab("CA frequency (%)") +
  scale_color_manual(values = c("black", "salmon"), breaks = c("ref", "alt")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(0,inactive.median,50,100), limits = c(0,100)) + 
  theme_bw() +
  theme(text = element_text(size = 18), legend.position=c(.25, .9), legend.title = element_blank(), legend.background = element_blank()) -> pl
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3b_", Sys.Date(), ".pdf"), width = 4.5, height = 4.5)
pl
dev.off()

# C
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
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3c_", Sys.Date(), ".pdf"), pl, width = 5.5, height = 4.5, onefile = FALSE)

# D left-to-right: ZIC3, KLF4, OCT4
partition.collapsing.dictionary = split(1:12,1:12)[c(2,9,4,6,7,8,3,10,5,11,12,1)]
patch.single.site.plots(
  interpretable.master.table = data.frame(Sample = "STKO", TFBS.cluster = "GenomicTile_19187323"), 
  rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictionary, 
  data.type = "F1_bait.capture", remove.TFBS.labels = TRUE
) -> pl
pl$chromatin.influence.df %>%
  filter(Sample == "STKO") %>%
  mutate(ChIP = NA) %>%
  process.CA.df_f1(cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  dplyr::select(CA_regulatory_R, CA_regulatory_A, width_regulatory_R, width_regulatory_A) # 63%-52% | 229bp-221bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3d_", Sys.Date(), ".png"), width = 30, height = 30, units = "cm", res = 300)
pl$pl
dev.off()

# F
inactive.median = as.integer(quantile(filter(CA_loci_df, chromHMM_annotation == "intergenic", ChIP_annotation == "unbound")$CA_regulatory, 0.50))
CA_loci_df_dTAG %>%
  filter(TF == "Sox2", ChIP_annotation == "bound") %>%
  dplyr::select(-c(chromHMM_annotation, ChIP_annotation, tot.read.count, width, CA, CA_regulatory_count, width_regulatory)) %>%
  spread(Sample, CA_regulatory) %>%
  dplyr::rename("untreated_2h" = "Sox2_2h_NT", "2h dTAG" = "Sox2_2h_T", "untreated_24h" = "Sox2_24h_NT", "24h dTAG" = "Sox2_24h_T") %>%
  na.omit() %>%
  # mutate(delta = `2h dTAG` - untreated) %>%
  arrange(untreated_2h) %>%
  mutate(rank = seq(nrow(.))) -> pl.df
pl.df %>%
  dplyr::select(rank, untreated_2h, untreated_24h, `2h dTAG`, `24h dTAG`) %>%
  gather(treatment, CA, untreated_2h, untreated_24h, `2h dTAG`, `24h dTAG`) %>%
  ggplot(aes(rank, CA, color = treatment)) +
  geom_point(size = .25, alpha = .1) +
  geom_smooth(method = "gam") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = inactive.median, fill = viridis::mako(n=1)[1], alpha = .25) +
  xlab("motifs rank by untreated CA frequency") + ylab("CA frequency (%)") +
  # scale_color_manual(values = c("black", "salmon"), breaks = c("untreated", "2h dTAG")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(0,inactive.median,50,100), limits = c(0,100)) + 
  theme_bw() +
  theme(text = element_text(size = 18), legend.position=c(.25, .9), legend.title = element_blank(), legend.background = element_blank()) -> pl
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3f_", Sys.Date(), ".pdf"), width = 4.5, height = 4.5)
pl
dev.off()

# G
CA_loci_df_dTAG %>%
  filter(ChIP_annotation == "bound") %>%
  full_join(., TFBS.cluster.composition.df, by = "TF.name") %>%
  filter(!is.na(cluster.id) & !is.na(Sample)) %>%
  group_by(Sample, cluster.id) %>% mutate(nr.sox2 = as.factor(sum(TF == "Sox2"))) %>% ungroup() %>%
  filter(nr.sox2 == 0 | TF == "Sox2") %>%
  dplyr::select(Sample, TF, TF.name, nr.sox2, CA_regulatory) %>%
  dplyr::rename(frequency = CA_regulatory) %>%
  mutate(Sample = factor(case_when(
    Sample == "Sox2_2h_NT" ~ "untreated_2h", 
    Sample == "Sox2_2h_T" ~ "2h dTAG", 
    Sample == "Sox2_24h_NT" ~ "untreated_24h", 
    Sample == "Sox2_24h_T" ~ "24h dTAG"
  ), levels = c("untreated_2h", "2h dTAG", "untreated_24h", "24h dTAG")
  )) %>%
  group_by(nr.sox2, Sample) %>% filter(n() > 100) %>% ungroup() %>%
  spread(Sample, frequency) %>%
  mutate(delta_2h = `2h dTAG` - untreated_2h, delta_24h = `24h dTAG` - untreated_24h) %>%
  dplyr::select(-untreated_2h, -`2h dTAG`, -untreated_24h, -`24h dTAG`) %>%
  gather(time, delta, delta_2h, delta_24h) %>%
  mutate(time = factor(ifelse(time == "delta_2h", "2h", "24h"), levels = c("2h", "24h"))) -> pl.df
pl.df %>%
  ggplot(aes(nr.sox2, delta, fill = nr.sox2)) +
  geom_boxplot(outlier.shape = NA) +
  ggpubr::stat_compare_means(method = "wilcox", comparisons = list(c("0", "1"), c("1", "2")), label = "p.signif", label.y = c(14,21), tip.length = 0) +
  facet_wrap(~time, ncol = 2) +
  scale_fill_manual(values = colorspace::sequential_hcl(n = 7, palette = "Blues3", rev = TRUE)[c(3,5,7)], breaks = 0:2) +
  ylab("CA frequency % (2h - untreated)") + xlab("Nr Sox2 motifs") +
  scale_y_continuous(breaks = c(-25,0,25), limits = c(-40,40)) +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3g_", Sys.Date(), ".pdf"), pl, width = 3, height = 4.5)

# H
partition.collapsing.dictionary = split(1:12,1:12)[c(2,11,6,5,9,10,4,7,12,3,1,8)]
patch.single.site.plots(
  interpretable.master.table = data.frame(Sample = "Sox2_kd_NO_", TFBS.cluster = "GenomicTile_45319426"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictionary, 
  data.type = "Sox2_kd_bait.capture", remove.TFBS.labels = TRUE, reutrn.chromatin.influence.df = TRUE,
  plotting.TFBSs = plyranges::filter(TFBS.cluster.compositions, cluster.id == plyranges::filter(TFBS.cluster.compositions, absolute.idx == "TFBS_8267068")$cluster.id)
) -> pl
CA_loci_df_dTAG %>% filter(TFBS.cluster == "GenomicTile_45319426") # 56%-33% | 307bp-296bp
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3h_", Sys.Date(), ".png"), width = 30, height = 30, units = "cm", res = 300)
pl$pl
dev.off()