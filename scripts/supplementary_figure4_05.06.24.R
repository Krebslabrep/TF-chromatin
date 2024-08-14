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

chromatin.influence.df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-06-04_chromatin.influence_F1_PooledReplicates.df.qs")
chromHMM = Load.chromHMM(GenomicTiles = TRUE)
GenomicTiles = Load.GenomicTiles(tiles.width = 80, GeneAssociation = FALSE, CGI = FALSE)
TFBSs = Load.TFBSs()

cluster <- multidplyr::new_cluster(16)
chromatin.influence.df %>%
  filter(motif.change != "g.o.f.") %>%
  left_join(., chromHMM, by = "TFBS.cluster") %>%
  mutate(
    locus = factor(case_when(
      TF == "Ctcf" ~ "Ctcf",
      chromHMM_annotation == "C4_CREs_promoter" ~ "promoter",
      chromHMM_annotation == "C5_CREs_enhancer" ~ "enhancer",
      chromHMM_annotation == "C2_Inactive" ~ "inactive"
    ), levels = c("inactive", "enhancer", "promoter", "Ctcf"))) %>%
  left_join(., ChIP_thresholds_dictionary_lenient, by = "TF") %>%
  filter(!is.na(ChIP)) %>%
  mutate(isBound = case_when(
    ChIP >= ChIP_threshold ~ "strong",
    ChIP < ChIP_threshold & ChIP > 5 ~ "weak",
    ChIP <= 5 ~ "unbound"
  )) %>%
  filter(isBound == "strong") %>% #
  mutate(TF = ifelse(TF == "Smad2::Smad3::Smad4", "Smad", ifelse(TF == "Bach1::Mafk", "Mafk", TF))) %>%
  mutate(TF = ifelse(motif.change == "no", "unperturbed", as.character(TF))) %>%
  group_by(TF, locus, motif.change) %>%
  unnest(c(acc.width.distro, acc.read.count.distro.R, acc.read.count.distro.A)) %>%
  filter(acc.width.distro >= 100) %>%
  group_by(Sample, TF.name, TF, locus, motif.change, isBound, tot.R, tot.A) %>%
  summarise(R = sum(acc.read.count.distro.R), A = sum(acc.read.count.distro.A), .groups = "drop") %>%
  rowwise() %>%
  multidplyr::partition(cluster) %>%
  mutate(pval = fisher.test(matrix(c(A,tot.A-A,R,tot.R-R),ncol = 2), alternative = "two.sided")$p.value) %>%
  ungroup() %>%
  collect() %>%
  mutate(CA.delta = A/tot.A*100 - R/tot.R*100) -> CA_loci_df
remove(cluster)
CA_loci_df$p.adj = NA
CA_loci_df$p.adj[CA_loci_df$motif.change=="l.o.f."] = p.adjust(p = filter(CA_loci_df, motif.change=="l.o.f.")$pval, method = "BH")
CA_loci_df %<>% mutate(change = ifelse(p.adj <= 0.05 & CA.delta < 0, "down", ifelse(p.adj <= 0.05 & CA.delta > 0, "up", "no")))

# Check if the percentage of change=="down" changes with isBound categories. Can I use the weak as well?
# CA_loci_df %>%
#   filter(motif.change == "l.o.f.") %>%
#   mutate(p.adj = ifelse(-log10(p.adj)>25, 0, p.adj)) %>%
#   ggplot() +
#   geom_point(aes(CA.delta, -log10(p.adj), color = change)) +
#   facet_grid(isBound~locus) +
#   scale_color_manual(values = colorspace::diverge_hcl(palette = "Berlin", n = 3), breaks = c("up", "no", "down")) +
#   theme_bw()
# 
# CA_loci_df %>%
#   filter(motif.change == "l.o.f.") %>%
#   group_by(locus, isBound) %>% summarise(up = sum(change == "up") / n() * 100, no = sum(change == "no") / n() * 100, down = sum(change == "down") / n() * 100, .groups = "drop") %>%
#   gather(change, p, up, no, down) %>%
#   mutate(change = factor(change, levels = c("up", "no", "down"))) %>%
#   mutate(isBound = factor(isBound, levels = c("strong", "weak", "unbound"))) %>%
#   ggplot() +
#   geom_col(aes(isBound, p, fill = change)) +
#   scale_fill_manual(values = colorspace::diverge_hcl(palette = "Berlin", n = 3), breaks = c("up", "no", "down")) +
#   facet_grid(~locus) +
#   xlab(NULL) + ylab(NULL) +
#   theme_bw()
# 
# CA_loci_df %>%
#   filter(motif.change == "l.o.f.") %>%
#   mutate(isBound = factor(isBound, levels = c("unbound", "weak", "strong"))) %>%
#   ggplot() +
#   geom_density(aes(pval, color = isBound), size = 2, adjust = 1/10) +
#   scale_color_manual(values = viridis::inferno(n = 3, direction = -1), breaks = c("strong", "unbound", "weak")) +
#   scale_y_continuous(breaks = c(0,14)) +
#   scale_x_continuous(breaks = c(0,1)) +
#   xlab("p-value") +
#   theme_bw() +
#   theme(text = element_text(size = 18), legend.position = c(.5,.9), legend.title = element_blank(), legend.background = element_blank())

  
# A
CA_loci_df %>%
  mutate(motif.change = factor(ifelse(motif.change == "l.o.f.", "allele-specific", "no change"), levels = c("no change", "allele-specific"))) %>%
  ggplot() +
  geom_density(aes(pval, color = motif.change), size = 2, adjust = 1/10) +
  scale_color_manual(values = c("black", "salmon"), breaks = c("no change", "allele-specific")) +
  scale_y_continuous(breaks = c(0,7.5)) +
  scale_x_continuous(breaks = c(0,1)) +
  xlab("p-value") +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = c(.5,.9), legend.title = element_blank(), legend.background = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs4a_", Sys.Date(), ".pdf"), pl, width = 4, height = 4)

# B
CA_loci_df %>%
  filter(motif.change == "l.o.f.") %>%
  filter(change == "down") -> lof.motifs

SingleMoleculeFootprinting::Arrange_TFBSs_clusters(TFBSs[unique(lof.motifs$TF.name)]) -> lof.clusters
data.frame(
  TF.name = Reduce(c, lof.clusters$ClusterComposition)$absolute.idx,
  nr.motifs = rep(unname(lengths(lof.clusters$ClusterComposition)), lengths(lof.clusters$ClusterComposition))
) -> lof.cobinders

lof.motifs %>%
  filter(locus != "Ctcf") %>%
  left_join(., lof.cobinders, by = "TF.name") %>%
  mutate(nr.motifs = factor(ifelse(is.na(nr.motifs), 1, nr.motifs), levels = c(seq(max(nr.motifs, na.rm = TRUE))))) %>%
  # group_by(nr.motifs, locus) %>% filter(n()>10) %>% ungroup() %>%
  ggplot(aes(nr.motifs, CA.delta)) +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox", comparisons = list(c(1,2), c(2,3))) +
  geom_boxplot() +
  # facet_wrap(~locus, nrow = 1) +
  theme_bw() -> boxes
lof.motifs %>%
  filter(locus != "Ctcf") %>%
  left_join(., lof.cobinders, by = "TF.name") %>%
  mutate(nr.motifs = factor(ifelse(is.na(nr.motifs), 1, nr.motifs), levels = c(seq(max(nr.motifs, na.rm = TRUE))))) %>%
  # group_by(nr.motifs, locus) %>% filter(n()>10) %>% ungroup() %>%
  mutate(CA_R = R/tot.R*100, CA_A = A/tot.A*100) %>%
  dplyr::select(TF.name, locus, CA_R, CA_A, nr.motifs) %>%
  gather(Sample, CA, CA_R, CA_A) %>%
  ggplot(aes(nr.motifs, CA, fill = Sample)) +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox") +
  geom_boxplot() +
  # facet_wrap(~locus, nrow = 1) +
  theme_bw() -> boxes2
lof.motifs %>%
  filter(locus != "Ctcf") %>%
  left_join(., lof.cobinders, by = "TF.name") %>%
  mutate(nr.motifs = factor(ifelse(is.na(nr.motifs), 1, nr.motifs), levels = c(seq(max(nr.motifs, na.rm = TRUE))))) %>%
  # group_by(nr.motifs, locus) %>% filter(n()>10) %>% ungroup() %>%
  ggplot(aes(nr.motifs)) +
  geom_bar() +
  # facet_wrap(~locus, nrow = 1) +
  scale_y_log10() +
  theme_bw() -> bars
boxes / boxes2 / bars + plot_layout(axes = "collect_x")
