library(doParallel)
library(foreach)
library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(magrittr)
library(GGally)
source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source('/g/krebs/barzaghi/Rscripts/useful_functionsV1.r')
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
source("/g/krebs/barzaghi/Rscripts/IGV_plotting.R")
detach("package:plyranges")

setwd("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/")

chromatin.influence.df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-23_chromatin.influence_F1_PooledReplicates.df.qs")
GenomicTiles = Load.GenomicTiles(tiles.width = 80)
chromHMM = Load.chromHMM(GenomicTiles = TRUE)
SNPs = LoadSNPs()

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
  mutate(CA.R = (accessible.R/tot.R)*100, CA.A = (accessible.A/tot.A)*100, CA.delta = CA.A - CA.R) -> CA_loci_df

left_join(
  chromatin.influence.df,
  chromHMM,
  by = "TFBS.cluster"
) %>%
  mutate(
    CA.R = (accessible.R/tot.R)*100, CA.A = (accessible.A/tot.A)*100, CA.delta = CA.A - CA.R,
    locus = ifelse(
      TF == "Ctcf", "Ctcf", 
      ifelse(chromHMM_annotation == "C4_CREs_promoter", "promoter", 
             ifelse(chromHMM_annotation == "C5_CREs_enhancer", "enhancer", 
                    ifelse(chromHMM_annotation == "C2_Inactive", "inactive", NA))))
  ) %>%
  filter(!is.na(locus)) %>%
  mutate(locus = factor(locus, levels = c("inactive", "enhancer", "promoter", "Ctcf"))) %>%
  filter(motif.change != "g.o.f.") %>%
  left_join(., ChIP_thresholds_dictionary, by = "TF") %>%
  unnest(c(acc.width.distro, acc.read.count.distro.R, acc.read.count.distro.A)) %>%
  filter(acc.width.distro >= 100) %>%
  group_by(Sample, TFBS.cluster, TF.name, motif.change, TF, tot.R, tot.A, CA.R, CA.A, locus, ChIP, ChIP_threshold) %>%
  summarise(CA.wide.R = sum(acc.read.count.distro.R), CA.wide.A = sum(acc.read.count.distro.A), .groups = "drop") %>%
  mutate(CA.wide.R = CA.wide.R/tot.R*100, CA.wide.A = CA.wide.A/tot.A*100, CA.delta = CA.wide.A-CA.wide.R) %>%
  filter(!is.na(ChIP)) %>%
  mutate(isBound = ChIP>ChIP_threshold) %>%
  filter(isBound) %>%
  mutate(TF = ifelse(TF == "Smad2::Smad3::Smad4", "Smad", ifelse(TF == "Bach1::Mafk", "Mafk", TF))) -> CA_loci_df

# A
CA_loci_df %>%
  dplyr::select(Sample, TFBS.cluster, locus) %>%
  mutate(
    cast = plyranges::count_overlaps(IRanges::resize(GenomicTiles[TFBS.cluster], 300, "center"), SNPs$Cast),
    spret = plyranges::count_overlaps(IRanges::resize(GenomicTiles[TFBS.cluster], 300, "center"), SNPs$Spret)
  ) %>%
  mutate(count = ifelse(Sample == "CTKO", cast, spret)) %>%
  dplyr::select(-cast, -spret)-> pl.df
medians = c(median(filter(pl.df, Sample == "CTKO")$count), median(filter(pl.df, Sample == "STKO")$count))

pl.df %>%
  mutate(Sample = factor(ifelse(Sample == "CTKO", "Bl6 x Castaneus", "Bl6 x Spretus"), levels = c("Bl6 x Spretus", "Bl6 x Castaneus"))) %>%
  ggplot(aes(Sample, count, fill = Sample)) +
  geom_boxplot() +
  scale_fill_manual(values = c("sienna", "grey45"), breaks = c("Bl6 x Castaneus", "Bl6 x Spretus")) +
  scale_y_continuous(breaks = c(0,medians, 10, 20, 30)) +
  xlab(NULL) + ylab("Nr. SNPs per CRE") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = "none") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3a_", Sys.Date(), ".pdf"), pl, width = 5, height = 1.5)

# B
CA_loci_df %>%
  filter(motif.change == "no") %>%
  dplyr::select(Sample, TFBS.cluster, locus, TF.name, CA.delta) %>%
  distinct() %>%
  mutate(Sample = factor(ifelse(Sample == "CTKO", "Bl6 x Castaneus", "Bl6 x Spretus"), levels = c("Bl6 x Castaneus", "Bl6 x Spretus"))) %>%
  ggplot(aes(Sample, CA.delta, fill = Sample)) +
  geom_boxplot() +
  scale_fill_manual(values = c("sienna", "grey45"), breaks = c("Bl6 x Castaneus", "Bl6 x Spretus")) +
  xlab(NULL) +
  ylab("CA delta (%)") +
  theme_bw() + 
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), legend.position = "none") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3b_", Sys.Date(), ".pdf"), pl, height = 5, width = 1.5)

# C
CA_loci_df %>%
  filter(!is.na(ChIP)) %>%
  mutate(isBound = ChIP>ChIP_threshold) %>%
  filter(isBound) %>%
  mutate(TF = ifelse(TF == "Smad2::Smad3::Smad4", "Smad", ifelse(TF == "Bach1::Mafk", "Mafk", TF))) %>%
  filter(motif.change == "l.o.f.") %>%
  group_by(TF) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>% mutate(TF = factor(TF, levels = unique(TF))) -> pl.df
pl.df$n %>% min
pl.df$n %>% max
pl.df %>% ggplot() +
  geom_bar(aes(TF, n), color = "black", stat = "identity") +
  xlab(NULL) + ylab("allele-specific motifs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.5), text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3c_", Sys.Date(), ".pdf"), pl, height = 4, width = 5)

# D-L: see /g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/discarded_scripts/supplementary_figure3_21.02.23.R

# B
CA_loci_df %>%
  filter(TF == "Ctcf", tot.R >= 100, tot.A >= 100, motif.change == "l.o.f.") %>%
  # filter(Sample == "STKO") %>% I tried until rank 11. Nothing useful
  # arrange(CA.delta) -> ctcf.pl.df
  filter(TFBS.cluster == "GenomicTile_41709602", Sample == "CTKO") -> ctcf.pl.df #  GenomicTile_32578421 [STKO, partial]
log2(TFBSs["TFBS_7406270"]$Cast.absScore/TFBSs["TFBS_7406270"]$BL6.absScore) # -2.065426

partition.collapsing.dictionary = split(1:12,1:12)[c(10,11,12,9,8,1,6,2,7,5,3,4)]
patch.single.site.plots(
  interpretable.master.table = ctcf.pl.df, rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictionary, 
  data.type = "F1_bait.capture", remove.TFBS.labels = TRUE, patch.version = 2
) -> pl
max(pl$chromatin.influence.df$acc.width.distro[[1]]) # 198bp
pl$chromatin.influence.df %>%
  mutate(CA.R = (accessible.R/tot.R)*100, CA.A = (accessible.A/tot.A)*100) %>%
  unnest(c(acc.width.distro, acc.read.count.distro.R, acc.read.count.distro.A)) %>%
  filter(acc.width.distro >= 100) %>%
  group_by(TF.name, TF, tot.R, tot.A, CA.R, CA.A) %>%
  summarise(CA.wide.R = sum(acc.read.count.distro.R), CA.wide.A = sum(acc.read.count.distro.A), .groups = "drop") %>%
  mutate(CA.wide.R = CA.wide.R/tot.R*100, CA.wide.A = CA.wide.A/tot.A*100) %>%
  mutate(CA.delta = CA.wide.A - CA.wide.R) # R: 86%, A: 26%
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3b_", Sys.Date(), ".png"), width = 30, height = 30, units = "cm", res = 300)
pl$pl
dev.off()

RegionOfInterest = IRanges::resize(GenomicTiles["GenomicTile_41709602"], 9000, "center")
plot_genomic_track(
  # sampleSheet = "/g/krebs/barzaghi/analyses/nf-core_runs/130923_BL6_Spret_F1_ATAC/aln_merged_deduplicated/samplesheet_PooledReplicates.txt",
  # samples = "BL6_Spret_F1",
  sampleSheet = "/g/krebs/barzaghi/HTS/ATAC-seq/Heard_NatGen_2017/aln/Qinput_male.txt",
  samples = "male_R1",
  RegionOfInterest = RegionOfInterest,
  allelic = TRUE,
  tile.width = 200,
  tile.step = 25,
  max.y.lim = 10,
  color = c("black", "sienna"), 
  y.labs = c("Bl6 ATAC", "Cast ATAC"),
  delta = FALSE,
  normalise = FALSE
) -> ATAC_track
plot_genomic_track(
  sampleSheet = ChIP_data_dictionary$Ctcf,
  samples = "Ctcf",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 25,
  max.y.lim = 2000,
  color = "black", 
  y.labs = "Ctcf"
) -> Ctcf_track
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
p_final <- ATAC_track + Ctcf_track + MNase_track +
  plot_layout(ncol = 1, heights = c(1/2, 1/4, 1/4))
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3b_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()

# M
CA_loci_df %>%
  filter(TF == "Ctcf", motif.change == "no", tot.R > 100, tot.A > 100, CA.wide.R > 90, CA.wide.A > 90) %>%
  filter(TFBS.cluster == "GenomicTile_11157474", Sample == "CTKO") -> ctcf.pl.df
partition.collapsing.dictionary = split(1:12,1:12)[c(10,9,3,2,7,5,4,1,11,6,8)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_10204442", Sample = "CTKO"), 
  rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictionary, 
  data.type = "F1_bait.capture", remove.TFBS.labels = TRUE, patch.version = 2, reutrn.chromatin.influence.df = FALSE
) -> pl
ctcf.pl.df$CA.wide.R # 94%
ctcf.pl.df$CA.wide.A # 90%
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3m_", Sys.Date(), ".png"), width = 30, height = 30, units = "cm", res = 300)
pl$pl
dev.off()

RegionOfInterest = IRanges::resize(GenomicTiles["GenomicTile_10204442"], 9000, "center")
plot_genomic_track(
  # sampleSheet = "/g/krebs/barzaghi/analyses/nf-core_runs/130923_BL6_Spret_F1_ATAC/aln_merged_deduplicated/samplesheet_PooledReplicates.txt",
  # samples = "BL6_Spret_F1",
  sampleSheet = "/g/krebs/barzaghi/HTS/ATAC-seq/Heard_NatGen_2017/aln/Qinput_male.txt",
  samples = "male_R1",
  RegionOfInterest = RegionOfInterest,
  allelic = TRUE,
  tile.width = 1000,
  tile.step = 100,
  max.y.lim = 30,
  color = "black", 
  y.labs = c("Bl6 ATAC", "Cast ATAC"),
  delta = FALSE,
  normalise = FALSE
) -> ATAC_track
plot_genomic_track(
  sampleSheet = ChIP_data_dictionary$Ctcf,
  samples = "Ctcf",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 10,
  max.y.lim = 2000,
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
  y.labs = "MNAse"
) -> MNase_track
p_final <- ATAC_track + Ctcf_track + MNase_track +
  plot_layout(ncol = 1, heights = c(1/2, 1/4, 1/4))
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3m_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()







