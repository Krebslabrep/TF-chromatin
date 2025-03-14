library(doParallel)
library(foreach)
library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(magrittr)
library(GGally)
source("./scripts/functions/utils.r")
source("./scripts/functions/source_FootprintCharter_functions.r")
source("./scripts/functions/IGV_plotting.r")
detach("package:plyranges")

GenomicTiles = Load.GenomicTiles(tiles.width = 80)
chromHMM = Load.chromHMM(GenomicTiles = TRUE)
SNPs = LoadSNPs()
CytosinesToMask = Load.CytosinesToMask()

CA_loci_df_f1 = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-08-22_chromatin.influence_F1_PooledReplicates.df.qs") %>%
  process.CA.df_f1(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>% filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

CA_loci_df_dTAG = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2025-02-24_chromatin.influence_Sox2_dTAG_run2_run3_PooledReps.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

# A
CA_loci_df_f1 %>%
  dplyr::select(Sample, TFBS.cluster, chromHMM_annotation) %>%
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
  scale_y_continuous(breaks = c(0,medians,10,20,30,40)) +
  xlab(NULL) + ylab("Nr. SNPs per CRE") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = "none") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3a_", Sys.Date(), ".pdf"), pl, width = 5, height = 1.5)

# B
CA_loci_df_f1 %>%
  filter(ChIP_annotation == "bound") %>%
  filter(motif.change == "l.o.f.") %>%
  group_by(TF) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>% mutate(TF = factor(TF, levels = unique(TF))) -> pl.df
pl.df %>% ggplot() +
  geom_bar(aes(TF, n), color = "black", stat = "identity") +
  xlab(NULL) + ylab("allele-specific motifs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.5), text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3b_", Sys.Date(), ".pdf"), pl, height = 4, width = 5)

# data for C-I
# see /g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/discarded_scripts/supplementary_figure3_21.02.23.R
# sampleSheet = "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_all_bait_capture_data.txt"
# MySamples = readr::read_delim(sampleSheet, "\t")$SampleName %>% unique() %>% grep("_DE_", ., value = TRUE)
# Regions.of.Interest <- GRanges(Biostrings::seqinfo(BSgenome.Mmusculus.UCSC.mm10))[1:19]
# seqlevels(Regions.of.Interest) <- seqlevels(Regions.of.Interest)[1:19]
# # Parallelize the computation for each region of interest
# Reduce(c,
# parallel::mclapply(seq_along(Regions.of.Interest), function(i){
#   
#   RegionOfInterest <- Regions.of.Interest[i]
#   SingleMoleculeFootprinting::CallContextMethylation(
#     sampleSheet = sampleSheet, sample = MySamples, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, 
#     RegionOfInterest = RegionOfInterest, coverage = 20, ConvRate.thr = NULL, 
#     returnSM = FALSE, clObj = NULL
#   )
#   
# }, mc.preschedule = TRUE, mc.cores = 19)) -> Methylation
# 
# saveRDS(Methylation, "/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/Methylation.RDS")
Methylation = readRDS("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/Methylation.RDS")
Methylation = MaskSNPs2(Methylation, CytosinesToMask, MaskSMmat = FALSE, Experiment = "DE")

# C-D
Methylation %>%
  plyranges::select(grep("_Coverage", colnames(elementMetadata(.)), value = TRUE)) %>%
  as_tibble() %>%
  unite("coord", c("seqnames", "start")) %>%
  dplyr::select(-c(end, width, strand)) %>%
  gather(Sample, Coverage, -coord) %>%
  filter(!is.na(Coverage)) %>%
  mutate(Sample = gsub("_.*Barzaghi|_dp_rm.*$|TKO|I|V|Sample|1|2", "", Sample)) %>%
  group_by(coord, Sample) %>%
  summarise(Coverage = sum(Coverage, na.rm = TRUE), .groups = "drop") %>%
  spread(Sample, Coverage) %>%
  dplyr::select(-coord) -> Coverage_mat
Coverage_mat = log10(Coverage_mat)

png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3c-d_", Sys.Date(), ".png"), width = 14, height = 14, units = "in", res = 300)
pairs(Coverage_mat[,c("R_C", "A_C", "R_S", "A_S")], upper.panel = panel.cor, diag.panel = panel.hist, lower.panel = panel.jet, labels = c("R_C", "A_C", "R_S", "A_S"))
dev.off()

# E-G
Methylation %>%
  elementMetadata() %>%
  as.data.frame() %>%
  dplyr::select(grep("_MethRate", colnames(.))) -> MethRate_mat
colnames(MethRate_mat) = gsub("_.*Barzaghi|_dp_rm.*$|Sample1|Sample2|TKO", "", colnames(MethRate_mat))
MethRate_mat = MethRate_mat[,c(
  "R_CI", "A_CI", "R_CII", "A_CII", "R_CIII", "A_CIII", "R_CIV", "A_CIV", "R_SI", "A_SI", 
  "R_SII", "A_SII", "R_SIII", "A_SIII", "R_SIV", "A_SIV", "R_SVI", "A_SVI", "R_SVII", "A_SVII"
  )]

pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3e_", Sys.Date(), ".pdf"), width = 14, height = 14)
pairs(MethRate_mat[,c("R_CI", "R_CII")], upper.panel = panel.cor, diag.panel = panel.hist, lower.panel = panel.jet, labels = c("R_CI", "R_CII"))
dev.off()
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3f_", Sys.Date(), ".pdf"), width = 14, height = 14)
pairs(MethRate_mat[,c("A_CI", "A_CII")], upper.panel = panel.cor, diag.panel = panel.hist, lower.panel = panel.jet, labels = c("A_CI", "A_CII"))
dev.off()
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3g_", Sys.Date(), ".pdf"), width = 14, height = 14)
pairs(MethRate_mat[,c("A_SI", "A_SII")], upper.panel = panel.cor, diag.panel = panel.hist, lower.panel = panel.jet, labels = c("A_SI", "A_SII"))
dev.off()

# H-I
MethRate_mat %>%
  mutate(idx = seq(nrow(.))) %>%
  gather(sample, methrate, -idx) %>%
  separate(sample, into = c("allele", "sample")) %>%
  spread(allele, methrate) %>%
  mutate(delta = A-R) %>%
  select(-c(A,R)) %>%
  spread(sample, delta) %>%
  select(-idx) -> MethRate_mat_delta_delta

pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3h_", Sys.Date(), ".pdf"), width = 14, height = 14)
MethRate_mat[,grep("SVI|SVII", colnames(MethRate_mat), value = TRUE)] %>%
  transmute(SVI_delta = A_SVI - R_SVI, SVII_delta = A_SVII - R_SVII) %>%
  pairs(., upper.panel = panel.cor, diag.panel = panel.hist, lower.panel = panel.jet, labels = c("SVI", "SVII"))
dev.off()
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3i_", Sys.Date(), ".pdf"), width = 14, height = 14)
MethRate_mat[,grep("CIII|CIV", colnames(MethRate_mat), value = TRUE)] %>%
  transmute(CIII_delta = A_CIII - R_CIII, CIV_delta = A_CIV - R_CIV) %>%
  pairs(., upper.panel = panel.cor, diag.panel = panel.hist, lower.panel = panel.jet, labels = c("CIII", "CIV"))
dev.off()

# J
partition.collapsing.dictionary = split(1:12,1:12)[c(10,11,12,9,8,1,6,2,7,5,3,4)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_41709602", Sample = "CTKO"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictionary, 
  data.type = "F1_bait.capture", remove.TFBS.labels = TRUE
) -> pl
log2(TFBSs["TFBS_7406270"]$Cast.absScore/TFBSs["TFBS_7406270"]$BL6.absScore) # -2.065426
max(pl$chromatin.influence.df$acc.width.distro[[1]]) # 198bp
pl$chromatin.influence.df %>%
  mutate(CA.R = (accessible.R/tot.R)*100, CA.A = (accessible.A/tot.A)*100) %>%
  unnest(c(acc.width.distro, acc.read.count.distro.R, acc.read.count.distro.A)) %>%
  filter(acc.width.distro >= 100) %>%
  group_by(TF.name, TF, tot.R, tot.A, CA.R, CA.A) %>%
  summarise(CA.wide.R = sum(acc.read.count.distro.R), CA.wide.A = sum(acc.read.count.distro.A), .groups = "drop") %>%
  mutate(CA.wide.R = CA.wide.R/tot.R*100, CA.wide.A = CA.wide.A/tot.A*100) %>%
  mutate(CA.delta = CA.wide.A - CA.wide.R) # R: 86%, A: 26%
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3j_", Sys.Date(), ".png"), width = 30, height = 30, units = "cm", res = 300)
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
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3j_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()

# K
partition.collapsing.dictionary = split(1:12,1:12)[c(10,9,3,2,7,5,4,1,11,6,8)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_10204442", Sample = "CTKO"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictionary, 
  data.type = "F1_bait.capture", remove.TFBS.labels = TRUE, reutrn.chromatin.influence.df = FALSE
) -> pl
# R: 94%
# A: 90%
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3k_", Sys.Date(), ".png"), width = 30, height = 30, units = "cm", res = 300)
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
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3k_tracks_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()

# L
CA_loci_df_f1 %>%
  filter(TF == "Ctcf", ChIP_annotation == "bound") %>%
  mutate(freq_delta = CA_A - CA_R) %>%
  dplyr::select(Sample, TF.name, motif.change, freq_delta) %>%
  mutate(motif.change = factor(motif.change, levels = c("no", "l.o.f."))) %>%
  ggplot(aes(motif.change, freq_delta)) +
  geom_boxplot(fill = "blue", alpha = .8) +
  ggpubr::stat_compare_means(method = "wilcox", method.args = list(alternative = "less"), label.x = 1.5, label = "p.format") +
  scale_y_continuous(breaks = c(-75,-50,0,50), limits = c(-75,50)) +
  scale_x_discrete(labels = c(str_wrap("unaltered CTCF motifs", width = 10), str_wrap("allele-specific CTCF motifs", width = 10))) +
  xlab("") + ylab("CA frequency % (alt - ref)") +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3l_", Sys.Date(), ".pdf"), pl, width = 3.25, height = 5)
    
# M
CA_loci_df_f1 %>%
  filter(motif.change == "l.o.f.", ChIP_annotation == "bound", TF != "Ctcf") %>%
  mutate(ref = CA_regulatory_R, alt = CA_regulatory_A, delta = CA.delta) %>%
  arrange(ref) %>%
  mutate(rank = seq(nrow(.))) -> pl.df
pl.df %>%
  dplyr::select(rank, delta) %>%
  ggplot(aes(rank, delta, color = "delta")) +
  geom_point(size = .25, alpha = .1) +
  geom_smooth(method = "gam") +
  geom_hline(yintercept = -10) +
  xlab("motifs rank by ref CA frequency") + ylab("CA frequency % (alt - ref)") +
  scale_color_manual(values = c("blue"), breaks = c("delta")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(-75,-50,-10,0,50), limits = c(-75,50)) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = "none") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3m_", Sys.Date(), ".pdf"), pl, width = 3.5, height = 4.5)

# P
CA_loci_df_dTAG %>%
  filter(TF == "Sox2") %>%
  dplyr::select(TF.name, Sample, CA_regulatory) %>%
  spread(Sample, CA_regulatory) %>%
  mutate(dTAG_delta = Sox2_2h - Sox2_NT) %>%
  dplyr::select(-c(Sox2_2h, Sox2_NT)) %>%
  full_join(
    CA_loci_df_f1 %>%
      filter(TF == "Sox2", motif.change == "l.o.f.") %>%
      dplyr::rename("F1_delta" = "CA.delta") %>%
      dplyr::select(TF.name, F1_delta),
    by = "TF.name", multiple = "all"
  ) %>%
  na.omit() %>%
  ggplot(aes(F1_delta, dTAG_delta)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "grey", linetype = 2) +
  coord_fixed() +
  ylab("CA frequency % (2h - untreated)") +
  xlab("CA frequency % (alt - ref)") +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3p_", Sys.Date(), ".pdf"), pl, width = 3.5, height = 2.5)

# Q
CA_loci_df_dTAG %>%
  filter(TF == "Sox2", ChIP_annotation == "bound") %>%
  dplyr::select(-c(chromHMM_annotation, ChIP_annotation, tot.read.count, width, CA, CA_regulatory_count, width_regulatory)) %>%
  spread(Sample, CA_regulatory) %>%
  dplyr::rename("untreated" = "Sox2_NT", "2h dTAG" = "Sox2_2h") %>%
  na.omit() %>%
  mutate(delta = `2h dTAG` - untreated) %>%
  arrange(untreated) %>%
  mutate(rank = seq(nrow(.))) -> pl.df
pl.df %>%
  dplyr::select(rank, delta) %>%
  ggplot(aes(rank, delta, color = "delta")) +
  geom_point(size = .25, alpha = .1) +
  geom_smooth(method = "gam") +
  geom_hline(yintercept = -10) +
  xlab("motifs rank by untreated CA frequency") + ylab("CA frequency % (2h - untreated)") +
  scale_color_manual(values = c("blue"), breaks = c("delta")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(-75,-50,-10,0,50), limits = c(-75,50)) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = "none") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs3q_", Sys.Date(), ".pdf"), pl, width = 3.5, height = 4.5)
