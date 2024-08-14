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
chromatin.influence.df_dTAG = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-30_chromatin.influence_Oct4_Sox2_KD.df.qs")
chromHMM = Load.chromHMM(GenomicTiles = TRUE)
GenomicTiles = Load.GenomicTiles(tiles.width = 80, GeneAssociation = FALSE, CGI = FALSE)
TFBSs = Load.TFBSs()

CA_loci_df = process.CA.df_f1(chromatin.influence.df, cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient)
# qs::qsave(filter(CA_loci_df, motif.change == "l.o.f.", ChIP_annotation %in% c("bound")), paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/f1_lofs_smf_", Sys.Date(), ".qs"))
CA_loci_df_dTAG = process.CA.df(chromatin.influence.df_dTAG, cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient)

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

# C
table(
  filter(CA_loci_df, ChIP_annotation == "bound")$motif.change, 
  filter(CA_loci_df, ChIP_annotation == "bound")$change
  )

CA_loci_df %>%
  filter(motif.change == "l.o.f.", ChIP_annotation == "bound") %>%
  filter(change == "down") %>%
  group_by(TF) %>% filter(n() > 10) %>% ungroup() %>%
  group_by(TF) %>% mutate(TF = paste0(TF, " (", n(), ")")) %>% ungroup() %>%
  group_by(TF) %>% mutate(m = median(CA.delta, na.rm=TRUE)) %>% ungroup() %>% arrange(desc(m)) %>% mutate(TF = factor(TF, levels = unique(TF))) -> pl.df
pl.df %>% ggplot(aes(TF, CA.delta)) +
  geom_hline(yintercept = c(0, median(filter(CA_loci_df, motif.change == "no")$CA.delta)), color = c("grey", "red"), linetype = 2) +
  geom_boxplot() +
  xlab(NULL) + ylab(paste0(expression(Delta), " CA (%)")) +
  ylim(c(-80, 10)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.5), text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3c_", Sys.Date(), ".pdf"), pl, width = 4, height = 4)
distinct(pl.df, TF, m) # 14%, 23%, 30%
quantile(filter(pl.df, str_detect(TF, "Ctcf", negate = TRUE))$CA.delta)

# table(filter(CA_loci_df, motif.change == "l.o.f.")$change)
# table(filter(CA_loci_df, motif.change == "l.o.f.")$change) / nrow(filter(CA_loci_df, motif.change == "l.o.f.")) * 100
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3c_inlet_", Sys.Date(), ".pdf"), width = 5, height = 5)
pie(x = table(filter(CA_loci_df, motif.change == "l.o.f.", ChIP_annotation == "bound")$change), labels = c("allele-specific", "no change"), col = c("salmon", "black"))
dev.off()

# D
# inactive.median = as.integer(quantile(mutate(filter(CA_loci_df, motif.change == "l.o.f.", locus == "inactive"), CA.R = R/tot.R*100)$CA.R, 0.50))
# pl.df %>%
#   mutate(color = ifelse(str_detect(TF, "Ctcf"), "Ctcf", "other TFs")) %>%
#   ggplot(aes(R/tot.R*100, A/tot.A*100, color = color)) +
#   geom_abline() +
#   geom_point(size = 2) +
#   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = inactive.median, fill = viridis::mako(n=1)[1], alpha = .25) +
#   scale_color_manual(values = c("green3", "black"), breaks = c("Ctcf", "other TFs")) + 
#   scale_x_continuous(breaks = c(0,100), limits = c(0,100)) + 
#   scale_y_continuous(breaks = c(0,inactive.median,100), limits = c(0,100)) + 
#   xlab("CA at reference allele (%)") +
#   ylab("CA at alternative allele (%)") +
#   coord_fixed() +
#   theme_bw() +
#   theme(text = element_text(size = 18), legend.position=c(.25, .9), legend.title = element_blank(), legend.background = element_blank()) -> pl
# ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3d_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5)

inactive.median = as.integer(quantile(filter(CA_loci_df, motif.change == "l.o.f.", chromHMM_annotation == "intergenic", ChIP_annotation == "unbound")$CA_R, 0.50))
pl.df %>%
  mutate(ref = CA_regulatory_R, alt = CA_regulatory_A) %>%
  arrange(ref) %>%
  mutate(rank = seq(nrow(.))) %>%
  dplyr::select(rank, ref, alt) %>%
  gather(allele, CA, ref, alt) %>%
  ggplot(aes(rank, CA, color = allele)) +
  geom_point(size = 1, alpha = .5) +
  geom_smooth(se = TRUE, span = .5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = inactive.median, fill = viridis::mako(n=1)[1], alpha = .25) +
  xlab(NULL) + ylab("CA (%)") +
  scale_color_manual(values = c("black", "salmon"), breaks = c("ref", "alt")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(0,inactive.median,50,100), limits = c(0,100)) + 
  theme_bw() +
  theme(text = element_text(size = 18), legend.position=c(.2, .875), legend.title = element_blank(), legend.background = element_blank()) -> up
pl.df %>%
  mutate(ref = CA_regulatory_R, alt = CA_regulatory_A, delta = CA.delta) %>%
  arrange(ref) %>%
  mutate(rank = seq(nrow(.))) %>%
  dplyr::select(rank, delta) %>%
  ggplot(aes(rank, delta, color = "delta")) +
  geom_point(size = 1, alpha = .5) +
  geom_smooth(se = TRUE, span = 1) +
  xlab("CA rank") + ylab("CA delta (%)") +
  scale_color_manual(values = c("blue"), breaks = c("delta")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(-80,-40,0), limits = c(-80,0)) +
  theme_bw() +
  theme(text = element_text(size = 18), plot.margin = unit(c(0,0,0,0), "cm")) -> down
up / down + patchwork::plot_layout(heights = c(10/18, 8/18), axes = "collect") -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3d_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5)

# E
# left-to-right: FOXD3, STAT3, KLF4, KLF4
partition.collapsing.dictionary = split(1:12,1:12)[c(12,7,6,1,5,2,4,11,3,9,10,8)]
patch.single.site.plots(
  interpretable.master.table = data.frame(Sample = "CTKO", TFBS.cluster = "GenomicTile_17185595"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictionary, 
  data.type = "F1_bait.capture", remove.TFBS.labels = TRUE, patch.version = 2
) -> pl
max(pl$chromatin.influence.df$acc.width.distro[[1]]) # 331bp
pl$chromatin.influence.df[1,] %>%
  mutate(CA.R = (accessible.R/tot.R)*100, CA.A = (accessible.A/tot.A)*100) %>%
  unnest(c(acc.width.distro, acc.read.count.distro.R, acc.read.count.distro.A)) %>%
  filter(acc.width.distro >= 100) %>%
  group_by(TF.name, TF, tot.R, tot.A, CA.R, CA.A) %>%
  summarise(CA.wide.R = sum(acc.read.count.distro.R), CA.wide.A = sum(acc.read.count.distro.A), .groups = "drop") %>%
  mutate(CA.wide.R = CA.wide.R/tot.R*100, CA.wide.A = CA.wide.A/tot.A*100) %>%
  mutate(CA.delta = CA.wide.A - CA.wide.R) # R: 50%, A: 30%
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3e_", Sys.Date(), ".png"), width = 30, height = 30, units = "cm", res = 300)
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

# F
cluster <- multidplyr::new_cluster(16)
CA_loci_df_dTAG %>%
  filter(str_detect(Sample, "Sox2") & TF == "Sox2") %>%
  mutate(Sample = gsub("Sox2|_NO_", "", Sample)) %>%
  dplyr::select(-TF, -CA) %>%
  pivot_wider(names_from = "Sample", id_cols = c("TFBS.cluster", "TF.name", "chromHMM_annotation"), values_from = c("CA_regulatory_count", "tot.read.count")) %>%
  na.omit() %>%
  rowwise() %>%
  multidplyr::partition(cluster) %>%
  mutate(pval = fisher.test(matrix(c(CA_regulatory_count_2h,tot.read.count_2h-CA_regulatory_count_2h,CA_regulatory_count_NT,tot.read.count_NT-CA_regulatory_count_NT),ncol = 2), alternative = "less")$p.value) %>%
  ungroup() %>%
  collect() %>%
  mutate(CA.delta = CA_regulatory_count_2h/tot.read.count_2h*100 - CA_regulatory_count_NT/tot.read.count_NT*100) -> dtag.tested
remove(cluster)
dtag.tested %<>% mutate(padj = p.adjust(p = pval, method = "BH"))

# dtag.tested %>%
#   ggplot() +
#   geom_histogram(aes(pval))

# dtag.tested %>%
#   filter(padj < .05) %>%
#   mutate(NT = CA_regulatory_count_NT/tot.read.count_NT*100, `2h` = CA_regulatory_count_2h/tot.read.count_2h*100) %>%
#   ggplot(aes(NT, `2h`)) +
#   geom_point(size = 2) +
#   geom_abline(linetype = 2, color = "grey") + 
#   scale_x_continuous(breaks = c(0,100), limits = c(0,100)) + 
#   scale_y_continuous(breaks = c(0,100), limits = c(0,100)) + 
#   xlab("Sox2 motifs, untreated (CA %)") +
#   ylab("Sox2 motifs, 2h dTAG (CA %)") +
#   coord_fixed() +
#   theme_bw() +
#   theme(text = element_text(size = 18), legend.position=c(.25, .9), legend.title = element_blank(), legend.background = element_blank()) -> pl
# ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig3f_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5)
# quantile(filter(dtag.tested, padj < .05)$CA.delta)

dtag.tested %>%
  filter(padj < .05) -> pl.df

pl.df %>%
  mutate(NT = CA_regulatory_count_NT/tot.read.count_NT*100, `2h` = CA_regulatory_count_2h/tot.read.count_2h*100) %>%
  arrange(NT) %>%
  mutate(rank = seq(nrow(.))) %>%
  dplyr::select(rank, NT, `2h`) %>%
  gather(condition, CA, NT, `2h`) %>%
  ggplot(aes(rank, CA, color = condition)) +
  geom_point(size = 1, alpha = .5) +
  geom_smooth(se = TRUE, span = .5) +
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = inactive.median, fill = viridis::mako(n=1)[1], alpha = .25) +
  xlab(NULL) + ylab("CA (%)") +
  scale_color_manual(values = c("black", "salmon"), breaks = c("NT", "2h")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(0,inactive.median,50,100), limits = c(0,100)) + 
  theme_bw() +
  theme(text = element_text(size = 18), legend.position=c(.2, .875), legend.title = element_blank(), legend.background = element_blank()) -> up
pl.df %>%
  mutate(NT = CA_regulatory_count_NT/tot.read.count_NT*100, `2h` = CA_regulatory_count_2h/tot.read.count_2h*100, delta = `2h` - NT) %>%
  arrange(NT) %>%
  mutate(rank = seq(nrow(.))) %>%
  dplyr::select(rank, delta) %>%
  ggplot(aes(rank, delta, color = "delta")) +
  geom_point(size = 1, alpha = .5) +
  geom_smooth(se = TRUE, span = 1) +
  xlab("CA rank") + ylab("CA delta (%)") +
  scale_color_manual(values = c("blue"), breaks = c("delta")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(-80,-40,0), limits = c(-80,0)) +
  theme_bw() +
  theme(text = element_text(size = 18), plot.margin = unit(c(0,0,0,0), "cm")) -> down
up / down + patchwork::plot_layout(heights = c(10/18, 8/18), axes = "collect") -> pl

chromatin.influence.df_chem = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-30_chromatin.influence_chemical_inhibitions.df.qs")
CA_loci_df_chem = process.CA.df(chromatin.influence.df_chem, cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient)

cluster <- multidplyr::new_cluster(16)
CA_loci_df_chem %>%
  filter(Sample %in% c("DMSO", "BRM"), TF == "Sox2", ChIP_annotation %in% c("bound", "weak")) %>%
  dplyr::select(-TF, -CA) %>%
  pivot_wider(names_from = "Sample", id_cols = c("TFBS.cluster", "TF.name", "chromHMM_annotation"), values_from = c("CA_regulatory_count", "tot.read.count")) %>%
  na.omit() %>%
  rowwise() %>%
  multidplyr::partition(cluster) %>%
  mutate(pval = fisher.test(matrix(c(CA_regulatory_count_BRM,tot.read.count_BRM-CA_regulatory_count_BRM,CA_regulatory_count_DMSO,tot.read.count_DMSO-CA_regulatory_count_DMSO),ncol = 2), alternative = "less")$p.value) %>%
  ungroup() %>%
  collect() %>%
  mutate(CA.delta = CA_regulatory_count_BRM/tot.read.count_BRM*100 - CA_regulatory_count_DMSO/tot.read.count_DMSO*100) -> chem.tested
remove(cluster)
chem.tested %<>% mutate(padj = p.adjust(p = pval, method = "BH"))

chem.tested %>%
  filter(padj < .05) -> pl.df

pl.df %>%
  mutate(DMSO = CA_regulatory_count_DMSO/tot.read.count_DMSO*100, BRM = CA_regulatory_count_BRM/tot.read.count_BRM*100) %>%
  arrange(DMSO) %>%
  mutate(rank = seq(nrow(.))) %>%
  dplyr::select(rank, DMSO, BRM) %>%
  gather(condition, CA, DMSO, BRM) %>%
  ggplot(aes(rank, CA, color = condition)) +
  geom_point(size = 1, alpha = .5) +
  geom_smooth(se = TRUE, span = .5) +
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = inactive.median, fill = viridis::mako(n=1)[1], alpha = .25) +
  xlab(NULL) + ylab("CA (%)") +
  scale_color_manual(values = c("black", "salmon"), breaks = c("DMSO", "BRM")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(0,inactive.median,50,100), limits = c(0,100)) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position=c(.2, .875), legend.title = element_blank(), legend.background = element_blank()) -> up
pl.df %>%
  mutate(DMSO = CA_regulatory_count_DMSO/tot.read.count_DMSO*100, BRM = CA_regulatory_count_BRM/tot.read.count_BRM*100, delta = BRM - DMSO) %>%
  arrange(DMSO) %>%
  mutate(rank = seq(nrow(.))) %>%
  dplyr::select(rank, delta) %>%
  ggplot(aes(rank, delta, color = "delta")) +
  geom_point(size = 1, alpha = .5) +
  geom_smooth(se = TRUE, span = 1) +
  xlab("CA rank") + ylab("CA delta (%)") +
  scale_color_manual(values = c("blue"), breaks = c("delta")) +
  scale_x_continuous(breaks = c(0,nrow(pl.df)), limits = c(0,nrow(pl.df))) +
  scale_y_continuous(breaks = c(-80,-40,0), limits = c(-80,0)) +
  theme_bw() +
  theme(text = element_text(size = 18), plot.margin = unit(c(0,0,0,0), "cm")) -> down
up / down + patchwork::plot_layout(heights = c(10/18, 8/18), axes = "collect") -> pl


full_join(
  dtag.tested %>% dplyr::select(TF.name, CA.delta, padj) %>% dplyr::rename("dTAG.delta" = "CA.delta", "dTAG.padj" = "padj"),
  chem.tested %>% dplyr::select(TF.name, CA.delta, padj) %>% dplyr::rename("chem.delta" = "CA.delta", "chem.padj" = "padj")
) %>%
  filter(dTAG.padj < .05 | chem.padj < .05) %>%
  ggplot(aes(dTAG.delta, chem.delta)) +
  ggpointdensity::geom_pointdensity() +
  geom_abline(color = "grey", linetype = 2) + geom_vline(xintercept = 0, color = "grey", linetype = 2) + geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  viridis::scale_colour_viridis() +
  theme_bw()

full_join(
  dtag.tested %>% dplyr::select(TF.name, CA.delta, padj) %>% dplyr::rename("dTAG.delta" = "CA.delta", "dTAG.padj" = "padj"),
  chem.tested %>% dplyr::select(TF.name, CA.delta, padj) %>% dplyr::rename("chem.delta" = "CA.delta", "chem.padj" = "padj")
) %>%
  filter(dTAG.padj < .05 | chem.padj < .05) %>%
  dplyr::select(-c(dTAG.padj, chem.padj)) %>%
  gather(experiment, CA.delta, dTAG.delta, chem.delta) %>%
  ggplot(aes(CA.delta, fill = experiment)) +
  geom_density(alpha = .5) +
  theme_bw()
















cluster <- multidplyr::new_cluster(16)
CA_loci_df_dTAG %>%
  filter(
    (str_detect(Sample, "Sox2") & TF == "Sox2" & ChIP_annotation == "bound") |
    (str_detect(Sample, "Sox2") & TF == "Oct4" & ChIP_annotation == "bound") |
    (str_detect(Sample, "Sox2") & TF == "Rest" & ChIP_annotation == "bound") |
    (str_detect(Sample, "Sox2") & TF == "Ctcf" & ChIP_annotation == "bound")
      ) %>% 
  mutate(Sample = gsub("Sox2|_NO_", "", Sample)) %>%
  dplyr::select(-CA) %>%
  pivot_wider(names_from = "Sample", id_cols = c("TF", "TFBS.cluster", "TF.name", "chromHMM_annotation"), values_from = c("CA_regulatory_count", "tot.read.count")) %>%
  na.omit() %>%
  rowwise() %>%
  multidplyr::partition(cluster) %>%
  mutate(pval = fisher.test(matrix(c(CA_regulatory_count_2h,tot.read.count_2h-CA_regulatory_count_2h,CA_regulatory_count_NT,tot.read.count_NT-CA_regulatory_count_NT),ncol = 2), alternative = "less")$p.value) %>%
  ungroup() %>%
  collect() %>%
  mutate(CA.delta = CA_regulatory_count_2h/tot.read.count_2h*100 - CA_regulatory_count_NT/tot.read.count_NT*100) -> dtag.tested
remove(cluster)

cluster <- multidplyr::new_cluster(16)
CA_loci_df_chem %>%
  filter(
    (TF == "Sox2" & ChIP_annotation == "bound") |
    (TF == "Oct4" & ChIP_annotation == "bound") |
    (TF == "Rest" & ChIP_annotation == "bound") |
    (TF == "Ctcf" & ChIP_annotation == "bound")
  ) %>% 
  dplyr::select(-CA) %>%
  pivot_wider(names_from = "Sample", id_cols = c("TF", "TFBS.cluster", "TF.name", "chromHMM_annotation"), values_from = c("CA_regulatory_count", "tot.read.count")) %>%
  na.omit() %>%
  rowwise() %>%
  multidplyr::partition(cluster) %>%
  mutate(pval = fisher.test(matrix(c(CA_regulatory_count_BRM,tot.read.count_BRM-CA_regulatory_count_BRM,CA_regulatory_count_DMSO,tot.read.count_DMSO-CA_regulatory_count_DMSO),ncol = 2), alternative = "less")$p.value) %>%
  ungroup() %>%
  collect() %>%
  mutate(CA.delta = CA_regulatory_count_BRM/tot.read.count_BRM*100 - CA_regulatory_count_DMSO/tot.read.count_DMSO*100) -> chem.tested
remove(cluster)

full_join(
  dtag.tested %>% dplyr::select(TF, chromHMM_annotation, TF.name, CA.delta, pval) %>% dplyr::rename("dTAG.delta" = "CA.delta", "dTAG.p" = "pval"),
  chem.tested %>% dplyr::select(TF, chromHMM_annotation, TF.name, CA.delta, pval) %>% dplyr::rename("chem.delta" = "CA.delta", "chem.p" = "pval")
) %>%
  # filter(dTAG.p < .05 | chem.p < .05) %>%
  ggplot(aes(dTAG.delta, chem.delta)) +
  ggpointdensity::geom_pointdensity() +
  geom_abline(color = "grey", linetype = 2) + geom_vline(xintercept = 0, color = "grey", linetype = 2) + geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  facet_wrap(~TF) +
  viridis::scale_colour_viridis() +
  scale_y_continuous(breaks = seq(-60,60,10), labels = seq(-60,60,10)) +
  scale_x_continuous(breaks = seq(-60,60,10), labels = seq(-60,60,10)) +
  coord_fixed() +
  theme_bw()

full_join(
  dtag.tested %>% dplyr::select(TF, TF.name, CA.delta, pval) %>% dplyr::rename("dTAG.delta" = "CA.delta", "dTAG.p" = "pval"),
  chem.tested %>% dplyr::select(TF, TF.name, CA.delta, pval) %>% dplyr::rename("chem.delta" = "CA.delta", "chem.p" = "pval"),
  by = join_by(TF, TF.name)
) %>%
  # filter(dTAG.p < .05 | chem.p < .05) %>%
  dplyr::select(-c(dTAG.p, chem.p)) %>%
  gather(experiment, CA.delta, dTAG.delta, chem.delta) %>%
  ggplot(aes(TF, CA.delta, fill = experiment)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  geom_violin(draw_quantiles = c(.25, .5, .75)) +
  ggpubr::stat_compare_means(method = "wilcox", label = "p.signif") +
  scale_y_continuous(breaks = seq(-60,60,10), labels = seq(-60,60,10)) +
  xlab(NULL) +
  theme_bw() +
  theme(text = element_text(size = 18))

full_join(
  dtag.tested %>% filter(tot.read.count_2h >= 100 & tot.read.count_NT >= 100) %>% dplyr::select(TF, chromHMM_annotation, TF.name, CA.delta, pval) %>% dplyr::rename("dTAG.delta" = "CA.delta", "dTAG.p" = "pval"),
  chem.tested %>% filter(tot.read.count_BRM >= 100 & tot.read.count_DMSO >= 100) %>% dplyr::select(TF, chromHMM_annotation, TF.name, CA.delta, pval) %>% dplyr::rename("chem.delta" = "CA.delta", "chem.p" = "pval")
) %>%
  filter(dTAG.p < .05 & chem.p < .05) %>%
  ggplot(aes(dTAG.delta, chem.delta)) +
  ggpointdensity::geom_pointdensity() +
  geom_abline(color = "grey", linetype = 2) + geom_vline(xintercept = 0, color = "grey", linetype = 2) + geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  facet_wrap(~TF) +
  viridis::scale_colour_viridis() +
  scale_y_continuous(breaks = seq(-60,60,10), labels = seq(-60,60,10)) +
  scale_x_continuous(breaks = seq(-60,60,10), labels = seq(-60,60,10)) +
  coord_fixed() +
  theme_bw()

full_join(
  dtag.tested %>% filter(tot.read.count_2h >= 100 & tot.read.count_NT >= 100) %>% dplyr::select(TF, TF.name, CA.delta, pval) %>% dplyr::rename("dTAG.delta" = "CA.delta", "dTAG.p" = "pval"),
  chem.tested %>% filter(tot.read.count_BRM >= 100 & tot.read.count_DMSO >= 100) %>% dplyr::select(TF, TF.name, CA.delta, pval) %>% dplyr::rename("chem.delta" = "CA.delta", "chem.p" = "pval"),
  by = join_by(TF, TF.name)
) %>%
  filter(dTAG.p < .05 & chem.p < .05) %>%
  dplyr::select(-c(dTAG.p, chem.p)) %>%
  gather(experiment, CA.delta, dTAG.delta, chem.delta) %>%
  ggplot(aes(TF, CA.delta, fill = experiment)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  geom_violin(draw_quantiles = c(.25, .5, .75)) +
  ggpubr::stat_compare_means(method = "wilcox", label = "p.signif") +
  scale_y_continuous(breaks = seq(-60,60,10), labels = seq(-60,60,10)) +
  xlab(NULL) +
  theme_bw() +
  theme(text = element_text(size = 18))

full_join(
  dtag.tested %>% filter(tot.read.count_2h >= 100 & tot.read.count_NT >= 100) %>% dplyr::select(TF, TF.name, CA.delta, pval) %>% dplyr::rename("dTAG.delta" = "CA.delta", "dTAG.p" = "pval"),
  chem.tested %>% filter(tot.read.count_BRM >= 100 & tot.read.count_DMSO >= 100) %>% dplyr::select(TF, TF.name, CA.delta, pval) %>% dplyr::rename("chem.delta" = "CA.delta", "chem.p" = "pval"),
  by = join_by(TF, TF.name)
) %>%
  filter(dTAG.p < .05 & chem.p < .05) %>%
  dplyr::select(-c(dTAG.p, chem.p)) %>%
  gather(experiment, CA.delta, dTAG.delta, chem.delta) %>%
  ggplot(aes(TF, fill = experiment)) +
  geom_bar(position = "dodge") +
  xlab(NULL) +
  theme_bw() +
  theme(text = element_text(size = 18))




full_join(
  dtag.tested %>% filter(tot.read.count_2h >= 100 & tot.read.count_NT >= 100) %>% dplyr::select(TF, TF.name, CA.delta, pval) %>% dplyr::rename("dTAG.delta" = "CA.delta", "dTAG.p" = "pval"),
  chem.tested %>% filter(tot.read.count_BRM >= 100 & tot.read.count_DMSO >= 100) %>% dplyr::select(TF, TF.name, CA.delta, pval) %>% dplyr::rename("chem.delta" = "CA.delta", "chem.p" = "pval"),
  by = join_by(TF, TF.name)
) %>% 
  na.omit() %>%
  mutate(class = case_when(
    dTAG.p < .05 & chem.p < .05 ~ "both",
    dTAG.p < .05 & chem.p >= .05 ~ "dTAG",
    dTAG.p >= .05 & chem.p < .05 ~ "chem",
    dTAG.p >= .05 & chem.p >= .05 ~ "none"
  )) %>%
  ggplot(aes(TF, fill = class)) +
  geom_bar(position = "dodge") +
  xlab(NULL) +
  theme_bw() +
  theme(text = element_text(size = 18))


full_join(
  dtag.tested %>% dplyr::select(TF, TF.name, CA.delta, pval) %>% dplyr::rename("dTAG.delta" = "CA.delta", "dTAG.p" = "pval"),
  chem.tested %>% dplyr::select(TF, TF.name, CA.delta, pval) %>% dplyr::rename("chem.delta" = "CA.delta", "chem.p" = "pval"),
  by = join_by(TF, TF.name)
) -> target.motifs

TFBSs[target.motifs$TF.name]

library(DESeq2)
dds = readRDS("/g/krebs/barzaghi/analyses/nf-core_runs/180923_SWISNF_inhibition_ATACseq/results/bwa/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.rds")
intervals = import.bed("/g/krebs/barzaghi/analyses/nf-core_runs/180923_SWISNF_inhibition_ATACseq/results/bwa/merged_replicate/macs2/broad_peak/consensus/consensus_peaks.mRp.clN.bed")
intervals.chr = GRanges(seqnames = paste0("chr", seqnames(intervals)), ranges(intervals))
intervals.chr$name = intervals$name

dds_60min = DESeqDataSetFromMatrix(
  countData = dds@assays@data$counts[,grep("DMSO_60min_wt|BRM014-10uM_60min_wt", colnames(dds@assays@data$counts))],
  colData = dds@colData[grep("DMSO_60min_wt|BRM014-10uM_60min_wt", colnames(dds@assays@data$counts)),],
  design = ~ Group1
)
dds_60min <- DESeq(dds_60min)
resultsNames(dds_60min)
res_60min <- results(dds_60min, name="Group1_DMSO_vs_BRM014.10uM")

dds_24h = DESeqDataSetFromMatrix(
  countData = dds@assays@data$counts[,grep("DMSO_24h_wt|BRM014-10uM_24h_wt", colnames(dds@assays@data$counts))],
  colData = dds@colData[grep("DMSO_24h_wt|BRM014-10uM_24h_wt", colnames(dds@assays@data$counts)),],
  design = ~ Group1
)
dds_24h <- DESeq(dds_24h)
resultsNames(dds_24h)
res_24h <- results(dds_24h, name="Group1_DMSO_vs_BRM014.10uM")

overlaps = findOverlaps(TFBSs[target.motifs$TF.name], intervals.chr)
target.motifs$interval = NA
target.motifs$interval[queryHits(overlaps)] = intervals.chr[subjectHits(overlaps)]$name

target.motifs %>%
  filter(chem.p < .05 & dTAG.p < .05) %>%
  left_join(., res_24h %>% data.frame() %>% rownames_to_column("interval") %>% dplyr::select(interval, log2FoldChange) %>% dplyr::rename("l2fc_24h" = "log2FoldChange"), by = "interval") %>%
  left_join(., res_60min %>% data.frame() %>% rownames_to_column("interval") %>% dplyr::select(interval, log2FoldChange) %>% dplyr::rename("l2fc_60min" = "log2FoldChange"), by = "interval") %>%
  gather(time, l2fc, l2fc_24h, l2fc_60min) %>%
  ggplot(aes(TF, l2fc, fill = time)) +
  geom_violin(draw_quantiles = c(.25,.5,.75)) +
  viridis::scale_color_viridis() +
  xlab(NULL) +
  theme_bw() +
  theme(text = element_text(size = 18))

library(Rsubread)
samplesheet = "/g/krebs/barzaghi/analyses/nf-core_runs/180923_SWISNF_inhibition_ATACseq/results/bwa/merged_library/samplesheet.txt"
TFBSs[filter(target.motifs, chem.p < .05 & dTAG.p < .05)$TF.name] %>% 
  IRanges::resize(., 200, "center") %>%
  data.frame() %>%
  dplyr::select(absolute.idx, seqnames, start, end, strand) -> annotation.query
colnames(annotation.query) = c("GeneID", "Chr",	"Start",	"End",	"Strand")
featureCounts(
  files = grep("DMSO_24h_wt|BRM014-10uM_24h_wt|DMSO_60min_wt|BRM014-10uM_60min_wt", readr::read_delim(samplesheet)$FileName, value = TRUE), 
  annot.ext = annotation.query,
  isPairedEnd = TRUE, requireBothEndsMapped = TRUE
) -> fC

fC$counts %>%
  data.frame() %>%
  rownames_to_column("TF.name") %>%
  gather(sample, count, -TF.name) %>%
  filter(str_detect(sample, "REP3|REP4", negate = TRUE)) %>% 
  mutate(sample = gsub("_wt_REP.*.mLb.clN.sorted.bam", "", sample)) %>%
  separate("sample", c("sample", "time"), sep = "_") %>%
  group_by(TF.name, sample, time) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  spread(sample, count) %>%
  filter((BRM014.10uM + DMSO) > 10) %>%
  mutate(l2fc = log2((BRM014.10uM+1)/(DMSO+1))) %>%
  dplyr::select(-BRM014.10uM, -DMSO) %>%
  right_join(., filter(target.motifs, chem.p < .05 & dTAG.p < .05), by = "TF.name") %>%
  filter(!is.na(time)) %>%
  mutate(time = factor(time, levels = c("60min", "24h"))) %>%
  ggplot(aes(time, l2fc)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(method = "wilcox", label = "p.signif", comparisons = list(c("60min", "24h"))) +
  facet_grid(~TF) +
  theme_bw() +
  theme(text = element_text(size = 18))

fC$counts %>%
  data.frame() %>%
  rownames_to_column("TF.name") %>%
  gather(sample, count, -TF.name) %>%
  filter(str_detect(sample, "REP3|REP4", negate = TRUE)) %>% 
  mutate(sample = gsub("_wt_REP.*.mLb.clN.sorted.bam", "", sample)) %>%
  separate("sample", c("sample", "time"), sep = "_") %>%
  group_by(TF.name, sample, time) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  spread(sample, count) %>%
  filter((BRM014.10uM + DMSO) > 10) %>%
  mutate(l2fc = log2((BRM014.10uM+1)/(DMSO+1))) %>%
  dplyr::select(-BRM014.10uM, -DMSO) %>%
  spread(time, l2fc) %>%
  right_join(., filter(target.motifs, chem.p < .05 & dTAG.p < .05), by = "TF.name") %>%
  filter(!is.na(`24h`)) %>%
  ggplot(aes(`60min`, `24h`)) +
  ggpointdensity::geom_pointdensity() +
  geom_abline(color = "grey", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "grey", linetype = 2) +
  facet_grid(~TF) +
  viridis::scale_color_viridis() +
  theme_bw() +
  theme(text = element_text(size = 18))


