library(tidyverse)
library(magrittr)
library(DESeq2)
library(BSgenome.Mmusculus.UCSC.mm10)
source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source("/g/krebs/barzaghi/Rscripts/IGV_plotting.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
detach("package:plyranges")

PROseq = qs::qread("/g/krebs/barzaghi/analyses/10.05.23_F1_PROseq/2024-04-04_PROseq.qs") # "/g/krebs/barzaghi/analyses/10.05.23_F1_PROseq/PROseq_processing.R"
f1_lofs_smf = qs::qread("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/f1_lofs_smf_2024-05-09.qs")
GenomicTiles = Load.GenomicTiles(tiles.width = 80)

PROseq %<>%
  filter(TF.name %in% f1_lofs_smf$TF.name) %>%
  unnest(cols = c(rep, A, R, Run)) %>%
  group_by(TF.name, Sample, log2FoldChange, padj) %>%
  summarise(R = sum(R, na.rm=TRUE), A=sum(A, na.rm=TRUE), .groups = "drop") %>%
  filter(R > 10) %>%
  mutate(PROseq.change = ifelse(padj<.05 & log2FoldChange < 0, "down", ifelse(padj<.05 & log2FoldChange > 0, "up", "no")))

f1_lofs_smf %>%
  left_join(., dplyr::rename(PROseq, PROseq.r = R, PROseq.a = A, PROseq.padj = padj), by = c("Sample", "TF.name")) %>%
  filter(!is.na(PROseq.change)) -> f1_smf_proseq

table(f1_smf_proseq$change)[1]
table(f1_smf_proseq$change, f1_smf_proseq$PROseq.change)[1,1] / 
  (table(f1_smf_proseq$change, f1_smf_proseq$PROseq.change)[1,1] + table(f1_smf_proseq$change, f1_smf_proseq$PROseq.change)[1,3]) *
  100 # 90%

# A
f1_smf_proseq %>%
  filter(change == "down") %>%
  ggplot() +
  geom_point(aes(log2FoldChange, -log10(PROseq.padj), color = PROseq.change)) +
  scale_color_manual(values = c("red", "black", "blue"), breaks = c("down", "no", "up")) +
  scale_y_continuous(breaks = c(0,15)) +
  scale_x_continuous(breaks = c(-6,0,3)) +
  xlab("qPRO-seq log2(FC)") +
  ylab("-log10(p-adj)") +
  theme_bw() +
  theme(text = element_text(size = 18),  legend.position = c(.8,.88), legend.background = element_blank(), legend.title = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig6a_", Sys.Date(), ".pdf"), pl, width = 3.5, height = 3.5)

# B
f1_smf_proseq %>%
  mutate(accessibility.delta.bins = factor(ifelse(p.adj>.05, "n.s.", ifelse(CA.delta >= -30, "0-30%", "31-100%")), levels = c("n.s.", "0-30%", "31-100%"))) %>%
  group_by(accessibility.delta.bins) %>% filter(all(c("enhancer", "promoter") %in% locus)) %>% ungroup() -> pl.df

pl.df %>%
  filter(locus != "Ctcf") %>%
  group_by(locus, accessibility.delta.bins) %>%
  filter(n() > 10) %>%
  group_by(locus) %>%
  filter(all(c("n.s.", "0-30%", "31-100%") %in% accessibility.delta.bins)) %>%
  ungroup() %>%
  ggplot(aes(accessibility.delta.bins, log2FoldChange)) +
  geom_boxplot(aes(fill = locus)) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", size = 7, ref.group = "n.s.", vjust = 1) +
  facet_grid(~locus) +
  xlab("CA frequency delta (%)") + ylab("transcriptional change (l2fc)") +
  scale_fill_manual(values = viridis::mako(n=4)[2:3], breaks = c("enhancer", "promoter")) +
  scale_y_continuous(breaks = c(-6,0,4)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig6b_", Sys.Date(), ".pdf"), pl, width = 5, height = 4)
pl.df %>%
  filter(locus != "Ctcf") %>%
  group_by(locus, accessibility.delta.bins) %>%
  summarise(n = n(), .groups = "drop") %>% # detach("package:plyranges")
  spread(accessibility.delta.bins, n)

# C
partition.collapsing.dictionary = split(1:12,1:12)[c(4,12,7,10,9,3,2,8,1,5,11,6)]
pl.df %>%
  filter(TF.name == "TFBS_2708927", Sample == "STKO") %>%
  mutate(TFBS.cluster = "GenomicTile_16920562") %>%
  patch.single.site.plots(
    interpretable.master.table = ., rank = 1, k = 12,
    pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictionary, 
    data.type = "F1_bait.capture", remove.TFBS.labels = TRUE, patch.version = 2
  ) -> pl
pl$pl
max(pl$chromatin.influence.df[2,]$acc.width.distro[[1]]) # 297bp
pl$chromatin.influence.df[2,] %>%
  mutate(CA.R = (accessible.R/tot.R)*100, CA.A = (accessible.A/tot.A)*100) %>%
  unnest(c(acc.width.distro, acc.read.count.distro.R, acc.read.count.distro.A)) %>%
  filter(acc.width.distro >= 100) %>%
  group_by(TF.name, TF, tot.R, tot.A, CA.R, CA.A) %>%
  summarise(CA.wide.R = sum(acc.read.count.distro.R), CA.wide.A = sum(acc.read.count.distro.A), .groups = "drop") %>%
  mutate(CA.wide.R = CA.wide.R/tot.R*100, CA.wide.A = CA.wide.A/tot.A*100) %>%
  mutate(CA.delta = CA.wide.A - CA.wide.R) # R: 66%, A: 40%
filter(pl.df, TF.name == "TFBS_2708927", Sample == "STKO")$log2FoldChange # -0.9711348 (PROseq)
filter(pl.df, TF.name == "TFBS_2708927", Sample == "STKO")$PROseq.padj # 0.002205783 (PROseq)
filter(pl.df, TF.name == "TFBS_2708927", Sample == "STKO")$p.adj # 5.942256e-05 (SMF)
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig6c_", Sys.Date(), ".png"), width = 30, height = 30, units = "cm", res = 300)
pl$pl
dev.off()

RegionOfInterest = IRanges::resize(TFBSs["TFBS_2708927"], width = 9000, fix = "center")
plot_genomic_track(
  sampleSheet = "/g/krebs/barzaghi/analyses/nf-core_runs/130923_BL6_Spret_F1_ATAC/aln_merged_deduplicated/samplesheet_PooledReplicates.txt",
  samples = "BL6_Spret_F1",
  # sampleSheet = "/g/krebs/barzaghi/HTS/ATAC-seq/Heard_NatGen_2017/aln/Qinput_male.txt",
  # samples = "male_R1",
  RegionOfInterest = RegionOfInterest,
  allelic = TRUE,
  tile.width = 1000,
  tile.step = 100,
  max.y.lim = 500,
  color = c("black", "grey45"), 
  y.labs = c("Bl6 ATAC", "Spret ATAC"),
  delta = FALSE,
  normalise = FALSE
) -> ATAC_track
plot_genomic_track(
  sampleSheet = ChIP_data_dictionary$Klf4,
  samples = "Klf4",
  RegionOfInterest = RegionOfInterest,
  allelic = FALSE,
  tile.width = 200,
  tile.step = 10,
  max.y.lim = 800,
  color = "black", 
  y.labs = "Klf4"
) -> Klf4_track
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
plot_genomic_track(
  sampleSheet = "/g/krebs/barzaghi/HTS/PROseq/Qinput_all_final.txt",
  samples = "Spret",
  RegionOfInterest = RegionOfInterest,
  allelic = TRUE, 
  tile.width = 4,
  tile.step = 2,
  max.y.lim = 50,
  color = c("black", "grey45"), 
  y.labs = c("Bl6 PROseq", "Spret PROseq"),
  delta = FALSE,
  normalise = FALSE, 
  plot.coordinates = TRUE
) -> PROseq_track
# Not aligned with QuasR
# plot_genomic_track( 
#   sampleSheet = "/g/krebs/barzaghi/HTS/RNAseq/OutDir/SPRET_TKO/Qinput.txt",
#   samples = c("BL6", "Spretus"),
#   RegionOfInterest = RegionOfInterest,
#   allelic = TRUE, 
#   tile.width = 4,
#   tile.step = 2,
#   max.y.lim = 50,
#   color = c("black", "grey45"), 
#   y.labs = c("Bl6 PROseq", "Spret PROseq"),
#   delta = FALSE,
#   normalise = FALSE, 
#   plot.coordinates = TRUE
# ) -> RNAseq_track
p_final <- ATAC_track + Klf4_track + MNase_track + PROseq_track +
  plot_layout(ncol = 1, heights = c(1/3, 1/6, 1/6, 1/3))
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig6b_tracks_", Sys.Date(), ".pdf"), width = 9, height = 5)
p_final
dev.off()

plot_genomic_track(
  sampleSheet = "/g/krebs/barzaghi/HTS/PROseq/Qinput_all_final.txt",
  samples = "Spret",
  RegionOfInterest = IRanges::resize(RegionOfInterest, 500, "center"),
  allelic = TRUE, 
  tile.width = 4,
  tile.step = 2,
  max.y.lim = 50,
  color = c("black", "grey45"), 
  y.labs = c("Bl6 PROseq", "Spret PROseq"),
  delta = FALSE,
  normalise = FALSE, 
  plot.coordinates = TRUE
) -> PROseq_track
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig6b_tracks_inlet_", Sys.Date(), ".pdf"), width = 3, height = 2)
PROseq_track
dev.off()
