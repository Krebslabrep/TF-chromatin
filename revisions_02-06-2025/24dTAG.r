library(tidyverse)
library(magrittr)
library(dplyr)
library(QuasR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SingleMoleculeFootprinting)
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/utils.r")
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/source_FootprintCharter_functions.r")
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/IGV_plotting.r")
detach("package:plyranges")

chromHMM = Load.chromHMM(GenomicTiles = TRUE)
GenomicTiles = Load.GenomicTiles(tiles.width = 80, GeneAssociation = FALSE, CGI = FALSE)
TFBSs = Load.TFBSs()

CA_loci_df_dTAG = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2025-02-24_chromatin.influence_Sox2_dTAG_run2_run3_PooledReps.df.qs") %>%
  process.CA.df(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>%
  filter(
    (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
      (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
      (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
  ) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF))

sox2.motifs = TFBSs[unique(filter(CA_loci_df_dTAG, TF == "Sox2", ChIP_annotation == "bound")$TF.name)] %>%
  IRanges::resize(., 500, "center")

# data.frame(
#   FileName = list.files("/g/krebs/barzaghi/analyses/nf-core_runs/190923_deWit_EMBO_2023_ATACseq/results/bwa/merged_replicate", ".bam$", full.names = TRUE)
#   ) %>%
#   filter(str_detect(FileName, "SOX2")) %>%
#   filter(str_detect(FileName, "NT|2H|24H")) %>%
#   mutate(SampleName = c("24H", "2H", "NT")) %>%
#   readr::write_tsv(., "deWit_sox2_dTAG_samplesheet.txt")

seqlevels(sox2.motifs) = gsub("^_", "", gsub("chr", "", seqlevels(sox2.motifs)))
seqnames(sox2.motifs) = gsub("^_", "", gsub("chr", "", seqnames(sox2.motifs)))
wg = GRanges(BSgenome.Mmusculus.UCSC.mm10@seqinfo)[1:19]
seqlevels(wg) = gsub("^_", "", gsub("chr", "", seqlevels(wg)))
seqnames(wg) = gsub("^_", "", gsub("chr", "", seqnames(wg)))

prj = qAlign("deWit_sox2_dTAG_samplesheet.txt", genome = "BSgenome.Mmusculus.UCSC.mm10", 
             projectName = "prj", paired = "fr", aligner = "Rbowtie", 
             bisulfite = "no")

counts = qCount(proj = prj, query = sox2.motifs)
counts %<>%data.frame() %>% dplyr::select(-width)
counts.wg = alignmentStats(prj)[,"mapped"]
counts.norm = counts/counts.wg*1e6

# Whole-genome
counts.norm %>%
  arrange(NT) %>%
  mutate(rank = seq_along(NT)) %>%
  gather(treatment, cpm, -rank) %>%
  mutate(treatment = case_when(
    treatment == "NT" ~ "untreated",
    treatment == "X2H" ~ "2h dTAG",
    treatment == "X24H" ~ "24h dTAG"
  )) -> pl.df
pl.df %>%  ggplot(aes(rank, cpm, color = treatment)) +
  geom_point(size = .25, alpha = .1) +
  geom_smooth(method = "gam") +
  xlab("motifs rank by untreated cpm") + ylab("cpm") +
  scale_color_manual(values = c("grey", "salmon", "darkred"), breaks = c("untreated", "2h dTAG", "24h dTAG")) +
  scale_x_continuous(breaks = c(0,max(pl.df$rank)), limits = c(0,max(pl.df$rank))) +
  # scale_y_continuous(breaks = c(0,inactive.median,50,100), limits = c(0,100)) + 
  theme_bw() +
  theme(text = element_text(size = 18), legend.position=c(.25, .9), legend.title = element_blank(), legend.background = element_blank()) -> pl
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/24h_sox2_dTAG_", Sys.Date(), ".pdf"), width = 4.5, height = 4.5)
pl
dev.off()

# 
counts.norm %>%
  mutate(l2FC_2h = log2((X2H+1)/(NT+1)), l2FC_24h = log2((X24H+1)/(NT+1))) %>%
  gather(treatment_time, l2FC, l2FC_2h, l2FC_24h) %>%
  mutate(treatment_time = gsub("l2FC_", "", treatment_time)) %>%
  ggplot(aes(l2FC, fill = treatment_time)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c("salmon", "darkred"), breaks = c("2h", "24h")) +
  theme_bw()

counts.norm %>%
  gather(treatment, cpm) %>%
  mutate(treatment = case_when(
    treatment == "NT" ~ "untreated",
    treatment == "X2H" ~ "2h dTAG",
    treatment == "X24H" ~ "24h dTAG"
  )) %>%
  ggplot(aes(cpm, fill = treatment)) +
  geom_density(alpha = .5) +
  scale_color_manual(values = c("grey", "salmon", "darkred"), breaks = c("untreated", "2h dTAG", "24h dTAG")) +
  theme_bw()

counts.norm %>%
  ggplot(aes(X2H, X24H)) +
  geom_abline(color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "grey", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  ggpointdensity::geom_pointdensity() +
  viridis::scale_color_viridis() +
  xlab("accessibility at 2h") + ylab("accessibility at 24h") + 
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl

counts.norm %>%
  mutate(l2FC_2h = log2((X2H+1)/(NT+1)), l2FC_24h = log2((X24H+1)/(NT+1))) %>%
  ggplot() +
  geom_abline(color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "grey", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  ggpointdensity::geom_pointdensity(aes(l2FC_2h, l2FC_24h)) +
  viridis::scale_color_viridis() +
  xlab("l2FC (2h - untreated)") + ylab("l2FC (24h - untreated)") + 
  theme_bw() +
  theme(text = element_text(size = 18), legend.position=c(.25, .75), legend.background = element_blank()) -> pl
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/24h_sox2_dTAG_scatter_", Sys.Date(), ".pdf"), width = 4.5, height = 4.5)
pl
dev.off()

# Single locus
# counts.norm %>%
#   filter(X24H > 0, X2H > 0, NT > 3) %>%
#   arrange(abs(X24H-X2H)) %>%
#   head(20)
RegionOfInterest = TFBSs["TFBS_9518804"] %>% IRanges::resize(., 5000, "center")
seqlevels(RegionOfInterest) = gsub("chr", "", seqlevels(RegionOfInterest))
seqnames(RegionOfInterest) = gsub("chr", "", seqnames(RegionOfInterest))

plot_genomic_track(
  sampleSheet = "/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/deWit_sox2_dTAG_samplesheet.txt",
  samples = "NT",
  RegionOfInterest = RegionOfInterest,
  normalize = TRUE, 
  normalization.factor = as.integer(counts.wg["NT:genome"]),
  max.y.lim = 4,
  y.labs = "NT"
) + 
  geom_segment(y = 2.25, yend = 2.25, x = start(RegionOfInterest)+200, xend = start(RegionOfInterest)+1200) +
  geom_text(data = data.frame(x=NA), aes(x = start(RegionOfInterest)+1400, y = 2.25, label = "1 kb"), inherit.aes = FALSE) +
  # ggtitle("enhancer") +
  theme(plot.title = element_text(hjust = 0.5)) -> NT_track

plot_genomic_track(
  sampleSheet = "/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/deWit_sox2_dTAG_samplesheet.txt",
  samples = "2H",
  RegionOfInterest = RegionOfInterest,
  normalize = TRUE, 
  normalization.factor = as.integer(counts.wg["2H:genome"]),
  max.y.lim = 4,
  y.labs = "2H"
) -> `2H_track`

plot_genomic_track(
  sampleSheet = "/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/deWit_sox2_dTAG_samplesheet.txt",
  samples = "24H",
  RegionOfInterest = RegionOfInterest,
  normalize = TRUE, 
  normalization.factor = as.integer(counts.wg["24H:genome"]),
  max.y.lim = 4,
  y.labs = "24H",
  plot.coordinates = TRUE
) -> `24H_track`

p_final <- NT_track + `2H_track` + `24H_track` + 
  plot_layout(ncol = 1, heights = c(1/3, 1/3, 1/3))

pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/24h_sox2_dTAG_locus_", Sys.Date(), ".pdf"), width = 9, height = 4)
p_final
dev.off()