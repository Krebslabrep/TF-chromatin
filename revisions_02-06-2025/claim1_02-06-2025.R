library(tidyverse)
library(magrittr)
library(QuasR)
library(ComplexHeatmap)
library(SingleMoleculeFootprinting)
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/utils.r")
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/source_FootprintCharter_functions.r")
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/IGV_plotting.r")
detach("package:plyranges")

setwd("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/")

chromHMM = Load.chromHMM(GenomicTiles = TRUE)
chromHMM_gr = Load.chromHMM(GenomicTiles = FALSE)
GenomicTiles = Load.GenomicTiles(tiles.width = 80)
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

###################
##### Heatmap #####
###################
regions = Load.chromHMM(GenomicTiles = FALSE) %>%
  filter(chromHMM %in% c("intergenic", "enhancer", "promoter", "Ctcf")) %>%
  IRanges::resize(., 5000, "center")
regions.list = split(regions, regions$chromHMM)
lapply(seq_along(regions.list), function(i){
  plyranges::write_bed(regions.list[[i]], file = paste0("revisions_02-06-2025/", names(regions.list)[i], ".bed"))
})
# /g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/deepTools_computeMatrix.sh
# /g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/deepTools_plotHeatmap.sh

###################
#### Sox2 LCR #####
###################
# PMID: 25512558 | ~104-111kb downstream of Sox2 chr3:34647066-34655389
Sox2_LCR_tracks = function(RegionOfInterest){
  
  plot_genomic_track(
    sampleSheet = "/g/krebs/barzaghi/HTS/DHS/Qinput.txt",
    samples = "DNAse-seq_mESC_TKO",
    RegionOfInterest = RegionOfInterest,
    max.y.lim = 800,
    y.labs = "DNase"
  ) + 
    geom_segment(y = 700, yend = 700, x = start(RegionOfInterest)+200, xend = start(RegionOfInterest)+1200) +
    geom_text(data = data.frame(x=NA), aes(x = start(RegionOfInterest)+1400, y = 700, label = "1 kb"), inherit.aes = FALSE) +
    ggtitle("Sox2 LCR") +
    theme(plot.title = element_text(hjust = 0.5)) -> DHS_track
  
  plot_genomic_track(
    sampleSheet = ChIP_data_dictionary$H3K27Ac,
    samples = "H3K27Ac_mESC_DNMT_TKO",
    RegionOfInterest = RegionOfInterest,
    max.y.lim = 500,
    y.labs = "H3K27Ac"
  ) -> H3K27Ac_track
  
  plot_genomic_track(
    sampleSheet = ChIP_data_dictionary$H3K4me3,
    samples = "H3K4me3_mESC_NT",
    RegionOfInterest = RegionOfInterest,
    max.y.lim = 150,
    y.labs = "H3K4me3"
  ) -> H3K4me3_track
  
  plot_genomic_track(
    sampleSheet = ChIP_data_dictionary$H3K4me1,
    samples = "H3K4me1_mESC_NT",
    RegionOfInterest = RegionOfInterest,
    max.y.lim = 60,
    y.labs = "H3K4me1",
  ) -> H3K4me1_track
  
  chromHMM_gr %>%
    plyranges::filter_by_overlaps(., RegionOfInterest) %>%
    plyranges::filter(chromHMM == "enhancer") %>% 
    data.frame() %>%
    dplyr::select(start, end) %>%
    ggplot() +
    geom_rect(
      aes(xmin = start, xmax = end, ymin = 0, ymax = 1), inherit.aes = FALSE
    ) +
    scale_x_continuous(labels = NULL) +
    scale_y_continuous(limits = c(0,1), breaks = c(0,1), labels = c("", ""), expand = expansion(mult = 0, add = 0)) +
    theme_bw() +
    theme(
      strip.background = element_blank(), strip.text = element_blank(), axis.title = element_blank(),
      legend.position = "none",
      axis.line.y = element_line(linewidth = 1),
      axis.ticks.length = unit(0, "mm"),
      panel.grid = element_blank(), panel.border = element_blank()
    ) -> enhancers_track
  
  # GenomicTiles[CA_loci_df$TFBS.cluster]["GenomicTile_10308748"] %>%
  GenomicTiles[filter(CA_loci_df, chromHMM_annotation == "enhancer")$TFBS.cluster] %>%
    plyranges::filter_by_overlaps(., RegionOfInterest) %>%
    data.frame() %>%
    dplyr::select(start, end) %>%
    ggplot() +
    geom_rect(
      aes(xmin = start, xmax = end, ymin = 0, ymax = 1), inherit.aes = FALSE
    ) +
    scale_x_continuous(
      limits = c(start(RegionOfInterest), end(RegionOfInterest)), 
      breaks = c(start(RegionOfInterest), end(RegionOfInterest)), 
      labels = format(c(start(RegionOfInterest), end(RegionOfInterest)), nsmall=1, big.mark=","), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,1), breaks = c(0,1), labels = c("", ""), expand = expansion(mult = 0, add = 0)) +
    theme_bw() +
    theme(
      strip.background = element_blank(), strip.text = element_blank(), axis.title = element_blank(),
      legend.position = "none",
      axis.line.y = element_line(linewidth = 1),
      axis.ticks.length = unit(0, "mm"),
      panel.grid = element_blank(), panel.border = element_blank()
    ) -> SMF_covered_enhancers_track
  
  p_final <- DHS_track + H3K27Ac_track + H3K4me1_track + H3K4me3_track + enhancers_track + SMF_covered_enhancers_track +
    plot_layout(ncol = 1, heights = c(rep(96/4/100, 4), rep(2/100, 2)))
  
  return(p_final)
  
}

# Sox2 gene
RegionOfInterest_Sox2_gene = IRanges::resize(GRanges("chr3", IRanges(34647066, 34655389)), 20000, "center")
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/SOX2gene_tracks_", Sys.Date(), ".pdf"), width = 4, height = 5)
Sox2_LCR_tracks(RegionOfInterest = RegionOfInterest_Sox2_gene)
dev.off()

library(plotgardener)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/SOX2gene_", Sys.Date(), ".png"), width = 10, height = 5, units = "cm", res = 300)
plotGenes(
  assembly = "mm10",
  chrom = "chr3", chromstart = 34650405, chromend = 34652461, strandLabels = FALSE, fontsize = 0
)
dev.off()

# LCR
RegionOfInterest_LCR = IRanges::resize(GRanges("chr3", IRanges(34647066+104000, 34655389+111000)), 50000, "center")
pdf(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/SOX2LCR_tracks_", Sys.Date(), ".pdf"), width = 10, height = 5)
Sox2_LCR_tracks(RegionOfInterest = RegionOfInterest_LCR)
dev.off()

# Single loci
enhancer = data.frame(Sample = "SMF_MM_TKO_DE_", TFBS.cluster = "GenomicTile_10308748")
partition.collapsing.dictinary = split(1:12,1:12)[c(5,12,9,6,3,1,2,7,4,10,11,8)]
patch.single.site.plots(
  interpretable.master.table = enhancer, rank = 1, k = 16,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE, deduplicate = FALSE
) -> pl
process.CA.df(x = mutate(pl$chromatin.influence.df, ChIP = NA), cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) #58% # SMAD2/3/4, ESRRB, SOX2, OCT4
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/SOX2LCR_GenomicTile_10308748_", Sys.Date(), ".png"), width = 32, height = 25, units = "cm", res = 300)
pl$pl
dev.off()

enhancer = data.frame(Sample = "SMF_MM_TKO_DE_", TFBS.cluster = "GenomicTile_10308797")
partition.collapsing.dictinary = split(1:9,1:9)[c(8,9,7,4,1,2,6,5,3)]
patch.single.site.plots(
  interpretable.master.table = enhancer, rank = 1, k = 16,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_bait.capture", remove.TFBS.labels = TRUE, deduplicate = FALSE
) -> pl
process.CA.df(x = mutate(pl$chromatin.influence.df, ChIP = NA), cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) # 27% #SOX2, OCT4, SOX2, KLF4
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/SOX2LCR_GenomicTile_10308797_", Sys.Date(), ".png"), width = 32, height = 25, units = "cm", res = 300)
pl$pl
dev.off()

CA_loci_df %>%
  filter(chromHMM_annotation == "enhancer") %>%
  mutate(LCR = ifelse(TFBS.cluster %in% names(plyranges::filter_by_overlaps(GenomicTiles, RegionOfInterest_LCR)), "Sox2 LCR", "other enhancers")) %>%
  ggplot(aes(LCR, CA_regulatory), fill = "transparent") +
  geom_boxplot() +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox", label.x = 1.5) +
  xlab(NULL) +
  ylab("CA frequency %") +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/enhancers_boxplot_", Sys.Date(), ".pdf"), pl, width = 3, height = 4.5)

###################
####### CTCF ######
###################
inactive.median = as.integer(quantile(filter(CA_loci_df, chromHMM_annotation == "intergenic", ChIP_annotation == "unbound")$CA_regulatory, 0.50))
CA_loci_df_f1 %>%
  filter(motif.change == "l.o.f.", ChIP_annotation == "bound") %>%
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
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/CTCF_", Sys.Date(), ".pdf"), pl, width = 4.5, height = 4.5)
