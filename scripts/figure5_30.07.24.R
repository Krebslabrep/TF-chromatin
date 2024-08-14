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

qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-06-04_chromatin.influence_F1_PooledReplicates.df.qs") %>%
  process.CA.df_f1(x = ., cre.annotation = chromHMM, chip.thr = ChIP_thresholds_dictionary_lenient) %>% filter(
  (TF == "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation == "repressed") |
    (TF != "Rest" & ChIP_annotation %in% c("bound") & chromHMM_annotation %in% c("Ctcf", "enhancer", "bivalent_promoter", "promoter", "transcription")) |
    (ChIP_annotation %in% c("unbound") & chromHMM_annotation %in% c("intergenic", "repressed", "heterochromatin"))
) %>%
  mutate(TF = ifelse(ChIP_annotation == "unbound", "unbound", TF)) -> CA_loci_df_f1

# A
TFBSs[unique(filter(CA_loci_df_f1, ChIP_annotation == "bound")$TF.name)] %>%
  SingleMoleculeFootprinting::Arrange_TFBSs_clusters(., add.single.TFs = FALSE, max_cluster_width = 300) -> TFBS.clusters
TFBS.cluster.compositions = unlist(TFBS.clusters$ClusterComposition)
TFBS.cluster.compositions$cluster.id = gsub("\\..*", "", names(TFBS.cluster.compositions))

CA_loci_df_f1 %>% 
  filter(ChIP_annotation == "bound") %>%
  full_join(., TFBS.cluster.compositions %>%
    data.frame() %>%
    dplyr::select(absolute.idx, cluster.id) %>%
    dplyr::rename("TF.name" = "absolute.idx"),
    by = "TF.name"
) %>%
  mutate(cluster.id = ifelse(is.na(cluster.id), TF.name, cluster.id)) %>%
  na.omit() %>%
  group_by(cluster.id) %>% mutate(nr.motifs = n()) %>% ungroup() %>%
  filter(motif.change == "l.o.f.") %>%
  group_by(nr.motifs) %>% filter(n() > 30) %>% ungroup() %>%
  dplyr::select(Sample, TF, TF.name, nr.motifs, CA_regulatory_R, CA_regulatory_A, width_regulatory_R, width_regulatory_A) %>%
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
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75), alpha = .5, linetype = 1) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75), alpha = .5, linetype = 1) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = nr.motifs), max.overlaps = 100, show.legend = FALSE) +
  scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
  scale_x_continuous(breaks = c(100,200,300), limits = c(100,300)) +
  xlab("width (bp)") + ylab("frequency (%)") +
  scale_color_manual(values = c("black","salmon"), breaks = c("ref", "alt")) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = c(.8,.87), legend.title = element_blank(), legend.background = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig5a_", Sys.Date(), ".pdf"), pl, width = 3, height = 3.3)

# B
CA_loci_df_f1 %>% 
  filter(ChIP_annotation == "bound") %>%
  full_join(., TFBS.cluster.compositions %>%
              data.frame() %>%
              dplyr::select(absolute.idx, cluster.id) %>%
              dplyr::rename("TF.name" = "absolute.idx"),
            by = "TF.name"
  ) %>%
  mutate(cluster.id = ifelse(is.na(cluster.id), TF.name, cluster.id)) %>%
  na.omit() %>%
  group_by(cluster.id) %>% mutate(nr.motifs = n()) %>% ungroup() %>%
  filter(motif.change == "l.o.f.") %>%
  group_by(nr.motifs) %>% filter(n() > 30) %>% ungroup() %>%
  dplyr::select(Sample, TFBS.cluster, TF, TF.name, nr.motifs, CA_regulatory_R, CA_regulatory_A, width_regulatory_R, width_regulatory_A) %>%
  pivot_longer(cols = c(CA_regulatory_R, CA_regulatory_A, width_regulatory_R, width_regulatory_A)) %>%
  mutate(name = gsub("_regulatory", "", gsub("CA", "frequency", name))) %>%
  separate("name", into = c("name", "allele"), sep = "_") %>%
  spread(name, value) %>%
  mutate(allele = factor(ifelse(allele == "R", "ref", "alt"), levels = c("ref", "alt"))) %>%
  pivot_wider(id_cols = c("Sample", "TFBS.cluster", "TF", "TF.name", "nr.motifs"), names_from = "allele", values_from = c("frequency", "width")) %>%
  na.omit() %>%
  filter(((frequency_alt - frequency_ref) < -5) & ((width_alt - width_ref) < -10) & nr.motifs > 1) -> pl.df


lapply(seq(nrow(pl.df)), function(i){
  
  patch.single.site.plots(
    interpretable.master.table = pl.df, rank = i, k = 12,
    pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = NULL, 
    data.type = "F1_bait.capture", remove.TFBS.labels = FALSE, plotting.TFBSs = NULL, 
    reference.genome = BSgenome.Mmusculus.UCSC.mm10, deduplicate = FALSE, reutrn.chromatin.influence.df = FALSE
  ) -> pl
  png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/2024-08-01_SingleSites/", i, ".png"), width = 30, height = 30, units = "cm", res = 300)
  print(pl$pl)
  dev.off()
  
})

"GenomicTile_52339787"
"GenomicTile_16569692"

partition.collapsing.dictinary = NULL
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_5335377", Sample = "STKO"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "F1_bait.capture", remove.TFBS.labels = FALSE, plotting.TFBSs = NULL, 
  reference.genome = BSgenome.Mmusculus.UCSC.mm10, deduplicate = FALSE, reutrn.chromatin.influence.df = FALSE
) -> pl
pl$pl

roi = IRanges::resize(GenomicTiles["GenomicTile_5335379"], 1000, "center")
CallContextMethylation(
  sampleSheet = "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_all_bait_capture_data_pooledReplicates.txt", 
  sample = c("A_CTKO_dp_rm_DE_", "R_CTKO_dp_rm_DE_"), 
  genome = BSgenome.Mmusculus.UCSC.mm10, RegionOfInterest = roi, 
  returnSM = FALSE
) %>%
  MaskSNPs2(., CytosinesToMaks = CytosinesToMask, MaskSMmat = FALSE, Experiment = "DE") %>%
  PlotAvgSMF(., RegionOfInterest = roi, TFBSs = plyranges::filter_by_overlaps(TFBSs, roi), SNPs = plyranges::filter_by_overlaps(SNPs$Cast, roi))

CallContextMethylation(
  sampleSheet = "/g/krebs/barzaghi/HTS/SMF/MM/2024-05-22-2222T2JNX/Qinput_Oct4_dTAG.txt", 
  sample = c("Oct4NT_NO_", "Oct42h_NO_"), 
  genome = BSgenome.Mmusculus.UCSC.mm10, RegionOfInterest = roi, 
  returnSM = FALSE
)$DGCHN %>%
  PlotAvgSMF(., RegionOfInterest = roi, TFBSs = plyranges::filter_by_overlaps(TFBSs, roi), SNPs = plyranges::filter_by_overlaps(SNPs$Spret, roi))


# B
TFBSs[unique(filter(CA_loci_df_f1, ChIP_annotation == "unbound")$TF.name)] %>%
  SingleMoleculeFootprinting::Arrange_TFBSs_clusters(., add.single.TFs = FALSE, max_cluster_width = 300) -> TFBS.clusters
TFBS.cluster.compositions = unlist(TFBS.clusters$ClusterComposition)
TFBS.cluster.compositions$cluster.id = gsub("\\..*", "", names(TFBS.cluster.compositions))

CA_loci_df_f1 %>%
  inner_join(., TFBS.cluster.compositions %>%
               data.frame() %>%
               dplyr::select(absolute.idx, cluster.id) %>%
               dplyr::rename("TF.name" = "absolute.idx"),
             by = "TF.name"
  ) %>%
  mutate(TF = "unbound") %>%
  group_by(cluster.id, TF) %>% mutate(nr.homotypic.cobinders = sum(TF == TF)) %>% ungroup() %>%
  filter(motif.change == "l.o.f.") %>%
  group_by(TF, nr.homotypic.cobinders) %>% filter(n() > 20) %>% ungroup() %>%
  dplyr::select(Sample, TF, TF.name, nr.homotypic.cobinders, CA_regulatory_R, CA_regulatory_A, width_R, width_A) %>%
  pivot_longer(cols = c(CA_regulatory_R, CA_regulatory_A, width_R, width_A)) %>%
  mutate(name = gsub("CA_regulatory", "frequency", name)) %>%
  separate("name", into = c("name", "allele"), sep = "_") %>%
  spread(name, value) %>%
  mutate(allele = factor(ifelse(allele == "R", "ref", "alt"), levels = c("ref", "alt"))) %>%
  group_by(TF, allele, nr.homotypic.cobinders) %>%
  summarise(
    width_0.25 = quantile(width, 0.25, na.rm=TRUE), 
    width_0.50 = quantile(width, 0.50, na.rm=TRUE), 
    width_0.75 = quantile(width, 0.75, na.rm=TRUE), 
    frequency_0.25 = quantile(frequency, 0.25, na.rm=TRUE), 
    frequency_0.50 = quantile(frequency, 0.50, na.rm=TRUE), 
    frequency_0.75 = quantile(frequency, 0.75, na.rm=TRUE), 
    .groups = "drop"
  ) %>%
  ggplot(aes(width_0.50, frequency_0.50, color = allele)) +
  geom_errorbar(aes(ymin = frequency_0.25, ymax = frequency_0.75), alpha = .5, linetype = 1) +
  geom_errorbarh(aes(xmin = width_0.25, xmax = width_0.75), alpha = .5, linetype = 1) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = nr.homotypic.cobinders), max.overlaps = 100, show.legend = FALSE) +
  scale_y_continuous(breaks = c(0,50,100), limits = c(0,100)) +
  scale_x_continuous(breaks = c(50,100,150,200), limits = c(50,200)) +
  xlab("width (bp)") + ylab("frequency (%)") +
  scale_color_manual(values = c("black","salmon"), breaks = c("ref", "alt")) +
  facet_grid(~TF) +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = c(.8,.87), legend.title = element_blank(), legend.background = element_blank()) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/fig5b_", Sys.Date(), ".pdf"), pl, width = 3, height = 3.3)

# C
TFBSs[unique(filter(CA_loci_df_f1, ChIP_annotation != "unbound")$TF.name)] %>%
  SingleMoleculeFootprinting::Arrange_TFBSs_clusters(., add.single.TFs = FALSE, max_cluster_width = 300) -> TFBS.clusters
TFBS.cluster.compositions = unlist(TFBS.clusters$ClusterComposition)
TFBS.cluster.compositions$cluster.id = gsub("\\..*", "", names(TFBS.cluster.compositions))

CA_loci_df_f1 %>%
  inner_join(., TFBS.cluster.compositions %>%
               data.frame() %>%
               dplyr::select(absolute.idx, cluster.id) %>%
               dplyr::rename("TF.name" = "absolute.idx"),
             by = "TF.name"
  ) %>%
  group_by(cluster.id, TF) %>% mutate(nr.homotypic.cobinders = sum(TF == TF)) %>% ungroup() %>%
  filter(motif.change == "l.o.f." & ChIP_annotation == "bound" & TF == "Klf4") %>%
  group_by(TF, nr.homotypic.cobinders) %>% filter(n() > 20) %>% ungroup() %>%
  dplyr::select(Sample, TF, TF.name, nr.homotypic.cobinders, CA_regulatory_R, CA_regulatory_A, width_R, width_A) %>%
  pivot_longer(cols = c(CA_regulatory_R, CA_regulatory_A, width_R, width_A)) %>%
  mutate(name = gsub("CA_regulatory", "frequency", name)) %>%
  separate("name", into = c("name", "allele"), sep = "_") %>%
  spread(name, value) %>%
  mutate(allele = factor(ifelse(allele == "R", "ref", "alt"), levels = c("ref", "alt"))) %>%
  group_by(TF, allele, nr.homotypic.cobinders) %>%
  mutate(width_0.50 = quantile(width, 0.50, na.rm=TRUE), frequency_0.50 = quantile(frequency, 0.50, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(frequency_residual = frequency - frequency_0.50, width_residual = width - width_0.50) %>%
  dplyr::select(Sample, TF.name, nr.homotypic.cobinders, allele, frequency_residual, width_residual) %>%
  pivot_wider(id_cols = c(Sample, TF.name, nr.homotypic.cobinders), names_from = allele, values_from = c(frequency_residual, width_residual)) %>%
  filter(nr.homotypic.cobinders %in% c(3,4) & abs(frequency_residual_ref) < 5 & abs(frequency_residual_alt) < 5 & abs(width_residual_ref) < 20 & abs(width_residual_alt) < 20)

CA_loci_df_f1 %>% filter(TF.name == "TFBS_8564625")

partition.collapsing.dictinary = NULL
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "GenomicTile_48661719", Sample = "CTKO"), rank = 1, k = 12,
  pool.replicates = TRUE, resize.size = 500, partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "F1_bait.capture", remove.TFBS.labels = TRUE, reference.genome = BSgenome.Mmusculus.UCSC.mm10, deduplicate = FALSE, plotting.TFBSs = NULL, reutrn.chromatin.influence.df = TRUE
) -> pl
pl$pl

# D


