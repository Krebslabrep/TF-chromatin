library(tidyverse)
source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
source("/g/krebs/barzaghi/Rscripts/IGV_plotting.R")

chromatin.influence.df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2024-05-21_chromatin.influence_Sonmezer_SeparateReplicates.df.qs")
chromHMM = Load.chromHMM(GenomicTiles = TRUE)

# B
partition.collapsing.dictinary = split(1:16,1:16)[c(1,9,4,10,13,5,6,2,14,15,12,3,16,11,7,8)]
patch.single.site.plots(
  interpretable.master.table = data.frame(TFBS.cluster = "AMP_56", Sample = "amplicon_DE_data"),
  rank = 1, k = 16, pool.replicates = TRUE, resize.size = 500, 
  partition.collapsing.dict = partition.collapsing.dictinary, 
  data.type = "WT_amplicon", remove.TFBS.labels = TRUE, patch.version = 2
) -> pl
png(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs1b_", Sys.Date(), ".png"), width = 32, height = 25, units = "cm", res = 300)
pl$pl
dev.off()

# C
chromatin.influence.df %>%
  left_join(., chromHMM, by = "TFBS.cluster") %>%
  mutate(
    locus = factor(case_when(
      TF == "Ctcf" ~ "Ctcf",
      chromHMM_annotation == "C4_CREs_promoter" ~ "promoter",
      chromHMM_annotation == "C5_CREs_enhancer" ~ "enhancer",
      chromHMM_annotation == "C2_Inactive" ~ "inactive"
    ), levels = c("inactive", "enhancer", "promoter", "Ctcf"))) %>%
  left_join(., ChIP_thresholds_dictionary, by = "TF") %>%
  mutate(isBound = ChIP>ChIP_threshold) %>%
  filter((locus == "inactive" & !isBound) | (locus %in% c("promoter", "enhancer", "Ctcf") & isBound)) %>%
  unnest(c(acc.width.distro, acc.read.count.distro)) %>%
  # filter(acc.width.distro >= 100) %>%
  # group_by(Sample, TF.name, TF, tot.read.count, locus) %>%
  # summarise(CA.wide = sum(acc.read.count.distro), .groups = "drop") %>%
  # mutate(CA = CA.wide / tot.read.count * 100) -> pl.df
  mutate(CA = (accessible.read.count/tot.read.count)*100) -> pl.df
  
pl.df %>%
  dplyr::select(Sample, TF.name, CA) %>%
  filter(!is.na(Sample)) %>%
  spread(Sample, CA) %>%
  na.omit() %>%
  ggplot(aes(SMF_MM_TKO_R1, SMF_MM_TKO_R2)) +
  ggpointdensity::geom_pointdensity(adjust = 3) +
  ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", size = 4.5) +
  viridis::scale_color_viridis() +
  scale_x_continuous(breaks = c(0,100), limits = c(0,100)) + 
  scale_y_continuous(breaks = c(0,100), limits = c(0,100)) + 
  xlab("Replicate 1 (CA %)") +
  ylab("Replicate 2 (CA %)") +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs1c_", Sys.Date(), ".png"), pl, width = 5.5, height = 5.5, dpi = 300)

# D
chromatin.influence.df %>%
  left_join(., chromHMM, by = "TFBS.cluster") %>%
  mutate(
    locus = factor(case_when(
      TF == "Ctcf" ~ "Ctcf",
      chromHMM_annotation == "C4_CREs_promoter" ~ "promoter",
      chromHMM_annotation == "C5_CREs_enhancer" ~ "enhancer",
      chromHMM_annotation == "C2_Inactive" ~ "inactive"
    ), levels = c("inactive", "enhancer", "promoter", "Ctcf"))) %>%
  left_join(., ChIP_thresholds_dictionary, by = "TF") %>%
  mutate(isBound = ChIP>ChIP_threshold) %>%
  filter((locus == "inactive" & !isBound) | (locus %in% c("promoter", "enhancer", "Ctcf") & isBound)) %>%
  mutate(width = unlist(pmap(list(acc.width.distro, acc.read.count.distro), function(x, w){weighted.mean(x, w)}))) %>%
  dplyr::select(Sample, TF.name, width) -> pl.df

pl.df %>%
  filter(!is.na(Sample)) %>%
  spread(Sample, width) %>%
  na.omit() %>%
  ggplot(aes(SMF_MM_TKO_R1, SMF_MM_TKO_R2)) +
  ggpointdensity::geom_pointdensity(adjust = 3) +
  ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", size = 4.5) +
  viridis::scale_color_viridis() +
  scale_x_continuous(breaks = c(0,100,500), limits = c(0,500)) +
  scale_y_continuous(breaks = c(0,100,500), limits = c(0,500)) +
  xlab("Replicate 1 (width bp)") +
  ylab("Replicate 2 (width bp)") +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs1d_", Sys.Date(), ".png"), pl, width = 5.5, height = 5.5, dpi = 300)




# RESUME HERE: work your way down to the remaining panels





















chromatin.influence.df = qs::qread("/g/krebs/barzaghi/analyses/17.02.23_chromatin_influence/2023-11-09_chromatin.influence_Sonmezer_PooledReplicates.df.qs")
GenomicTiles = Load.GenomicTiles(tiles.width = 80)
TFBSs = Load.TFBSs()














##############
# RESUME HERE: I just threw the right code here but haven't remade the figures, need to refine.
#############

# C
CA_loci_df %>% 
  mutate(isBound = ChIP>ChIP_threshold) %>%
  filter((locus == "inactive" & !isBound) | (locus %in% c("promoter", "enhancer", "Ctcf") & isBound)) -> CA_loci_df_filtered

query = IRanges::resize(TFBSs[CA_loci_df_filtered$TF.name], 200, "center")
clObj = makeCluster(16)
QuasR::qCount(proj = GetQuasRprj("/g/krebs/barzaghi/HTS/DHS/Qinput.txt", BSgenome.Mmusculus.UCSC.mm10)[1:2], query = query, clObj = clObj) -> qCounts
stopCluster(clObj)

left_join(
  CA_loci_df_filtered, 
  qCounts %>%
    data.frame() %>%
    dplyr::select(-width) %>%
    rownames_to_column("TF.name"), 
  by = "TF.name"
) %>%
  mutate(DNAse.coverage = DNAse.seq_mESC_TKO) %>%
  ggplot() +
  geom_violin(aes(locus, log2(DNAse.coverage+1)), fill = "grey", draw_quantiles = c(.25, .5, .75), scale = "width") +
  xlab(NULL) +
  theme_bw() +
  theme(text = element_text(size = 18)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/figs1b_", Sys.Date(), ".pdf"), pl, width = 7, height = 5)








