library(dplyr)
library(tidyverse)
library(magrittr)
library(QuasR)
library(BSgenome.Mmusculus.UCSC.mm10)
setwd("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/")
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/utils.r")
source("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/scripts/functions/source_FootprintCharter_functions.r")

chromHMM = Load.chromHMM(GenomicTiles = TRUE)
TFBSs = Load.TFBSs()
chip_qinput_pooledReps = "/g/krebs/barzaghi/HTS/ChIP-seq/2025-07-15_AV233002_Doyle_2450544031/Qinput_wasp_indel.txt"
chipped_tfs = c("Ctcf", "Klf4", "Oct4", "Sox2", "Nrf1")

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

CA_loci_df_f1 %>% 
  filter(ChIP_annotation == "bound") %>%
  full_join(., TFBS.cluster.composition.df, by = "TF.name") %>% #filter(!is.na(cluster.id)) %>% 
  mutate(cluster.id = ifelse(is.na(cluster.id), TF.name, cluster.id), nr.motifs = ifelse(is.na(nr.motifs), 1, nr.motifs)) %>%
  filter(nr.motifs == 1) %>%
  filter(Sample == "STKO") %>%
  filter(TF %in% chipped_tfs) %>%
  dplyr::select(TF.name, TF, motif.change, CA_regulatory_R, CA_regulatory_A) -> f1_motifs_smf

# F1s ChIP-seq
qprj = qAlign(
  chip_qinput_pooledReps, 
  genome = "BSgenome.Mmusculus.UCSC.mm10", 
  paired = "fr", 
  snpFile = "/g/krebs/barzaghi/DB/genetic_variation/raw_SNP_files/SPRET/SPRET_EiJ.mgp.v5.snps.dbSNP142_QuasRformat.vcf.gz"
)
Reduce(rbind, lapply(chipped_tfs, function(tf){
  
  query = TFBSs[filter(f1_motifs_smf, TF %in% tf)$TF.name] %>%
    IRanges::resize(., 2000, "center")
  
  print(paste0(tf, ": ", length(query), " motifs"))
  
  qprj_tf = qprj[unique(grep(tf, qprj@alignments$SampleName, value = FALSE, ignore.case = TRUE))]
  
  clObj = makeCluster(16)
  qcounts = qCount(qprj_tf, query, clObj = clObj)
  stopCluster(clObj)
  
  qcounts %>%
    data.frame() %>%
    rownames_to_column("TF.name") %>%
    dplyr::select(-width) %>%
    gather(sample, count, -TF.name) %>%
    mutate(TF = tf) -> df
  
  return(df)
  
})) -> ChIP_counts
qs::qsave(ChIP_counts, "ChIP_counts.qs")
ChIP_counts = qs::qread("ChIP_counts.qs")

ChIP_counts %>%
  separate(sample, c("sample", "allele"), sep = "_") %>%
  spread(allele, count) %>%
  dplyr::select(-U, -sample) %>%
  dplyr::rename("ChIPseq_count_R" = "R", "ChIPseq_count_A" = "A") %>%
  left_join(., f1_motifs_smf, by = c("TF.name", "TF")) %>%
  group_by(TF, motif.change) %>% mutate(threshold = quantile(ChIPseq_count_R, 0.05)) %>% ungroup() %>%
  mutate(motif.set = factor(case_when(
    motif.change == "no" ~ "conserved motifs",
    motif.change == "l.o.f." & ChIPseq_count_A <= threshold & ChIPseq_count_R > threshold ~ "motifs with loss of TF binding",
    motif.change == "l.o.f." ~ "allele-specific motifs"
  ), levels = c("conserved motifs", "allele-specific motifs", "motifs with loss of TF binding"))) %>%
  filter(TF %in% c("Ctcf", "Klf4", "Sox2")) %>%
  ggplot(aes(ChIPseq_count_R+1, ChIPseq_count_A+1, color = motif.set, alpha = motif.set, size = motif.set)) +
  geom_point() +
  geom_abline(color = "grey", linetype = 2) +
  facet_wrap(~TF, scales = "free") +
  scale_color_manual(values = c("grey", "black", "salmon"), breaks = c("conserved motifs", "allele-specific motifs", "motifs with loss of TF binding")) +
  scale_alpha_manual(values = c(.25, 1, 1), breaks = c("conserved motifs", "allele-specific motifs", "motifs with loss of TF binding")) +
  scale_size_manual(values = c(1, 2, 2), breaks = c("conserved motifs", "allele-specific motifs", "motifs with loss of TF binding")) +
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  theme(text = element_text(size = 16))-> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/chip_seq_f1_scatter_", Sys.Date(), ".pdf"), pl, width = 12, height = 4)

ChIP_counts %>%
  separate(sample, c("sample", "allele"), sep = "_") %>%
  spread(allele, count) %>%
  dplyr::select(-U, -sample) %>%
  dplyr::rename("ChIPseq_count_R" = "R", "ChIPseq_count_A" = "A") %>%
  left_join(., f1_motifs_smf, by = c("TF.name", "TF")) %>%
  group_by(TF, motif.change) %>% mutate(threshold = quantile(ChIPseq_count_R, 0.05)) %>% ungroup() %>%
  mutate(motif.set = factor(case_when(
    motif.change == "no" ~ "conserved motifs",
    motif.change == "l.o.f." & ChIPseq_count_A <= threshold & ChIPseq_count_R > threshold ~ "motifs with loss of TF binding"
  ), levels = c("conserved motifs", "allele-specific motifs", "motifs with loss of TF binding")))  %>%
  filter(!is.na(motif.set)) -> tmp

f1_motifs_smf %>%
  mutate(
    ChIPseq_count_A = NA,
    ChIPseq_count_R = NA,
    threshold = NA,
    motif.set = factor("allele-specific motifs", levels = c("conserved motifs", "allele-specific motifs", "motifs with loss of TF binding"))
    ) %>%
  filter(!TF.name %in% tmp$TF.name) -> tmp2

inactive.median = as.integer(quantile(filter(CA_loci_df, chromHMM_annotation == "intergenic", ChIP_annotation == "unbound")$CA_regulatory, 0.50))
rbind(tmp, tmp2) %>%
  gather(allele, CA_regulatory, CA_regulatory_R, CA_regulatory_A) %>%
  mutate(allele = gsub("CA_regulatory_", "", allele)) %>%
  mutate(allele = factor(allele, levels = c("R", "A"))) %>%
  filter(TF %in% c("Ctcf", "Klf4", "Sox2")) %>%
  ggplot(aes(motif.set, CA_regulatory, fill = allele)) +
  geom_boxplot(alpha = .5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = inactive.median, fill = viridis::mako(n=1)[1], alpha = .25) +
  ggpubr::stat_compare_means(method = "wilcox", label.x = 1.5, label = "p.format", comparisons = list(c(1,2), c(1,3))) +
  facet_wrap(~TF) +
  scale_x_discrete(labels = c(
    str_wrap("conserved motifs", width = 10), 
    str_wrap("allele-specific motifs", width = 10), 
    str_wrap("motifs with loss of TF binding", width = 10)
    )) +
  scale_fill_manual(values = c("black", "salmon"), breaks = c("R", "A")) +
  scale_y_continuous(breaks = c(0,inactive.median,50,100)) + 
  xlab("") + ylab("CA frequency % (alt - ref)") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> pl
ggplot2::ggsave(paste0("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/chip_seq_f1_", Sys.Date(), ".pdf"), pl, width = 12, height = 6)
 







