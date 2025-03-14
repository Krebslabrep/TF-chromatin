library(tidyverse)
library(GenomicRanges)
library(scales)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(plyranges)
library(rtracklayer)
library(ggbio)
library(BSgenome.Mmusculus.UCSC.mm10)

plot_genomic_track <- function(sampleSheet, samples, allelic = FALSE, RegionOfInterest, tile.width = 200, tile.step = 50, max.y.lim = NULL, color = "black", y.labs = "", delta = FALSE, normalise = FALSE, plot.coordinates = FALSE){
  
  proj = QuasR::qAlign(
    sampleFile = sampleSheet, genome = "BSgenome.Mmusculus.UCSC.mm10", 
    paired = "fr", aligner = "Rbowtie", snpFile = if(allelic){"/g/krebs/barzaghi/Tmp/1.pdf"}else{NULL})
  proj = proj[which(proj@alignments$SampleName %in% samples)]
  proj@aligner = "Rbowtie"
  target_ranges = plyranges::slide_ranges(RegionOfInterest, width = tile.width, step = tile.step)
  names(target_ranges) = paste0(names(target_ranges), "_", seq_along(target_ranges))
  counts = QuasR::qCount(proj = proj, query = target_ranges)
  
  if(is.null(max.y.lim)){max.y.lim = max(counts[,-1])}
  
  labeller.vector = y.labs
  names(labeller.vector) = grep("_U$", colnames(counts[,-1,drop=FALSE]), invert = TRUE, value = TRUE)
  
  cbind(
    target_ranges %>%
      IRanges::resize(., 1, "center") %>%
      as.data.frame() %>%
      dplyr::select(start),
    as.data.frame(counts[,-1,drop=FALSE])
  ) %>%
    gather(sample, score, -start) %>% 
    filter(case_when(isTRUE(allelic) ~ str_detect(sample, "_U$", negate = TRUE),
                     isFALSE(allelic) ~ score >= 0) # mock filter to let the case_when pass
           ) -> pl.df
  
  if(allelic & delta){
    
    pl.df %<>%
      mutate(allele = ifelse(str_detect(sample, "_A$"), "A", "R")) %>%
      mutate(sample = gsub("_A$|_R$", "", sample)) %>%
      spread(allele, score) %>%
      mutate(score = A-R)
    
    labeller.vector = c(labeller.vector, paste0("delta ", paste(y.labs, collapse = " ")))
    names(labeller.vector)[3] = unique(pl.df$sample)
    
  } else if (allelic & !delta) {
    
    factor.levels = c(unique(grep(".*_R$", pl.df$sample, value = TRUE)), unique(grep(".*_A$", pl.df$sample, value = TRUE)))
    pl.df %<>%
      mutate(sample = factor(sample, levels = factor.levels))
    
  }
  
  if(plot.coordinates){
    x.axis.breaks = c(start(RegionOfInterest), end(RegionOfInterest))
    x.axis.labels = c(start(RegionOfInterest), end(RegionOfInterest))
  } else {
    x.axis.breaks = c(start(RegionOfInterest), end(RegionOfInterest))
    x.axis.labels = c("", "")
  }
  
  pl.df %>%
    ggplot(aes(start, score)) +
    geom_vline(xintercept=start(IRanges::resize(RegionOfInterest, width = 1, fix = "center")), color = "grey", linetype = 2) +
    geom_area(aes(fill = sample)) +
    geom_text(data = data.frame(x=NA), aes(x = end(RegionOfInterest)-200, y = max.y.lim*.75, label = y.labs), hjust = "right", inherit.aes = FALSE) +
    facet_wrap(~sample, ncol = 1, strip.position = "right", labeller = labeller(sample = labeller.vector)) +
    scale_fill_manual(values = color) + 
    scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.labels, nsmall=1, big.mark=","), expand = c(0, 0)) +
    scale_y_continuous(limits = c(ifelse(delta, -max.y.lim, 0), max.y.lim), breaks = c(ifelse(delta, -max.y.lim, 0), max.y.lim), expand = expansion(mult = 0, add = 0)) +
    theme_bw() +
    theme(
      strip.background = element_blank(), strip.text = element_blank(), axis.title = element_blank(),
      legend.position = "none",
      axis.line.y = element_line(linewidth = 1),
      axis.ticks.length = unit(0, "mm"),
      panel.grid = element_blank(), panel.border = element_blank()
      )
  
}