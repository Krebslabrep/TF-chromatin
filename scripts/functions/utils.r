#' This masks samples uniformly...i.e. all samples will be masked for all SNPs
#'
MaskSNPs = function(Methylation, CytosinesToMask, MaskSMmat = FALSE, Experiment){
  
  if(Experiment == "DE"){
    message("Masking GRanges in DE mode")
    if (MaskSMmat){
      Cast_DE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]], CytosinesToMask[CytosinesToMask$DisruptedInCast])))
      Spret_DE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]], CytosinesToMask[CytosinesToMask$DisruptedInSpret])))
      all_DE_disrupted = c(Cast_DE_disrupted, Spret_DE_disrupted) %>% sort() %>% unique()
      if(length(all_DE_disrupted)>0){elementMetadata(Methylation[[1]])[all_DE_disrupted,] = NA}
      
      message("Masking SM matrix")
      DisruptedCoords = start(Methylation[[1]])[all_DE_disrupted]
      Methylation[[2]] = lapply(Methylation[[2]], function(mat){mat[,!colnames(mat) %in% DisruptedCoords, drop=FALSE]})
      
    } else {
      message("Skipping SM matrix")
      Cast_DE_disrupted = unique(queryHits(findOverlaps(Methylation, CytosinesToMask[CytosinesToMask$DisruptedInCast])))
      Spret_DE_disrupted = unique(queryHits(findOverlaps(Methylation, CytosinesToMask[CytosinesToMask$DisruptedInSpret])))
      all_DE_disrupted = c(Cast_DE_disrupted, Spret_DE_disrupted) %>% sort() %>% unique()
      if(length(all_DE_disrupted)>0){elementMetadata(Methylation)[all_DE_disrupted,] = NA}
    }
    
  } else if (Experiment == "NO") {
    message("Masking GRanges in NO mode")
    for(context in seq(2)){
      if (MaskSMmat){
        Cast_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]][[context]], CytosinesToMask[CytosinesToMask$DisruptedInCast])))
        Spret_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]][[context]], CytosinesToMask[CytosinesToMask$DisruptedInSpret])))
        all_SE_disrupted = c(Cast_SE_disrupted, Spret_SE_disrupted) %>% sort() %>% unique()
        if(length(all_SE_disrupted)>0){elementMetadata(Methylation[[1]][[context]])[all_SE_disrupted,] = NA}
        
        message("Masking SM matrix")
        DisruptedCoords = start(Methylation[[1]][[context]])[all_SE_disrupted]
        # CastMats = grep(Castaneus.string.match, grep("_MethRate$", colnames(elementMetadata(Methylation[[1]][[context]])), value = TRUE))
        # SpretMats = grep(Spretus.string.match, grep("_MethRate$", colnames(elementMetadata(Methylation[[1]][[context]])), value = TRUE))
        
        for(sample in seq_along(Methylation[[2]])){
          Coordinates_to_keep = !colnames(Methylation[[2]][[sample]][[context]]) %in% DisruptedCoords
          Methylation[[2]][[sample]][[context]] = Methylation[[2]][[sample]][[context]][,Coordinates_to_keep, drop=FALSE]
        }
 
      } else {
        message("Skipping SM matrix")
        Cast_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[context]], CytosinesToMask[CytosinesToMask$DisruptedInCast])))
        Spret_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[context]], CytosinesToMask[CytosinesToMask$DisruptedInSpret])))
        all_SE_disrupted = c(Cast_SE_disrupted, Spret_SE_disrupted) %>% sort() %>% unique()
        if(length(all_SE_disrupted)>0){elementMetadata(Methylation[[context]])[all_SE_disrupted,] = NA}
      }
    }
  } else {stop("Unrecognized Experiment")}
  
  return(Methylation)
  
}

#' This masks samples distinctly...only relevant SNPs to species
#'
MaskSNPs2 = function(Methylation, CytosinesToMaks, MaskSMmat = FALSE, Experiment){
  
  Castaneus.string.match = c("Cast|Barzaghi_C|_CTKO")
  Spretus.string.match = c("Spret|Barzaghi_S|STKmixed|guidosp|_STKO")
  
  if(Experiment == "DE"){
    message("Masking GRanges in DE mode")
    if (MaskSMmat){
      Cast_DE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]], CytosinesToMaks[CytosinesToMaks$DisruptedInCast])))
      Spret_DE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]], CytosinesToMaks[CytosinesToMaks$DisruptedInSpret])))
      CastCols = grep(Castaneus.string.match, colnames(elementMetadata(Methylation[[1]])))
      SpretCols = grep(Spretus.string.match, colnames(elementMetadata(Methylation[[1]])))
      if(length(CastCols)>0){elementMetadata(Methylation[[1]])[Cast_DE_disrupted,CastCols] = NA}
      if(length(SpretCols)>0){elementMetadata(Methylation[[1]])[Spret_DE_disrupted,SpretCols] = NA}
      
      message("Masking SM matrix")
      CastDisruptedCoords = start(Methylation[[1]])[Cast_DE_disrupted]
      SpretDisruptedCoords = start(Methylation[[1]])[Spret_DE_disrupted]
      CastMats = grep(Castaneus.string.match, grep("_MethRate$", colnames(elementMetadata(Methylation[[1]])), value = TRUE))
      SpretMats = grep(Spretus.string.match, grep("_MethRate$", colnames(elementMetadata(Methylation[[1]])), value = TRUE))
      
      Methylation[[2]][CastMats] = lapply(Methylation[[2]][CastMats], function(mat){mat[,!colnames(mat) %in% CastDisruptedCoords, drop=FALSE]})
      Methylation[[2]][SpretMats] = lapply(Methylation[[2]][SpretMats], function(mat){mat[,!colnames(mat) %in% SpretDisruptedCoords, drop=FALSE]})
      
    } else {
      message("Skipping SM matrix")
      Cast_DE_disrupted = unique(queryHits(findOverlaps(Methylation, CytosinesToMaks[CytosinesToMaks$DisruptedInCast])))
      Spret_DE_disrupted = unique(queryHits(findOverlaps(Methylation, CytosinesToMaks[CytosinesToMaks$DisruptedInSpret])))
      CastCols = grep(Castaneus.string.match, colnames(elementMetadata(Methylation)))
      SpretCols = grep(Spretus.string.match, colnames(elementMetadata(Methylation)))
      if(length(CastCols)>0){elementMetadata(Methylation)[Cast_DE_disrupted,CastCols] = NA}
      if(length(SpretCols)>0){elementMetadata(Methylation)[Spret_DE_disrupted,SpretCols] = NA}
    }
    
  } else if (Experiment == "NO") {
    message("Masking GRanges in NO mode")
    for(context in seq(2)){
      if (MaskSMmat){
        Cast_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]][[context]], CytosinesToMaks[CytosinesToMaks$DisruptedInCast])))
        Spret_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]][[context]], CytosinesToMaks[CytosinesToMaks$DisruptedInSpret])))
        CastCols = grep(Castaneus.string.match, colnames(elementMetadata(Methylation[[1]][[context]])))
        SpretCols = grep(Spretus.string.match, colnames(elementMetadata(Methylation[[1]][[context]])))
        if(length(CastCols)>0){elementMetadata(Methylation[[1]][[context]])[Cast_SE_disrupted,CastCols] = NA}
        if(length(SpretCols)>0){elementMetadata(Methylation[[1]][[context]])[Spret_SE_disrupted,SpretCols] = NA}
        
        message("Masking SM matrix")
        CastDisruptedCoords = start(Methylation[[1]][[context]])[Cast_SE_disrupted]
        SpretDisruptedCoords = start(Methylation[[1]][[context]])[Spret_SE_disrupted]
        CastMats = grep(Castaneus.string.match, grep("_MethRate$", colnames(elementMetadata(Methylation[[1]][[context]])), value = TRUE))
        SpretMats = grep(Spretus.string.match, grep("_MethRate$", colnames(elementMetadata(Methylation[[1]][[context]])), value = TRUE))
        
        for(sample in CastMats){
          Coordinates_to_keep = !colnames(Methylation[[2]][[sample]][[context]]) %in% CastDisruptedCoords
          Methylation[[2]][[sample]][[context]] = Methylation[[2]][[sample]][[context]][,Coordinates_to_keep, drop=FALSE]
        }
        
        for(sample in SpretMats){
          Coordinates_to_keep = !colnames(Methylation[[2]][[sample]][[context]]) %in% SpretDisruptedCoords
          Methylation[[2]][[sample]][[context]] = Methylation[[2]][[sample]][[context]][,Coordinates_to_keep, drop=FALSE]
        }
      } else {
        message("Skipping SM matrix")
        Cast_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[context]], CytosinesToMaks[CytosinesToMaks$DisruptedInCast])))
        Spret_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[context]], CytosinesToMaks[CytosinesToMaks$DisruptedInSpret])))
        CastCols = grep(Castaneus.string.match, colnames(elementMetadata(Methylation[[context]])))
        SpretCols = grep(Spretus.string.match, colnames(elementMetadata(Methylation[[context]])))
        if(length(CastCols)>0){elementMetadata(Methylation[[context]])[Cast_SE_disrupted,CastCols] = NA}
        if(length(SpretCols)>0){elementMetadata(Methylation[[context]])[Spret_SE_disrupted,SpretCols] = NA}
      }
    }
  } else {stop("Unrecognized Experiment")}
  
  return(Methylation)
  
}

LoadSNPs = function(){
  SNPs = qs::qread("/g/krebs/barzaghi/DB/genetic_variation/raw_SNP_files/SNPs.qs")
  return(SNPs)
}

Load.CytosinesToMask = function(){
  CytosinesToMask = qs::qread("/g/krebs/barzaghi/DB/CytosinesToMask.qs")
  return(CytosinesToMask)
}

#'
#' @param tiles.width 80 or 120
Load.GenomicTiles = function(tiles.width, GeneAssociation = FALSE, CGI = FALSE, CpG_density = FALSE){
  
  if(tiles.width == 80){
    if(GeneAssociation){
      GenomicTiles = readRDS("/g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/GeneExpressionAnnotated_GenomicTiles.RDS")
    } else{
      GenomicTiles = qs::qread("/g/krebs/barzaghi/analyses/10.05.2022_create.master.table/11.07.22_run.all.accessible.tiles/GenomicTiles_width80_step_40_EnrichmentRegions.qs")
    }
  }
  
  if(tiles.width == 120){
    GenomicTiles = readRDS("/g/krebs/barzaghi/DB/GenomicTiles_width120_step_60_EnrichmentRegions.rds")
  }
  
  if (CGI){
    CGI = annotatr::build_annotations(genome = "mm10", annotations = "mm10_cpg_islands")
    GenomicTiles$CGI = FALSE
    GenomicTiles$CGI[unique(queryHits(GenomicRanges::findOverlaps(GenomicTiles, CGI)))] = TRUE
  }
  
  if (CpG_density){
    sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, GenomicTiles)
    CpG_counts = vcountPattern("CG", subject = sequences)
    CpG_densities <- CpG_counts/lengths(sequences)
    GenomicTiles$CpG_density = CpG_densities
  }
  
  return(GenomicTiles)
  
}

#' @param kind stringent/lenient/Sonmezer/Chen/Barzaghi
#' @param GenVar.info none/GenVar/GenVarChange
#' @param idx.by ""/TF/Seqnames
#' @param reduced.overlapping reduces overlapping binding sites for the same TFs. Picks the highest scoring match regardless of strand
#' 
Load.TFBSs = function(kind = NULL, GenVar.info = NULL, idx.by = NULL, reduced.overlapping = NULL, custom.file = NULL){
  
  # I'm not using anymore anything different
  if(!is.null(custom.file)){
    warning("working exclusively with qs files")
    TFBSs = qs::qread(custom.file)
  } else {
    TFBSs = qs::qread("/g/krebs/barzaghi/DB/TFBSs.qs")
  }
  
  return(TFBSs)
  
}

Load.Sonmezer.amplicon.GRanges = function(){
  # These come from the supplementary of Sönmezer et al
  xlsx::read.xlsx("/g/krebs/barzaghi/analyses/Can_Amplicons/1-s2.0-S1097276520307930-mmc3.xlsx", sheetIndex = 1) %>%
    dplyr::select(ID, chr, start, end) %>%
    dplyr::rename(seqnames = "chr") %>%
    GRanges() -> all.plates.coords
  all.plates.coords$Forward.Primer.Sequence = xlsx::read.xlsx("/g/krebs/barzaghi/analyses/Can_Amplicons/1-s2.0-S1097276520307930-mmc3.xlsx", sheetIndex = 1)$Forward.Primer.Sequence
  all.plates.coords$Forward.Primer.Sequence.1 = xlsx::read.xlsx("/g/krebs/barzaghi/analyses/Can_Amplicons/1-s2.0-S1097276520307930-mmc3.xlsx", sheetIndex = 1)$Forward.Primer.Sequence.1
  xlsx::read.xlsx("/g/krebs/barzaghi/analyses/Can_Amplicons/1-s2.0-S1097276520307930-mmc2.xlsx", sheetIndex = 1) -> all.plates.annotation
  
  all.plates.coords$Plate_ID = all.plates.annotation$Plate_ID
  all.plates.coords$TargetType = all.plates.annotation$TargetType
  all.plates.coords$TargetTF = all.plates.annotation$TargetTF
  all.plates.coords$TargetSite = all.plates.annotation$TargetSite
  
  names(all.plates.coords) = paste0("AMP_", seq_along(all.plates.coords))
  
  return(all.plates.coords)
}

Load.chromHMM = function(GenomicTiles = FALSE){
  
  if(!GenomicTiles){
    chromHMM = rtracklayer::import.bed("/g/krebs/barzaghi/DB/mESC_E14_12_dense.annotated.bed")
    chromHMM %<>%
      plyranges::mutate(chromHMM = case_when(
        name == "1_Insulator" ~ "Ctcf",
        name == "2_Intergenic" ~ "intergenic",
        name == "3_Heterochromatin" ~ "heterochromatin",
        name == "4_Enhancer" ~ "enhancer",
        name == "5_RepressedChromatin" ~ "repressed",
        name == "6_BivalentChromatin" ~ "bivalent_promoter",
        name == "7_ActivePromoter" ~ "promoter",
        name == "8_StrongEnhancer" ~ "enhancer",
        name == "9_TranscriptionTransition" ~ "transcription",
        name == "10_TranscriptionElongation" ~ "transcription",
        name == "11_WeakEnhancer" ~ "enhancer",
        name == "12_LowSignal/RepetitiveElements" ~ "low_signal"
      ))
  } else {
    chromHMM = qs::qread("/g/krebs/barzaghi/DB/mESC_E14_12_dense.annotated_GenomicTiles.qs")
  }
  
  return(chromHMM)
  
}

# data.frame(
#   FileName = c(
#     grep(".*bam$", list.files("/g/krebs/barzaghi/analyses/nf-core_runs/261023_Banp_ChIPseq/bwa/mergedLibrary", "tko_Banp.*bam", full.name = TRUE, recursive = TRUE), value = TRUE)
#     # grep(".*bam$", list.files("/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/Zic3/", "", full.name = TRUE, recursive = TRUE), value = TRUE)
#     # grep(".*bam$", list.files("/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/ESC/", "", full.name = TRUE, recursive = TRUE), value = TRUE)
#     ),
#   SampleName = "Banp"
#   ) %>% readr::write_delim(., file = "/g/krebs/barzaghi/Tmp/QuasR_input_Banp_tko_ChIPseq.txt", delim = "\t")

# N.b.: for each TF I picked the 5 largest files so that the estimation of unbound would be comparable across TFs
ChIP_data_dictionary = list(
  Oct4 = "/g/krebs/barzaghi/Tmp/QuasR_input_Oct4_ChIPnexus.txt",
  Ctcf = "/g/krebs/barzaghi/Tmp/QuasR_input_Ctcf_ChIPseq.txt",
  Yy1 = "/g/krebs/barzaghi/Tmp/QuasR_input_Yy1_ChIPseq.txt",
  Rest = "/g/krebs/barzaghi/Tmp/QuasR_input_Rest_ChIPseq.txt",
  Klf4 = "/g/krebs/barzaghi/Tmp/QuasR_input_Klf4_ChIPnexus.txt",
  Sox2 = "/g/krebs/barzaghi/Tmp/QuasR_input_Sox2_ChIPnexus.txt",
  Nanog = "/g/krebs/barzaghi/Tmp/QuasR_input_Nanog_ChIPnexus.txt",
  Esrrb = "/g/krebs/barzaghi/Tmp/QuasR_input_Esrrb_ChIPseq.txt",
  Nrf1 = "/g/krebs/barzaghi/Tmp/QuasR_input_Nrf1_ChIPseq.txt",
  Zic3 = "/g/krebs/barzaghi/Tmp/QuasR_input_Zic3_ChIPseq.txt",
  `Bach1::Mafk` = "/g/krebs/barzaghi/Tmp/QuasR_input_Bach1::Mafk_ChIPseq.txt",
  Banp = "/g/krebs/barzaghi/Tmp/QuasR_input_Banp_tko_ChIPseq.txt",
  Nfya = "/g/krebs/barzaghi/Tmp/QuasR_input_Nfya_ChIPseq.txt",
  `Smad2::Smad3::Smad4` = "/g/krebs/barzaghi/Tmp/QuasR_input_Smad2::Smad3::Smad4_ChIPseq.txt",
  Sox21 = NULL, 
  Tbp = NULL, 
  Tgif1 = NULL, 
  Tgif2 = NULL, 
  Tp53 = NULL,
  Foxd3 = "/g/krebs/barzaghi/Tmp/QuasR_input_Foxd3_ChIPseq.txt",
  E2f1 = "/g/krebs/barzaghi/Tmp/QuasR_input_E2f1_ChIPseq.txt",
  Stat3 = "/g/krebs/barzaghi/Tmp/QuasR_input_Stat3_ChIPseq.txt",
  Myc = "/g/krebs/barzaghi/Tmp/QuasR_input_Myc_ChIPnexus.txt",
  Otx2 = "/g/krebs/barzaghi/Tmp/QuasR_input_Otx2_ChIPseq.txt",
  Zfx = "/g/krebs/barzaghi/Tmp/QuasR_input_Zfx_ChIPseq.txt",
  H3K27Ac = "/g/krebs/barzaghi/Tmp/QuasR_input_H3K27Ac_ChIPseq.txt",
  H3K4me3 = "/g/krebs/barzaghi/Tmp/QuasR_input_H3K4me3_ChIPseq.txt"
)

ChIP_thresholds_dictionary = data.frame(
  `Bach1::Mafk` = 50,
  Banp = 50,
  Ctcf = 150,
  E2f1 = 50,
  Esrrb = 50,
  Foxd3 = 50,
  Klf4 = 100,
  Myc = 40,
  Nanog = NA,
  Nfya = 50,
  Nrf1 = 50,
  Oct4 = 100,
  Otx2 = 100,
  Rest = 150,
  `Smad2::Smad3::Smad4` = 50,
  Sox2 = 40,
  Sox21 = NA, 
  Stat3 = 10,
  Tbp = NA, 
  Tgif1 = NA, 
  Tgif2 = NA, 
  Tp53 = NA,
  Yy1 = 150,
  Zfx = 40,
  Zic3 = 50
) %>% gather(TF, ChIP_threshold) %>%
  mutate(TF = gsub("\\.\\.", "::", TF))

ChIP_thresholds_dictionary_lenient = data.frame(
  `Bach1::Mafk` = 10,
  Banp = 50,
  Ctcf = 75,
  E2f1 = 15,
  Esrrb = 20,
  Foxd3 = 35,
  Klf4 = 50,
  Myc = 40,
  Nanog = NA,
  Nfya = 10,
  Nrf1 = 40,
  Oct4 = 50,
  Otx2 = 20,
  Rest = 15,
  `Smad2::Smad3::Smad4` = 50,
  Sox2 = 30,
  Sox21 = NA, 
  Stat3 = 10,
  Tbp = NA, 
  Tgif1 = NA, 
  Tgif2 = NA, 
  Tp53 = NA,
  Yy1 = 30,
  Zfx = 10,
  Zic3 = 10
) %>% gather(TF, ChIP_threshold) %>%
  mutate(TF = gsub("\\.\\.", "::", TF))

jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

panel.jet <- function(...) {
  smoothScatter(..., nrpoints=0, add=TRUE, colramp=jet.colors) }

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="grey", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 2.5/strwidth(txt)
  text(0.5, 0.5, txt)
}
