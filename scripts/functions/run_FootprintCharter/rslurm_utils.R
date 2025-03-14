# Rscript_path = "/g/funcgen/bin/Rscript-4.2.0"
Rscript_path = "/g/easybuild/x86_64/Rocky/8/haswell/software/R/4.2.2-foss-2022b/bin/Rscript"
lib_paths = c( # the first and last are already there
  # "/g/krebs/barzaghi/R-lib/x86_64-pc-linux-gnu/4.2.2",
  "/g/easybuild/x86_64/Rocky/8/haswell/software/R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2",
  "/g/easybuild/x86_64/Rocky/8/haswell/software/arrow-R/11.0.0.3-foss-2022b-R-4.2.2"
  # "/g/easybuild/x86_64/Rocky/8/haswell/software/R/4.2.2-foss-2022b/lib64/R/library"
)

#### UnsupervisedClustering ####

rslurm_UnsupervisedClustering = function(indexes, sampleSheet, Load.TFBSs.args, genetic.variation, methylation, strong.binders.only, tryCatchLog.write.error.dump.file){
  
  source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
  source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
  
  sampleSheet = sampleSheet
  MySamples = list(NO = grep("NO", readr::read_delim(sampleSheet, "\t", show_col_types = FALSE)$SampleName, value = TRUE) %>% unique(),
                   DE = grep("DE", readr::read_delim(sampleSheet, "\t", show_col_types = FALSE)$SampleName, value = TRUE) %>% unique())
  if(genetic.variation){
    SNPs = LoadSNPs()
    CytosinesToMask = Load.CytosinesToMask()
  } else {
    SNPs = NULL
    CytosinesToMask = NULL
  }
  
  TFBSs = Load.TFBSs(kind = Load.TFBSs.args$kind, GenVar.info = Load.TFBSs.args$GenVar.info, idx.by = Load.TFBSs.args$idx.by)
  if (strong.binders.only){
    TFBSs.for.CallingWindows = TFBSs[TFBSs$TF %in% strong.binders()]
  } else {
    TFBSs.for.CallingWindows = TFBSs
  }
  GenomicTiles = Load.GenomicTiles(tiles.width = 80, GeneAssociation = FALSE, CGI = FALSE)
  seqlevels(GenomicTiles) = seqlevels(TFBSs)
  MethylationCallingWindows = SingleMoleculeFootprinting::Create_MethylationCallingWindows(
    TFBS_cluster_coordinates = GenomicTiles, 
    genomic.seqlenghts = GenomeInfoDb::seqlengths(BSgenome.Mmusculus.UCSC.mm10), 
    fix.window.size = TRUE, max.window.size = 50)
  CurrentGenomicTiles = plyranges::filter_by_overlaps(
    GenomicTiles, MethylationCallingWindows[unlist(indexes)], minoverlap = 80
    #       N.b. minoverlap ensures that the MCWs computed within <Unsupervised.clustering.MultiSite.wrapper> are exact and that each GenomicTile is analyised exactly once
    )
  # <-> #
  
  # # PREPARE INPUT #
  TFBS.clusters = list()
  TFBS.clusters$ClusterCoordinates = CurrentGenomicTiles
  TFBS.clusters$ClusterComposition = NULL
  Unsupervised.clustering.MultiSite.wrapper(
    sampleSheet = sampleSheet,
    Samples = MySamples[[methylation]],
    pool.samples = TRUE,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    TFBS.clusters = TFBS.clusters,
    MethylationCallingWindows = MethylationCallingWindows[unlist(indexes)],
    SNPs = SNPs,
    TFBSs = TFBSs, 
    CytosinesToMask = CytosinesToMask,
    species.string.match = "_C",
    params = bait.capture.parameters.list,
    cores = 8, 
    tryCatchLog.write.error.dump.file = tryCatchLog.write.error.dump.file
  ) -> UnsupervisedClusteringResults
  return(UnsupervisedClusteringResults)

}

launch_rslurm_UnsupervisedClustering = function(working.directory, sampleSheet, Load.TFBSs.args, genetic.variation=TRUE, methylation = "DE", strong.binders.only, additional.jobname.info, NR.NODES = 1600, tryCatchLog.write.error.dump.file){
  
  setwd(working.directory)
  
  message("Loading TFBSs ... ")
  TFBSs = Load.TFBSs(kind = Load.TFBSs.args$kind, GenVar.info = Load.TFBSs.args$GenVar.info, idx.by = Load.TFBSs.args$idx.by)
  if (strong.binders.only){
    TFBSs.for.CallingWindows = TFBSs[TFBSs$TF %in% strong.binders()]
  } else {
    TFBSs.for.CallingWindows = TFBSs
  }
  GenomicTiles = Load.GenomicTiles(tiles.width = 80, GeneAssociation = FALSE, CGI = FALSE)
  seqlevels(GenomicTiles) = seqlevels(TFBSs)
  MethylationCallingWindows = SingleMoleculeFootprinting::Create_MethylationCallingWindows(
    TFBS_cluster_coordinates = GenomicTiles, 
    genomic.seqlenghts = GenomeInfoDb::seqlengths(BSgenome.Mmusculus.UCSC.mm10), 
    fix.window.size = TRUE, 
    max.window.size = 50
    )
  
  l = length(MethylationCallingWindows)
  params.df = data.frame(NODE = seq(NR.NODES))
  data.frame(start.index = ceiling(seq(1,l, l/NR.NODES)), end.index = ceiling(seq(l/NR.NODES,l,l/NR.NODES))) -> df
  params.df$indexes = lapply(seq(NR.NODES), function(i){seq(df$start.index[i], df$end.index[i])})
  params.df$sampleSheet = sampleSheet
  params.df$Load.TFBSs.args = list(Load.TFBSs.args)
  params.df$genetic.variation = genetic.variation
  params.df$methylation = methylation
  params.df$strong.binders.only = strong.binders.only
  params.df$NODE = NULL
  params.df$tryCatchLog.write.error.dump.file = tryCatchLog.write.error.dump.file
  
  slurm.options <- list(time = '03:00:00', mem = "200G", qos = "normal", error = "slurm.%N.%j.err")#, partition = "bigmem")
  
  jobname = paste0("UnsupervisedClustering_", Load.TFBSs.args$kind, ifelse(strong.binders.only, "_strong.binders.only", ""), additional.jobname.info)
  rslurm::slurm_apply(
    rslurm_UnsupervisedClustering, params.df, 
    nodes = NR.NODES, cpus_per_node = 8, jobname = jobname, rscript_path = Rscript_path,
    preschedule_cores = TRUE, slurm_options = slurm.options, libPaths = lib_paths,
    submit = FALSE) -> slurm.object
  
  return(slurm.object)
  
}

rslurm_UnsupervisedClustering_PRA = function(indexes, sampleSheet, Load.TFBSs.args, strong.binders.only, tryCatchLog.write.error.dump.file){
  
  source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
  source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
  
  sampleSheet = sampleSheet
  MySamples = list(NO = grep("NO", readr::read_delim(sampleSheet, "\t", show_col_types = FALSE)$SampleName, value = TRUE) %>% unique(),
                   DE = grep("DE", readr::read_delim(sampleSheet, "\t", show_col_types = FALSE)$SampleName, value = TRUE) %>% unique())
  
  TFBSs = Load.TFBSs(kind = Load.TFBSs.args$kind, GenVar.info = Load.TFBSs.args$GenVar.info, idx.by = Load.TFBSs.args$idx.by)
  GenomicTiles = qs::qread("/g/krebs/barzaghi/analyses/05.04.23_Unsupervised_on_PRA/CRE_autonomy_final_ranges.qs")
  seqlevels(GenomicTiles) = seqlevels(TFBSs)
  MethylationCallingWindows = SingleMoleculeFootprinting::Create_MethylationCallingWindows(
    TFBS_cluster_coordinates = GenomicTiles, 
    genomic.seqlenghts = GenomeInfoDb::seqlengths(BSgenome.Mmusculus.UCSC.mm10), 
    fix.window.size = TRUE, max.window.size = 50)
  CurrentGenomicTiles = plyranges::filter_by_overlaps(
    GenomicTiles, MethylationCallingWindows[unlist(indexes)], minoverlap = min(width(GenomicTiles))
    #       N.b. minoverlap ensures that the MCWs computed within <Unsupervised.clustering.MultiSite.wrapper> are exact and that each GenomicTile is analyised exactly once
    #       WARNING: the above comment isn't true here. But beucase it's amplicons this shouldn't be a problem
  )
  # <-> #
  
  # # PREPARE INPUT #
  TFBS.clusters = list()
  TFBS.clusters$ClusterCoordinates = CurrentGenomicTiles
  TFBS.clusters$ClusterComposition = NULL
  
  amplicon.parameters.list$CallContextMethylation[["coverage"]] = 20
  amplicon.parameters.list$SortReadsByTFCluster[["coverage"]] = 20
  amplicon.parameters.list[["unsupervised.clustering.coverage"]] = 20
    
  Unsupervised.clustering.MultiSite.wrapper(
    sampleSheet = sampleSheet,
    Samples = MySamples$DE,
    pool.samples = TRUE,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    MethylationCallingWindows = MethylationCallingWindows[unlist(indexes)],
    TFBS.clusters = TFBS.clusters,
    SNPs = NULL,
    TFBSs = TFBSs,
    CytosinesToMask = NULL,
    species.string.match = "_C",
    deduplicate = TRUE,
    params = amplicon.parameters.list,
    cores = 8,
    tryCatchLog.write.error.dump.file = tryCatchLog.write.error.dump.file
  ) -> UnsupervisedClusteringResults
  
  return(UnsupervisedClusteringResults)
  
}

launch_rslurm_UnsupervisedClustering_PRA = function(working.directory, sampleSheet, Load.TFBSs.args, strong.binders.only, additional.jobname.info, NR.NODES = 64, tryCatchLog.write.error.dump.file){
  
  setwd(working.directory)
  
  message("Loading TFBSs ... ")
  TFBSs = Load.TFBSs(kind = Load.TFBSs.args$kind, GenVar.info = Load.TFBSs.args$GenVar.info, idx.by = Load.TFBSs.args$idx.by)
  GenomicTiles = qs::qread("/g/krebs/barzaghi/analyses/05.04.23_Unsupervised_on_PRA/CRE_autonomy_final_ranges.qs")
  seqlevels(GenomicTiles) = seqlevels(TFBSs)
  MethylationCallingWindows = SingleMoleculeFootprinting::Create_MethylationCallingWindows(
    TFBS_cluster_coordinates = GenomicTiles, 
    genomic.seqlenghts = GenomeInfoDb::seqlengths(BSgenome.Mmusculus.UCSC.mm10), 
    fix.window.size = TRUE, 
    max.window.size = 50
  )
  
  l = length(MethylationCallingWindows)
  params.df = data.frame(NODE = seq(NR.NODES))
  data.frame(start.index = ceiling(seq(1,l, l/NR.NODES)), end.index = ceiling(seq(l/NR.NODES,l,l/NR.NODES))) -> df
  params.df$indexes = lapply(seq(NR.NODES), function(i){seq(df$start.index[i], df$end.index[i])})
  params.df$sampleSheet = sampleSheet
  params.df$Load.TFBSs.args = list(Load.TFBSs.args)
  params.df$strong.binders.only = strong.binders.only
  params.df$NODE = NULL
  params.df$tryCatchLog.write.error.dump.file = tryCatchLog.write.error.dump.file
  
  slurm.options <- list(time = '03:00:00', mem = "200G", qos = "normal", error = "slurm.%N.%j.err")
  
  jobname = paste0("UnsupervisedClustering_PRA_", Load.TFBSs.args$kind, ifelse(strong.binders.only, "_strong.binders.only", ""), additional.jobname.info)
  rslurm::slurm_apply(
    rslurm_UnsupervisedClustering_PRA, params.df, 
    nodes = NR.NODES, cpus_per_node = 8, jobname = jobname, rscript_path = Rscript_path,
    preschedule_cores = TRUE, slurm_options = slurm.options, libPaths = lib_paths, 
    submit = FALSE) -> slurm.object
  
  return(slurm.object)
  
}

rslurm_UnsupervisedClustering_PRA_mutants = function(indexes, sampleSheet, Load.TFBSs.args, strong.binders.only, tryCatchLog.write.error.dump.file){
  
  library(BSgenome.Mmusculus.KrebsEMBL.MutantsBatch1VB)
  source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
  source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
  
  sampleSheet = sampleSheet
  MySamples = list(NO = grep("NO", readr::read_delim(sampleSheet, "\t", show_col_types = FALSE)$SampleName, value = TRUE) %>% unique(),
                   DE = grep("DE", readr::read_delim(sampleSheet, "\t", show_col_types = FALSE)$SampleName, value = TRUE) %>% unique())
  
  TFBSs = qs::qread("/g/krebs/barzaghi/analyses/08.02.24_RMCE_mutants/TFBSs.qs")
  GenomicTiles = qs::qread("/g/krebs/barzaghi/analyses/08.02.24_RMCE_mutants/MutantsBatch1VB_ranges.qs")
  seqlevels(GenomicTiles) = seqlevels(TFBSs)
  MethylationCallingWindows = SingleMoleculeFootprinting::Create_MethylationCallingWindows(
    TFBS_cluster_coordinates = GenomicTiles, 
    genomic.seqlenghts = GenomeInfoDb::seqlengths(BSgenome.Mmusculus.KrebsEMBL.MutantsBatch1VB), 
    fix.window.size = TRUE, max.window.size = 50)
  CurrentGenomicTiles = plyranges::filter_by_overlaps(
    GenomicTiles, MethylationCallingWindows[unlist(indexes)], minoverlap = min(width(GenomicTiles))
    #       N.b. minoverlap ensures that the MCWs computed within <Unsupervised.clustering.MultiSite.wrapper> are exact and that each GenomicTile is analyised exactly once
    #       WARNING: the above comment isn't true here. But beucase it's amplicons this shouldn't be a problem
  )
  # <-> #
  
  # # PREPARE INPUT #
  TFBS.clusters = list()
  TFBS.clusters$ClusterCoordinates = CurrentGenomicTiles
  TFBS.clusters$ClusterComposition = NULL
  
  bait.capture.parameters.list$CallContextMethylation[["coverage"]] = 20
  bait.capture.parameters.list$SortReadsByTFCluster[["coverage"]] = 20
  bait.capture.parameters.list[["unsupervised.clustering.coverage"]] = 20
  
  Unsupervised.clustering.MultiSite.wrapper(
    sampleSheet = sampleSheet,
    Samples = MySamples$DE,
    pool.samples = TRUE,
    genome = BSgenome.Mmusculus.KrebsEMBL.MutantsBatch1VB,
    MethylationCallingWindows = MethylationCallingWindows[unlist(indexes)],
    TFBS.clusters = TFBS.clusters,
    SNPs = NULL,
    TFBSs = TFBSs,
    CytosinesToMask = NULL,
    species.string.match = "_C",
    deduplicate = TRUE,
    params = bait.capture.parameters.list,
    cores = 8,
    tryCatchLog.write.error.dump.file = tryCatchLog.write.error.dump.file
  ) -> UnsupervisedClusteringResults
  
  return(UnsupervisedClusteringResults)
  
}

launch_rslurm_UnsupervisedClustering_PRA_mutants = function(working.directory, sampleSheet, Load.TFBSs.args, strong.binders.only, additional.jobname.info, NR.NODES = 64, tryCatchLog.write.error.dump.file){
  
  library(BSgenome.Mmusculus.KrebsEMBL.MutantsBatch1VB)
  
  setwd(working.directory)
  
  message("Loading TFBSs ... ")
  TFBSs = qs::qread("/g/krebs/barzaghi/analyses/08.02.24_RMCE_mutants/TFBSs.qs")
  GenomicTiles = qs::qread("/g/krebs/barzaghi/analyses/08.02.24_RMCE_mutants/MutantsBatch1VB_ranges.qs")
  seqlevels(GenomicTiles) = seqlevels(TFBSs)
  MethylationCallingWindows = SingleMoleculeFootprinting::Create_MethylationCallingWindows(
    TFBS_cluster_coordinates = GenomicTiles, 
    genomic.seqlenghts = GenomeInfoDb::seqlengths(BSgenome.Mmusculus.KrebsEMBL.MutantsBatch1VB), 
    fix.window.size = TRUE, 
    max.window.size = 50
  )
  
  l = length(MethylationCallingWindows)
  params.df = data.frame(NODE = seq(NR.NODES))
  data.frame(start.index = ceiling(seq(1,l, l/NR.NODES)), end.index = ceiling(seq(l/NR.NODES,l,l/NR.NODES))) -> df
  params.df$indexes = lapply(seq(NR.NODES), function(i){seq(df$start.index[i], df$end.index[i])})
  params.df$sampleSheet = sampleSheet
  params.df$Load.TFBSs.args = list(Load.TFBSs.args)
  params.df$strong.binders.only = strong.binders.only
  params.df$NODE = NULL
  params.df$tryCatchLog.write.error.dump.file = tryCatchLog.write.error.dump.file
  
  slurm.options <- list(time = '03:00:00', mem = "200G", qos = "normal", error = "slurm.%N.%j.err")
  
  jobname = paste0("UnsupervisedClustering_PRA_mutants_", Load.TFBSs.args$kind, ifelse(strong.binders.only, "_strong.binders.only", ""), additional.jobname.info)
  rslurm::slurm_apply(
    rslurm_UnsupervisedClustering_PRA_mutants, params.df, 
    nodes = NR.NODES, cpus_per_node = 8, jobname = jobname, rscript_path = Rscript_path,
    preschedule_cores = TRUE, slurm_options = slurm.options, libPaths = lib_paths, 
    submit = FALSE) -> slurm.object
  
  return(slurm.object)
  
}

retrieve_UnsupervisedClusteringResults = function(UnsupervisedClusteringResults_path, working.directory, Load.TFBSs.args, data.type, strong.binders.only, additional.jobname.info, NR.NODES = 200, NR.CORES = 8, tryCatchLog.write.error.dump.file){
  
  setwd(working.directory)
  source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
  
  list.of.files = list.files(UnsupervisedClusteringResults_path, "results_", full.names = TRUE)
  list.of.files.split = base::split(list.of.files, f = rep(seq(NR.NODES), each = length(list.of.files)/NR.NODES))
  names(list.of.files.split) = NULL
  
  data.frame(
    list.of.files = list.of.files,
    idx = rep(seq(NR.NODES), length(list.of.files)/NR.NODES)
  ) %>%
    group_by(idx) %>%
    summarise(list.of.files = list(list.of.files), .groups = "drop") %>%
    select(-idx) %>%
    mutate(
      partition.total.coverage = 0,
      single.site.coverage.thr = 20,
      data.type = data.type,
      cores = NR.CORES
    ) -> params.df
  
  jobname = paste0("retrieve_UnsupervisedClusteringResults_", Load.TFBSs.args$kind, ifelse(strong.binders.only, "_strong.binders.only", ""), additional.jobname.info)
  slurm.options <- list(time = '01:00:00', mem = "200G", qos = "normal", error = "slurm.%N.%j.err")#, partition = "bigmem")
  rslurm::slurm_apply(
    create.master.table, params.df, 
    nodes = NR.NODES, cpus_per_node = NR.CORES, jobname = jobname, rscript_path = Rscript_path,
    preschedule_cores = TRUE, slurm_options = slurm.options, libPaths = lib_paths, 
    submit = FALSE) -> slurm.object

  return(slurm.object)

}