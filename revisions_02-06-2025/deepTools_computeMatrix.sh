#!/bin/bash
#BATCH -A krebs                # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 16                       # number of cores
#SBATCH --mem 128G                  # memory pool for all cores
#SBATCH -t 2-00:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o /g/krebs/barzaghi/bash_scripts/slurm_err_out/slurm.%N.%j.out          # STDOUT
#SBATCH -e /g/krebs/barzaghi/bash_scripts/slurm_err_out/slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=guido.barzaghi@embl.de # send-to address

module load deepTools/3.5.5-foss-2022b
cd /g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025

computeMatrix scale-regions -S /g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025/SRR1973530_71ac45d8de4f.bw \
  /g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K27Ac/NT/PMID_24905168/signal/H3K27Ac_mESC_NT_SRX499119_rpgc.bigwig \
	/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K4me1/NT/PMID_22170606/signal/H3K4me1_mESC_NT_SRX080175_rpgc.bigwig \
	/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/mESC/H3K4me3/NT/PMID_28212747/signal/H3K4me3_mESC_NT_SRX1342333_rpgc.bigwig \
	/g/krebs/DWH/public/MusMusculus/sequencing/ChIP-seq/ESC/CTCF/WT/UniBind_mESC_ChIPseq/signal/ESC_CTCF_SRX160828_rpgc.bigwig \
	-R enhancer.bed promoter.bed intergenic.bed Ctcf.bed \
	--skipZeros -o deepTools_matrix.gz -m 5000 -bs 200 -p 16
