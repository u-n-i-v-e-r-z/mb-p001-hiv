# HEADER ====================================================== =
# GEN2I
# Wed Nov  2 17:04:29 2022
# R version : R version 4.0.5 (2021-03-31)
# System : x86_64, linux-gnu
# ============================================================= =

# DESC  ======================================================= =
# Data : 
# - HIC explorer TSV data
#
# Description :
# Turn TSV adjacencyy matrix to MxM contact matrices
#
# ============================================================= =

# SETSEED  ==================================================== =
set.seed(123)
# ============================================================= =

# SOURCE  =======================================================
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(InteractionSet))
# ===============================================================


# SETSEED  ======================================================
set.seed(123)
options(scipen = 100, digits = 4)
# ===============================================================

# GLOBAL PARAMETERS =============================================
xtna.mat <- function(mat=NULL,win=20){
  na_mat <- matrix(NA,dim(mat)[1]+2*win,dim(mat)[1]+2*win)
  na_mat[(win+1):((win)+dim(mat)[1]),(win+1):((win)+dim(mat)[1])] <- mat
  return(na_mat)
}
# ===============================================================

#  READ I/O  ====================================================

fai_lst <- readRDS("data/others/fai_lst.rds")

# NOTRUN
tsv <- data.table::fread("path/to/HiC_OE_KR.50kb.tsv") # Input from HiCExplorer adjacency matrices
res <- 50000
genome <- "hg19"
seqlevelstyle <- "Ensembl"
mtdna <- F

# ===============================================================


# RUN ===========================================================
# SET INPUT
# ************************************************************
tsv <- data.table::fread(tsv)

# Create genome ranges based on genome.FAI (stored in g2i lib)
# ************************************************************
gnm_sqfo <- fai_lst[[seqlevelstyle]][[genome]]
mtdna_chr <- names(gnm_sqfo)[GenomeInfoDb::isCircular(gnm_sqfo)]
if ( !mtdna ) gnm_sqfo <- gnm_sqfo[GenomeInfoDb::seqnames(gnm_sqfo)[!GenomeInfoDb::isCircular(gnm_sqfo)]]
genome_gr <- GenomicRanges::GRanges(gnm_sqfo)

# Load tiled genome at given resolut
# ***********************************
genome_t_gr <- g2i:::loadranges(opt$bed, genome = genome, "Ensembl")

# Turn TSV into InteractionSet object
# ***********************************
tsv[, V2 := V2+1] # Add an offset as it is 0-based
tsv[, V5 := V5+1] # Add an offset as it is 0-based
tsv <- if (!mtdna) tsv[V1!=mtdna_chr & V4!=mtdna_chr,] else tsv # Filter mtdna
h_gi <- he.df2gi(tsv)
ensembldb::seqlevelsStyle(h_gi) <- "Ensembl"

# Replace regions with full genome (HicExplorer does not output bin where there is a missing value)
# *************************************************************************************************
GenomeInfoDb::seqinfo(h_gi) <- gnm_sqfo[seqlevels(h_gi)]
h_gi <- trim(h_gi)
h_intra_gi <- h_gi[InteractionSet::intrachr(h_gi)]
InteractionSet::replaceRegions(h_intra_gi) <- GenomicRanges::GRanges(genome_t_gr)
h_intra_r_gi <- InteractionSet::regions(h_intra_gi)

# Built contact matrix from TSV file
# **********************************
cm_l <- list()
s_chr <- as.character(unique(seqnames(h_intra_r_gi)))

invisible(sapply(s_chr,
                 function(chr){
                   h_intra_r_chr_gi <- h_intra_r_gi[seqnames(h_intra_r_gi)==chr]
                   h_intra_r_chr_cm <- InteractionSet::inflate(h_intra_gi,
                                                               h_intra_r_chr_gi,h_intra_r_chr_gi,
                                                               fill=h_intra_gi$norm.freq,swap=T,sparse=T)
                   
                   cm_l[[chr]] <<- xtna.mat(as.matrix(as.matrix(h_intra_r_chr_cm)))
                   getKDiag(cm_l[[chr]],-1:1) <<- NA
                   if (dim(cm_l[[chr]])[1] != (length(subset(genome_t_gr, seqnames==chr))))
                     message("Dimension of your matrix does not correspond to tiled genome. Check your tiles and matrix")
                 }
))

HDF5Array::setHDF5DumpFile(out_h5)
for(chr in names(cm_l)){
  cm <- cm_l[[chr]]
  HDF5Array::setHDF5DumpName(chr)
  cm_c <-  as.matrix(cm_l[[chr]])
  invisible(as(cm_c, "HDF5Matrix"))
}

# ===============================================================

