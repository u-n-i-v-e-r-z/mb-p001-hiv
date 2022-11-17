# HEADER ====================================================== =
# GEN2I
# Wed Oct 26 11:33:12 2022
# R version : R version 4.0.5 (2021-03-31)
# System : x86_64, linux-gnu
# ============================================================= =

# DESC  ======================================================= =
# Orientation of Integration Sites based on nearest gene strand.
# ============================================================= =


# SOURCE  ===================================================== =
library(rtracklayer)
library(GenomicRanges)
source("libraries/hiv_lib.R")
# ============================================================= =

# SETSEED  ==================================================== =
set.seed(123)
# ============================================================= =


#  READ I/O  ================================================== =
genes_gr <- setSeqinfo(rtracklayer::import("data/genes/Homo_sapiens.GRCh37.85.chr.genes.bed"),"hg19", fai_lst = fai_lst)
is <- rtracklayer::import("path/to/is_of_interest.bed")
# ============================================================= =


# RUN ========================================================= =
is_go <- rtracklayer::import(is)
GenomicRanges::strand(is_go) <- "*"
nrst <- GenomicRanges::nearest(is_go, genes_gr)
GenomicRanges::strand(is_go) <- GenomicRanges::strand(genes_gr[nrst])
rtracklayer::export.bed(is_go,paste0(licht_go_repo,"/",out_nm))
# ============================================================= =