# HEADER ====================================================== =
# GEN2I
# Wed Nov  2 17:19:01 2022
# R version : R version 4.0.5 (2021-03-31)
# System : x86_64, linux-gnu
# ============================================================= =

# DESC  ======================================================= =
# Data : 
# - HiC MxM matrices
# - GenomicRranges of bait and anchor data
# - Distance min and max between couples 
# - Constraint (like a TAD structure)
#
# Description :
# Extract APA.
#
# ============================================================= =

# SOURCE  ===================================================== =
suppressWarnings(library(GenomicRanges))
suppressWarnings(library(rtracklayer))
suppressWarnings(library(dplyr))
suppressWarnings(library(InteractionSet))
suppressWarnings(library(data.table))
suppressWarnings(library("doParallel"))
suppressWarnings(library(ggplot2))
source("libraries/hiv_lib.R")
# ============================================================= =

# SETSEED  ==================================================== =
set.seed(123)
# ============================================================= =

# GLOBAL PARAMETERS =========================================== =
# ============================================================= =

#  READ I/O  ================================================== =
bait <- import("data/is/before_cART_RNAp_GO.bed")
anchor <- import("data/is/before_cART_RNAp_GO.bed")
bait.nm <- "before_cART_RNAp_GO"
anchor.nm <- "before_cART_RNAp_GO"

cm <- load_h5_matrix(list(Act = "data/contactMatrices/Act_merged.bwa_mem.50kb.obs_exp.h5"))
# ============================================================= =


# RUN ========================================================= =

# Prepare tile genome ----
`%ni%` <-  Negate('%in%')
gnm_sqfo <- rtracklayer::SeqinfoForBSGenome("hg19")
seqlevelsStyle(gnm_sqfo) <- "Ensembl"
gnm_sqfo <- gnm_sqfo[GenomeInfoDb::seqnames(gnm_sqfo)[!GenomeInfoDb::isCircular(gnm_sqfo)]]
mtdna_chr <- names(gnm_sqfo)[GenomeInfoDb::isCircular(gnm_sqfo)]
genome_gr <- GRanges(gnm_sqfo)%>%
  GenomeInfoDb::dropSeqlevels(mtdna_chr,pruning.mode = "coarse")

genome_t_gr <- tileGenome(seqlengths(genome_gr), tilewidth = 50e3, cut.last.tile.in.chrom = TRUE)

genome_t_gr$bin <- 1:length(genome_t_gr)
genome_t_gr$chr_bin <- unlist(lapply(split(genome_t_gr,seqnames(genome_t_gr)),
                                     function(l)1:length(l)))
genome_t_gr$name <- with(genome_t_gr,paste(seqnames,start,end,strand,chr_bin,sep="_"))

# APA° global parameters ----
win_apa <- 21L
res <- "50kb"
res_int <- as.numeric(gsub("(.+)kb","\\1",res))
dmin <- 21
dmax <- 80
xt_apa <- (win_apa-1)/2 # number of bin of the apa has to be odd to center the point
# Coordinate of central pixel
i_CP <- (win_apa-1)/2+1;j_CP <- (win_apa-1)/2+1 # Center win_apa interaction of size (1win_apa;1win_apa)
# Coordinate of 3x3 square centered on central pixel
m_size <- (win_apa-1)/2 # Middle pixel value - 1
i_CS <- m_size:(m_size+2);j_CS <- m_size:(m_size+2)
# Coordinate of Upper left 3x3 square from 3x3 square at center
i_ULS <- i_CS-3;j_ULS <- j_CS-3
# Coordinate of  Upper Right 3x3 square from 3x3 square at center (Compartment)
i_URS <- i_CS-3;j_URS <- j_CS+3
# Coordinate of  Bottom Right 3x3 square from 3x3 square at center (Compartment)
i_BRS <- i_CS+3;j_BRS <- j_CS+3
# Coordinate of  Bottom Right 3x3 square from 3x3 square at center (Compartment)
i_BLS <- i_CS+3;j_BLS <- j_CS-3


# Make overlaps with bait and anchor bed files
fol1 <- findOverlaps(genome_t_gr,bait)
fol2 <- findOverlaps(genome_t_gr,anchor)

bait_gr <- bait
anchor_gr <- anchor


bait_tile <- genome_t_gr[unique(fol1@from)]
anchor_tile <- genome_t_gr[unique(fol2@from)]

bait_tile$l_idx_a <- lapply(split(fol1, fol1@from), subjectHits)
anchor_tile$l_idx_b <- lapply(split(fol2, fol2@from), subjectHits)


pw_dt <- pairwiseInteractionDT(anc_go = anchor_tile,
                               bait_go =  bait_tile,
                               hic_res = 50e3,
                               apa_size =  15, 
                               d_min = dmin, 
                               d_max = dmax)


# Prepare coordinate to exctract apa submatrix (rows and cols)
pw_dt$cbw_a  <- sapply(1:nrow(pw_dt),function(idx){(as.numeric(pw_dt[idx,"anc_chrom_bin"])-xt_apa):(as.numeric(pw_dt[idx,"anc_chrom_bin"])+xt_apa)},simplify = F)
pw_dt$cbw_b  <- sapply(1:nrow(pw_dt),function(idx){(as.numeric(pw_dt[idx,"bait_chrom_bin"])-xt_apa):(as.numeric(pw_dt[idx,"bait_chrom_bin"])+xt_apa)},simplify = F)


# Loop to exctract apa
numCores <- detectCores()
registerDoParallel(numCores)  # use multicore, set to the number of our cores

pw_dt[anc_chrom_bin < bait_chrom_bin, uniq := paste0(anc_chrom, "_", anc_chrom_bin, "_", bait_chrom_bin)]
pw_dt[anc_chrom_bin > bait_chrom_bin, uniq := paste0(anc_chrom, "_", bait_chrom_bin, "_", anc_chrom_bin)]

apa_struct <- foreach::foreach(row=1:nrow(pw_dt), .export = "mp") %dopar% {
  anc_chrom_bin <- pw_dt[row,"anc_chrom_bin"]
  bait_chrom_bin <- pw_dt[row,"bait_chrom_bin"]
  b_a <- pw_dt[row,"anc_chrom_bin"]
  b_b <- pw_dt[row,"bait_chrom_bin"]
  chr_a <- pw_dt[row,"anc_chrom"] 
  chr_b <- pw_dt[row,"bait_chrom"]
  cbw_a <- unlist(pw_dt[row,"cbw_a"])
  cbw_b <- unlist(pw_dt[row,"cbw_b"])
  
  if(anc_chrom_bin < bait_chrom_bin){
    rev <- T
  }else{
    rev <- F
  }
  # Compute the APA matrix and return it in myMAT as a list of matrix
  cm_mat <- cm[["Act"]][[as.character(unlist(chr_a))]][cbw_b+20,cbw_a+20] # +20 comes from NA extend performed on all contact matrices to capture interacting couples at chromosomes borders
  cm_mat[!is.finite(cm_mat) | cm_mat==0] <- NA
  cmq_mat <- matrix(ntile(as.matrix(cm_mat),21*21),dim(cm_mat)[1],dim(cm_mat)[2])
  
  
  if(rev){
    cm_mat <- matrix(rev(cm_mat),21,21)
    cmq_mat <- matrix(rev(cmq_mat),21,21)
  }
  
  l_ <- list()
  l_[["cp_r"]] <- cm_mat[i_CP,j_CP]
  l_[["cs_r"]] <- cm_mat[i_CS,j_CS]
  l_[["uls_r"]] <- cm_mat[i_ULS,j_ULS]
  l_[["urs_r"]] <- cm_mat[i_URS,j_URS]
  l_[["bls_r"]] <- cm_mat[i_BLS,j_BLS]
  l_[["brs_r"]] <- cm_mat[i_BRS,j_BRS]
  
  
  raw_df <- data.frame(chr_a=chr_a,chr_b=chr_b,b_a=b_a,b_b=b_b,anc_chrom_bin=anc_chrom_bin,bait_chrom_bin=bait_chrom_bin)
  for(col in names(l_)){
    raw_df[,col] <- mean(as.array(l_[[col]]),na.rm=T)
  }
  
  l_ <- list()
  l_[["cp_r"]] <- cmq_mat[i_CP,j_CP]
  l_[["cs_r"]] <- cmq_mat[i_CS,j_CS]
  l_[["uls_r"]] <- cmq_mat[i_ULS,j_ULS]
  l_[["urs_r"]] <- cmq_mat[i_URS,j_URS]
  l_[["bls_r"]] <- cmq_mat[i_BLS,j_BLS]
  l_[["brs_r"]] <- cmq_mat[i_BRS,j_BRS]
  quant_df <-data.frame(chr_a=chr_a,chr_b=chr_b,b_a=b_a,b_b=b_b,anc_chrom_bin=anc_chrom_bin,bait_chrom_bin=bait_chrom_bin)
  for(col in names(l_)){
    quant_df[,col] <- mean(as.array(l_[[col]]),na.rm=T)
  }
  
  bin_df <- data.frame(chr_a=chr_a,chr_b=chr_b,b_a=b_a,b_b=b_b,anc_chrom_bin=anc_chrom_bin,bait_chrom_bin=bait_chrom_bin,rev=rev)
  
  list(bin_df=bin_df,
       scm_m=as.matrix(cm_mat),scmq_m=cmq_mat,
       raw_df=raw_df,
       quant_df=quant_df)
}
stopImplicitCluster()



f_apa_struct <- list(bin_df=do.call("rbind",lapply(apa_struct,"[[","bin_df")),
                     raw_df=do.call("rbind",lapply(apa_struct,"[[","raw_df")),
                     quant_df=do.call("rbind",lapply(apa_struct,"[[","quant_df")),
                     scm_lm=lapply(apa_struct,"[[","scm_m"),
                     scmq_lm=lapply(apa_struct,"[[","scmq_m"),
                     a_gr=bait_gr,
                     b_gr=anchor_gr,
                     t_gr=genome_t_gr)


apa_mat_p <- list()
apa_mat_p[["MEAN"]] <- apply(simplify2array(f_apa_struct$scm_lm), 1:2, mean, na.rm = TRUE)
apa_mat_p[["QUANTMEAN"]] <- apply(simplify2array(f_apa_struct$scmq_lm), 1:2, mean,na.rm=T)


for (type in c("QUANTMEAN","MEAN")){
  
  m2t <- matrix2tibble(apa_mat_p[[type]], 21, 21)
  names(m2t) <- c(paste0(bait.nm,".Bait") , paste0(anchor.nm,".Anchor"), "values")  # Bait  Anchor counts
  min <- min(min(m2t$values), 0)
  max <- max(max(m2t$values), 100)
  
  
  
  gg <- ggplot2::ggplot(m2t, ggplot2::aes_string(x=paste0(anchor.nm,".Anchor"), y=paste0(bait.nm,".Bait"), fill="values")) +
    ggplot2::geom_tile() +  ggtitle(paste0("Condition : Act\n ",
                                           "HiC metrics : OE/KR")) +
    ggplot2::theme(axis.text.y = element_blank(),axis.text.x = element_blank()) +
    scale_y_discrete(breaks=(21+1)/2) +
    scale_x_discrete(breaks=(21+1)/2)
  
  ggs <- gg + ggplot2::scale_fill_gradientn(colors = scico::scico(6, palette = "lajolla"), breaks = round(seq(min, max, length.out = 5),digits = 3),
                                            limits = c(min, max))
  ggns <- gg + ggplot2::scale_fill_gradientn(colors = scico::scico(100, palette = "lajolla"))
  
}



# ============================================================= =