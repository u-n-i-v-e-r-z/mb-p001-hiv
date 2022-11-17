# HEADER ====================================================== =
# Gen2i
# Heurteau Alexandre
# Tue Jan 25 14:14:42 2022
# R version : R version 4.0.5 (2021-03-31)
# System : x86_64, linux-gnu
# ============================================================= =

# DESC  ======================================================= =
# ============================================================= =

# SOURCE  ===================================================== =
library(data.table)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
source("libraries/hiv_lib.R")
# ============================================================= =

# SETSEED  ==================================================== =
set.seed(123)
# ============================================================= =

# SETWD  ==================================================== =
setwd("resources/mat_and_met/")
# ============================================================= =

# GLOBAL PARAMETERS =========================================== =
tad_res <- "50kb"
tad_res_kb <- as.numeric(gsub("kb","",tad_res))*1000
fill_by <- "rname"
# ============================================================= =


#  READ I/O  ================================================== =
# Tiled genome
bin_tot <- readRDS(paste0("data/others/genome_tiled_",tad_res,".rds"))
# FAI genome seqinfo
fai_lst <- readRDS("data/others/fai_lst.rds")
# Ranges features of interest
Go <- readRDS("data/f1/f1b/f1b.rds")
# TADS
tad_cd4a1_gr <- GRanges(fread(paste0("data/tads/Act_merged.bwa_mem.",tad_res,"_boundaries.bed"), select = c(1,2,3),col.names = c("seqnames","start","end")))
tad_cd4a1_gr$name <- "tad_cd4a1"
tad_cd4q1_gr <-  GRanges(fread(paste0("data/tads/NAct_merged.bwa_mem.",tad_res,"_boundaries.bed"), select = c(1,2,3),col.names = c("seqnames","start","end")))
tad_cd4q1_gr$name <- "tad_cd4q1"
tad_borders_gr <- g2i:::setSeqinfo(c(tad_cd4a1_gr, tad_cd4q1_gr), "hg19")
# ============================================================= =


#  RUN  ================================================== =
# IS
out_df <- NULL
# Create working GRanges for a given virus activity
for (tad_nm in unique(tad_borders_gr$name)) {
  borders_gr <- subset(tad_borders_gr, name == tad_nm)
  dtn <- distanceToNearest(Go, borders_gr)
  bin_dst_unsign <- ceiling(GenomicRanges::mcols(dtn)$distance/tad_res_kb)
  sign <- ifelse(start(Go[dtn@from]) < start(borders_gr[dtn@to]), -1, 1) 
  distances <- bin_dst_unsign * sign
  out_df <- rbind.data.frame(out_df,
                             cbind.data.frame(fill_by=GenomicRanges::mcols(Go[dtn@from])[[fill_by]],
                                              tad_type=tad_nm,
                                              prof=as.vector(distances)))
}


df <- as.data.frame(out_df %>% group_by(fill_by, prof,tad_type) %>% count() %>% ungroup(prof) %>% mutate(nn = n*100/sum(n))) 

gg <- list()
gg[["histogram.TADborders_distribution"]] <-  ggplot2::ggplot(df, aes(x=prof, y = nn, fill = fill_by)) + 
  ggplot2::geom_bar(stat = "identity", position = "dodge2") +
  ggplot2::xlim(-10,10) +
  ggplot2::facet_wrap(~  tad_type) + 
  ggplot2::xlab(paste0("Distance from TAD border in  (x",tad_res,")")) + 
  ggplot2::ylab(paste0("Percentage of IS at a given distance from TAD borders")) 


# Associated proportion test ----
for (tad_nm in c("Act","NAct")){
  
  borders_gr <- GRanges(fread(paste0("data/tads/",tad_nm,"_merged.bwa_mem.",tad_res,"_boundaries.bed"), select = c(1,2,3),col.names = c("seqnames","start","end")))
  tad_borders_gr <- setSeqinfo(borders_gr, "hg19", fai_lst = fai_lst)
  
  # GET DISTANCES
  all_fish_df <- NULL
  out_df <- NULL
  dtn <- distanceToNearest(Go, borders_gr)
  bin_dst_unsign <- ceiling(mcols(dtn)$distance/tad_res_kb)
  sign <- ifelse(start(Go[dtn@from]) < start(borders_gr[dtn@to]), -1, 1) 
  distances <- bin_dst_unsign * sign
  out_df <- rbind.data.frame(out_df,
                             cbind.data.frame(donor=Go[dtn@from]$rname,
                                              tad_type=tad_nm,
                                              in_tadb = ifelse(abs(as.vector(distances))>0 & abs(as.vector(distances)) < 5,"ADJACENT_4BINS",
                                                               ifelse(as.vector(distances) == 0, "TAD_BORD_BIN", "INSIDE_TAD")),
                                              prof=as.vector(distances)))
  
  
  prop_df <- NULL
  prop_l <- list()
  
  prop_l[["ADJACENT_4BINS"]] <- length(tad_borders_gr)*8/length(bin_tot)
  prop_l[["TAD_BORD_BIN"]] <- length(tad_borders_gr)/length(bin_tot)
  prop_l[["INSIDE_TAD"]] <- (length(bin_tot)-length(tad_borders_gr))/length(bin_tot)
  
  
  # Compute enrichment for each features at TAD borders
  # This is performed by a test comparing the proportion of IS at TADborders
  # versus the expected proportion (given by the number of total bin in the genome at a given resolution divided
  # by the total number of borders)
  for (bin_type in c("ADJACENT_4BINS","TAD_BORD_BIN", "INSIDE_TAD")){
    for (nm in unique(Go$rname)){
      nr <- nrow(subset(out_df, donor == nm &  in_tadb == bin_type))
      pt <- prop.test(nr, nrow(subset(out_df, donor == nm)), prop_l[[bin_type]], correct = F,alternative = "greater")
      pv <- pt$p.value
      ml10pv <- -log10(pv)
      # deal with inf cases
      ml10pv <- ifelse(is.finite(ml10pv),ml10pv,sign(ml10pv)*500)
      sign <- ifelse(pv < 0.05,ifelse(pv < 0.01,ifelse(pv < 0.001,"***","**"),"*"),"NS")
      prop_df <- rbind.data.frame(prop_df,
                                  cbind.data.frame(COL = nm,
                                                   ROW = bin_type,
                                                   number = nr,
                                                   ml10pv=ml10pv,
                                                   pv=pv,
                                                   tad_type=tad_nm,
                                                   sign=sign))
    }
  }
  gg[[paste0("Proportion_test.",ifelse(tad_nm=="Act","Active","Quiescent"),"_Cells")]] <- prop_plot_pdf(prop_df)
}







# ============================================================= =
