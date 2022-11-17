# HEADER ====================================================== =
# GEN2I
# Wed Oct 26 10:29:57 2022
# R version : R version 4.0.5 (2021-03-31)
# System : x86_64, linux-gnu
# ============================================================= =

# DESC  ======================================================= =
# Description :
# Libraries for HIV paper
# ============================================================= =

# RUN ========================================================= =

# Graphics related functions ----
fisher.grad.grad <- function(df, grad_row, grad_col, verbose=F) {
  r_df <- NULL
  for(i in unique(df[[grad_row]])){
    for(j in unique(df[[grad_col]])){
      cr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(grad_row)==i & !!sym(grad_col)==j)))
      cnr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(grad_row)!=i & !!sym(grad_col)==j)))
      ncr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(grad_row)==i & !!sym(grad_col)!=j)))
      ncnr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(grad_row)!=i & !!sym(grad_col)!=j)))
      f_mat <- matrix(c(L_C=cr,L_nC=cnr,C_nL=ncr,nC_nL=ncnr),nrow=2,ncol=2)
      
      rownames(f_mat) <- c(i,paste0("no_",i))
      colnames(f_mat) <- c(j,paste0("no_",j))
      mosaicplot(f_mat)
      
      if(verbose) {
        print(i)
        print(j)
        print(f_mat)
      }
      ft <- fisher.test(f_mat)
      mypv <- formatC(as.numeric(ft$p.value, format = "e", digits = 2))
      lfc <- round(log2(ft$estimate),2)
      lfc <- ifelse(is.finite(lfc),lfc,sign(lfc)*10)
      sign <- ifelse(ft$p.value < 0.05,ifelse(ft$p.value < 0.01,ifelse(ft$p.value < 0.001,"***","**"),"*"),"NS")
      r_df <- rbind.data.frame(r_df,
                               cbind.data.frame(ROW=i,COL=j,
                                                pv=mypv,lfc=lfc,sign=sign))
    }
  }
  
  r_df$ROW <- factor(r_df$ROW, levels = sort(unique(r_df$ROW)))
  r_df$COL <- factor(r_df$COL, levels = rev(sort(unique(r_df$COL))))
  r_df
}

fisher_plot_asym_pdf <- function(DF, main="", norm=NULL, range=c(-2,2)){
  if(!is.null(norm)) DF <- DF %>% mutate(lfc=rangeMinMaxAsym(x = lfc,nmin = -2,nmax=0,pmin = 0,pmax = 2))
  colfunc <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))
  p <- ggplot(DF, aes(x = ROW, y = COL, fill=lfc)) +      geom_tile(col=NA) +
    theme_minimal()+
    coord_equal(ratio=1) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top",expand=c(-1,0)) +
    theme(axis.text.x=element_text(angle = 90,size=15, vjust=0, hjust=0.6,face = "bold",
                                   margin=margin(0,0,0,0)),
          axis.text.y=element_text(size=15, margin=margin(0,0,0,0),face = "bold"),
          legend.background = element_rect(color="grey",size=0.8, linetype="dashed"),
          plot.margin = unit(c(0,1,0,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    
    geom_text(aes(label =sign),size=4) +
    ggtitle(main) +
    scale_fill_gradientn(colours = colfunc(10) ,
                         guide = guide_colorbar(barwidth = 0.8,
                                                title = 'Log2 \nOddsRatio',
                                                label.theme = element_text(size=10, vjust=1, hjust=.5,face = "bold",angle = 0),
                                                barheight = 10,
                                                nbin = 20,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                ticks = FALSE),
                         breaks=c(min(range),0,max(range)),
                         labels=c(min(range),0,max(range)),
                         limits=c(min(range),max(range)))
  
  return(p)
}

prop_plot_pdf <- function(DF, main=""){
  range <- c(round(min(DF$ml10pv),1),round(max(DF$ml10pv),1))
  colfunc <- colorRampPalette(rev(brewer.pal(n = 9, name = "PuBuGn")))
  p <- ggplot(DF, aes(x = ROW, y = COL, fill=ml10pv)) +      geom_tile(col=NA) +
    theme_minimal()+
    coord_equal(ratio=1) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top",expand=c(-1,0)) +
    theme(axis.text.x=element_text(angle = 90,size=15, vjust=0, hjust=0.6,face = "bold",
                                   margin=margin(0,0,0,0)),
          axis.text.y=element_text(size=15, margin=margin(0,0,0,0),face = "bold"),
          plot.margin = unit(c(0,1,0,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    
    geom_text(aes(label = paste0(sign,"\nIn = ",number)),size=4) +
    ggtitle(main) +
    scale_fill_gradientn(colours = rev(colfunc(10)),
                         guide = guide_colorbar(barwidth = 0.8,
                                                title = '-log10(pvalue)',
                                                label.theme = element_text(size=10, vjust=1, hjust=.5,face = "bold",angle = 0)))
  
  return(p)
}



# Data transformation related functions ----
rangeMinMaxAsym <- function(x,nmin=-1,nmax=0,pmin=0,pmax=1)
{
  neg <- x[x<=0]
  pos <- x[x>=0]
  if(length(neg)){
    if(length(unique(neg))-1)
      ret_neg <- nmin+((neg- min(c(-0.0001,neg),na.rm=T))*(nmax-nmin)) /(max(c(-0.0001,neg),na.rm=T)-min(c(-0.0001,neg),na.rm=T))
    else
      ret_neg <- rep(nmin,length(neg))
  }else{
    ret_neg <- NULL
  }
  if(length(pos)){
    if(length(unique(pos))-1)
      ret_pos <- pmin+((pos- min(c(0.0001,pos),na.rm=T))*(pmax-pmin)) /(max(c(0.0001,pos),na.rm=T)-min(c(0.0001,pos),na.rm=T))
    else
      ret_pos <- rep(pmax,length(pos))
  }else{
    ret_pos <- NULL
  }
  idx_neg <- x<0
  idx_pos <- x>=0
  x[idx_neg] <- ret_neg
  x[idx_pos] <- ret_pos
  x
}# TODO example : rangeMinMax(-5:5,-1,1)



# Ranges related functions ----
setSeqinfo <- function(gr,genome,seqlevelstyle="Ensembl",fai_lst){
  gnm_sqfo <- fai_lst[[seqlevelstyle]][[genome]]
  if (is.null(gnm_sqfo)) fireError(paste0("The genome ",genome," is not available.","Contact developpers for them to add it to the database."))
  GenomeInfoDb::seqlevelsStyle(gnm_sqfo) <- seqlevelstyle
  GenomeInfoDb::seqlevelsStyle(gr) <- seqlevelstyle
  GenomeInfoDb::seqlevels(gr, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(gnm_sqfo)
  GenomeInfoDb::seqinfo(gr,pruning.mode=c("coarse")) <- gnm_sqfo
  GenomeInfoDb::seqlevels(gr) <-  GenomeInfoDb::seqlevelsInUse(gr)
  gr
}



# HIC related functions ----
# Turn TSV To InteractionSet object
he.df2gi <- function(df){
  row.regions <- GenomicRanges::GRanges(df$V1, IRanges(df$V2,df$V3))# interaction start
  col.regions <- GenomicRanges::GRanges(df$V4, IRanges(df$V5,df$V6))# interaction end
  gi <- InteractionSet::GInteractions(row.regions, col.regions)
  gi$norm.freq <- df$V7 # Interaction frequencies
  return(gi)
}


#' Create pairwise bin/bin interactions data frame
#' to use to dump HiC contacts
#'
#' @param anc_go Tile genome of anchor
#' @param bait_go Tile genome of Bait
#' @param apa_size Size of APA. Default 41
#' @param constraint_go Constraint GO. TAD for instance or reduce(tile_genome) for all interactions
#' @param hic_res
#' @param d_min distance minimum between bait & anc in bin
#' @param d_max distance max between bait & anc in bin
#'
#' @return
#' @export
pairwiseInteractionDT <- function(anc_go=NULL,
                                  bait_go=NULL,
                                  apa_size=41,
                                  d_min=4,
                                  d_max=100,
                                  hic_res=50e3){
  # chromosomes collection
  chrs <- unique(c(as.vector(seqnames(bait_go)),as.vector(seqnames(anc_go))))
  Go_cons <- NULL
  for (chr in chrs) {
    min_anc <- min(start(subset(anc_go, seqnames == chr)))
    max_anc <- max(end(subset(anc_go, seqnames == chr)))
    min_bait <- min(start(subset(bait_go, seqnames == chr)))
    max_bait <- max(end(subset(bait_go, seqnames == chr)))
    Go_cons <- rbind(Go_cons, data.frame(start = min(min_anc,min_bait), end = max(max_anc,max_bait) ,seqnames = chr))
  }
  constraint_go <- GenomicRanges::GRanges(Go_cons)
  names(constraint_go) <- NULL
  #Bait feature
  bait_r_go <- IRanges::resize(bait_go, 1, fix="center") #resize the feature to do the overlap on domains without multi hits
  bait_r_go$cb_b <- ceiling((GenomicRanges::start(bait_r_go)) / hic_res )
  bait_r_go$l_idx_b <- 1:length(bait_r_go)
  mol_dt <- data.table::data.table( data.frame(IRanges::mergeByOverlaps(constraint_go, bait_r_go) ))
  mol_dt$filt <- paste0(mol_dt$constraint_go.seqnames, "_", mol_dt$constraint_go.start)
  mol_b_dt <- mol_dt[ , c( "bait_r_go.name", "chr_bin", "filt", "bait_r_go.seqnames", "bait_r_go.strand" ) ]
  colnames(mol_b_dt) <- c("bait_name", "bait_chrom_bin", "constraint_id", "bait_chrom", "bait_strand")
  data.table::setkey(mol_b_dt, constraint_id)
  #Anchoring feature
  anc_r_go <- IRanges::resize(anc_go, 1, fix="center") #resize the feature to do the overlap on domains without multi hits
  anc_r_go$cb_b <- ceiling((GenomicRanges::start(anc_r_go)) / hic_res )
  anc_r_go$l_idx_b <- 1:length(anc_r_go)
  mol_dt <- data.table::data.table( data.frame(IRanges::mergeByOverlaps(constraint_go, anc_r_go) ))
  mol_dt$filt <- paste0(mol_dt$constraint_go.seqnames, "_", mol_dt$constraint_go.start)
  mol_a_dt <- mol_dt[ , c( "anc_r_go.name", "chr_bin", "filt", "anc_r_go.seqnames", "anc_r_go.strand" ) ]
  colnames(mol_a_dt) <- c("anc_name", "anc_chrom_bin", "constraint_id", "anc_chrom", "anc_strand")
  data.table::setkey(mol_a_dt, constraint_id)
  # Merge
  inter_dt <- mol_b_dt[mol_a_dt, allow.cartesian=TRUE,nomatch=0]
  
  if(nrow(inter_dt) == 0) return(NULL)
  
  # Mark the one to revert
  inter_dt$rev <- ifelse(inter_dt$anc_chrom_bin < inter_dt$bait_chrom_bin, T, F)
  # Add distance between int and anc in bin
  inter_dt$dst_bin <- abs(inter_dt$bait_chrom_bin - inter_dt$anc_chrom_bin)
  # Filter distance min and max
  inter_dt <- subset(inter_dt, dst_bin > d_min & dst_bin < d_max)
  # RETURN DATA.TABLE
  inter_dt
}


load_h5_matrix <- function(h5path){
  l_cm <- list()
  # Load H5 matrix
  if (!is.null(h5path)) {
    for(h5.nm in names(h5path)){
      # get h5 name info by h5path
      cm_info <- rhdf5::h5ls(h5path[[h5.nm]])
      l_cm[[h5.nm]] <- sapply(cm_info$name, function(nm) { HDF5Array::HDF5Array(h5path[[h5.nm]], nm)})
    }
  } else {
    # TODO : put loop here to load cm in l_cm
    mat_nm <- paste0(cond,"_",rep,"_",norm,"_",as.integer(res))
    filt_cm <- filterCmDB(cmname=list(in_=list(value=mat_nm,strict=T)))
    l_cm[["to_dev"]] <- dumpCM(filt_cm)
  }
  l_cm
}

xtna.mat <- function(mat=NULL,win=20){
  na_mat <- matrix(NA,dim(mat)[1]+2*win,dim(mat)[1]+2*win)
  na_mat[(win+1):((win)+dim(mat)[1]),(win+1):((win)+dim(mat)[1])] <- mat
  return(na_mat)
}

matrix2tibble <- function(m, n_r=41, n_c = 41) {
  colnames(m) <- paste0(seq(1,n_c))
  rownames(m) <- paste0(seq(1,n_r))
  m %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Bait") %>%
    tidyr::pivot_longer(-c(Bait), names_to = "Anchor", values_to = "counts") %>%
    dplyr::mutate(Anchor= forcats::fct_relevel(Anchor,colnames(m))) %>%
    dplyr::mutate(Bait= forcats::fct_relevel(Bait,rev(rownames(m))))
}

# ============================================================= =