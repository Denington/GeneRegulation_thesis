#redoing 5.2.2 but as a log2(obs/exp) heatmap


ctcf_bed <- import.bedGraph('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/1_CTCFCD69negDPWTR1/1_CTCFCD69negDPWTR1_peaks.subpeaks.bedgraph')
rad21_bed <- import.bedGraph('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/Rad21_ChIPseq_in_CD69negDP/WT_CD69negDP/Rad21CD69negDPWTR1_mapped_sorted_RemoveDuplicates_peaks.subpeaks.bed')

ctcfOnly <- ctcf_bed[-queryHits(findOverlaps(ctcf_bed,rad21_bed,maxgap = 500))]
rad21Only <- rad21_bed[-subjectHits(findOverlaps(ctcf_bed,rad21_bed,maxgap = 500))]
both.marks <- ctcf_bed[queryHits(findOverlaps(ctcf_bed,rad21_bed,maxgap = 500))]


g <- conds.gr.list$dev
length(unique(queryHits(findOverlaps(g,ctcf_bed,maxgap = 250))))
length(unique(queryHits(findOverlaps(g,rad21_bed,maxgap = 250))))

getbars <- function(gd,direction,maxgap=500){
  ctcfOnlyL <- length(unique(queryHits(findOverlaps(gd,ctcfOnly,maxgap = maxgap))))
  rad21OnlyL <- length(unique(queryHits(findOverlaps(gd,rad21Only,maxgap = maxgap))))
  bothL <- length(unique(queryHits(findOverlaps(gd,both.marks,maxgap = maxgap))))
  neitherL <- length(gd) - ctcfOnlyL - rad21OnlyL - bothL
  
  neitherL ; ctcfOnlyL ; rad21OnlyL ; bothL
  
  z=data.frame(x=c('Neither','RAD21 only','CTCF only','both'),y=c(neitherL,rad21OnlyL,ctcfOnlyL,bothL),d=c(direction,direction,direction,direction))
  return(z)
}




test <- permTest(A=g[g$log2FoldChange < 0 & g$padj < 0.05,],B=ctcfOnly,ntimes = 100,
                      randomize.function = randomizeRegions,evaluate.function = numOverlaps,genome='mm9')
observed = test$numOverlaps$observed
expected = mean(test$numOverlaps$permuted)

permTest(A=g[g$log2FoldChange < 0 & g$padj < 0.05,],B=ctcfOnly,ntimes = 100,randomize.function = resampleRegions,universe = g,
          evaluate.function = numOverlaps,genome='mm9',count.once=T)

permTest(A=g[g$log2FoldChange > 0 & g$padj < 0.05,],B=rad21Only,ntimes = 100,randomize.function = resampleRegions,universe = g,
         evaluate.function = numOverlaps,genome='mm9',count.once=T)

permTest(A=g[g$log2FoldChange > 0 & g$padj < 0.05,],B=both.marks,ntimes = 100,randomize.function = resampleRegions,universe = g,
         evaluate.function = numOverlaps,genome='mm9',count.once=T)

get_obsexp <- function(genes,mark){
  up <- genes[genes$padj < 0.05 & genes$log2FoldChange > 0]
  down <- genes[genes$padj < 0.05 & genes$log2FoldChange < 0]
  
  upPT <- permTest(A=up,B=mark,ntimes = 500,randomize.function = resampleRegions,universe = genes,
                   evaluate.function = numOverlaps,genome='mm9',count.once=T)
  downPT <- permTest(A=down,B=mark,ntimes = 500,randomize.function = resampleRegions,universe = genes,
                   evaluate.function = numOverlaps,genome='mm9',count.once=T)
  
  obsUp = upPT$numOverlaps$observed
  expUp = mean(upPT$numOverlaps$permuted)
  
  obsDown = downPT$numOverlaps$observed
  expDown = mean(downPT$numOverlaps$permuted)
  
  return(list(obsUp/expUp,obsDown/expDown))
}

dev_ctcf <- get_obsexp(conds.gr.list$dev,ctcfOnly)
dev_rad21 <- get_obsexp(conds.gr.list$dev,rad21Only)
dev_both <- get_obsexp(conds.gr.list$dev,both.marks)

dko_ctcf <- get_obsexp(conds.gr.list$dko,ctcfOnly)
dko_rad21 <- get_obsexp(conds.gr.list$dko,rad21Only)
dko_both <- get_obsexp(conds.gr.list$dko,both.marks)

ctcfko_ctcf <- get_obsexp(conds.gr.list$ctcf_ko,ctcfOnly)
ctcfko_rad21 <- get_obsexp(conds.gr.list$ctcf_ko,rad21Only)
ctcfko_both <- get_obsexp(conds.gr.list$ctcf_ko,both.marks)

rad21ko_ctcf <- get_obsexp(conds.gr.list$rad21_ko,ctcfOnly)
rad21ko_rad21 <- get_obsexp(conds.gr.list$rad21_ko,rad21Only)
rad21ko_both <- get_obsexp(conds.gr.list$rad21_ko,both.marks)


ctcfonlyOE <- log2(t(matrix(c(dev_ctcf[[1]],dko_ctcf[[1]],ctcfko_ctcf[[1]],rad21ko_ctcf[[1]],
       dev_ctcf[[2]],dko_ctcf[[2]],ctcfko_ctcf[[2]],rad21ko_ctcf[[2]]),ncol = 2 )))

rad21onlyOE <-log2(t(matrix(c(dev_rad21[[1]],dko_rad21[[1]],ctcfko_rad21[[1]],rad21ko_rad21[[1]],
       dev_rad21[[2]],dko_rad21[[2]],ctcfko_rad21[[2]],rad21ko_rad21[[2]]),ncol = 2 )))

bothOE <-log2(t(matrix(c(dev_both[[1]],dko_both[[1]],ctcfko_both[[1]],rad21ko_both[[1]],
           dev_both[[2]],dko_both[[2]],ctcfko_both[[2]],rad21ko_both[[2]]),ncol = 2 )))

allOE <- rbind(ctcfonlyOE,rad21onlyOE,bothOE)
rownames(allOE) <- c('CTCF only','CTCF only','RAD21 only','Rad21 only','Both','Both')
colnames(allOE) <- c('Developmental','DKO','CTCFKO','Both')

Heatmap(allOE,cluster_rows = F,cluster_columns = F,border = T,column_gap = unit(0, "mm"),rect_gp = gpar(col = "black", lwd = 1),
        row_names_side = 'left',column_names_side = 'bottom')

#would it make more sense to have it as dev down /up

ctcfOE <- c(ctcfko_ctcf[[1]], rad21ko_ctcf[[1]],dko_ctcf[[1]],dev_ctcf[[1]],  ctcfko_ctcf[[2]], rad21ko_ctcf[[2]],dko_ctcf[[2]],dev_ctcf[[2]])  
rad21OE <- c(ctcfko_rad21[[1]], rad21ko_rad21[[1]],dko_rad21[[1]],dev_ctcf[[1]],  ctcfko_rad21[[2]], rad21ko_rad21[[2]],dko_rad21[[2]],dev_rad21[[2]])  
bothOE <- c(ctcfko_both[[1]], rad21ko_both[[1]],dko_both[[1]],dev_ctcf[[1]],  ctcfko_both[[2]], rad21ko_both[[2]],dko_both[[2]],dev_both[[2]])  

x <- data.frame('CTCF only'=ctcfOE,'RAD21 only'=rad21OE,'Both'=bothOE)
rownames(x) <- c('ΔCTCF up','ΔRAD21 up','ΔCTCF/ΔRAD21 up','SP up','ΔCTCF down','ΔRAD21 down','ΔCTCF/ΔRAD21 down','SP down')
allOE <- log2(x)
Heatmap(allOE,cluster_rows = F,cluster_columns = F,border = T,column_gap = unit(0, "mm"),rect_gp = gpar(col = "black", lwd = 1),
        row_names_side = 'left',column_names_side = 'top',heatmap_legend_param = list(title = "log2(o/e)", at = seq(-0.5, 0.5)))
