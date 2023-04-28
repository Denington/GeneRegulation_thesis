#plotting histone marks/rad21 around genes of different conditions using ggplot
#not bothering with ggplot. if want nice facet structure, make a false plot with same dimensions
#can't false plot as can't save plot genomation :(
#still feel plots have the 2 mid points mixed up (more vis with heatMatrices)- fixed
#looks much more distinctive if dif genes have a logfold change cutoff


library(ggplot2)
library(genomation)
library(ggpubr)
library(GenomicRanges)
library(ggsci)
library(rtracklayer)
library(Cairo)
library(RColorBrewer)

theme_set(theme_classic(base_size = 6))

load('~/Documents/Projects/Thymocyte_HiC/Thesis/scripts/mm9_gene_expression.rdata.RData')
#take wtwt,wtko etc from insulation_differences_gene_expression

CairoPDF('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/ctcf_rad21.at.genes.pdf')

difexp_to_granges <- function(genechanges){
  changes <- genechanges #[wtwt$padj < 0.05,]
  changes.gr <- mm9_pc_genes.gr[mm9_pc_genes.gr$ens_id %in% genechanges$ensembl_gene_id ]
  changes.gr <- changes.gr[!duplicated(changes.gr$ens_id)]
  changes.gr <- promoters(changes.gr,upstream = 100,downstream = 100 ) #YES!
  
  changes <- changes[changes$ensembl_gene_id %in% changes.gr$ens_id,]
  changes <- changes[!duplicated(changes$ensembl_gene_id),]
  
  changes <- changes[order(changes$ensembl_gene_id),]
  changes.gr <- changes.gr[order(changes.gr$ens_id)]
  changes.gr <- promoters(changes.gr,upstream = 100,downstream = 100)
  
  changes.gr$log2FoldChange <- changes$log2FoldChange
  changes.gr$padj  <- changes$padj
  return(changes.gr)
}
#should also get a 'disreg in both' to see if the marks are particularly disreged at genes disreged in >1 cond
setwd('~/Documents/Projects/Thymocyte_HiC/DEseq2_in_vivo_Development_YaRNAseq/')

wtko <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_DKO_CD69nDPWT1-3_vs_CD69nDPDKO1-3.csv',header=T)
wtko <- wtko[!is.na(wtko$padj),]

wtwt <- read.csv('04_CD69nDPWT1-3_vs_CD69pCD4SPWT1-2.csv',header=T)
wtwt <- wtwt[!is.na(wtwt$padj),] 

setwd('~/Documents/Projects/Thymocyte_HiC/DEseq2_in_vivo_Development_YaRNAseq/')

wtctcf <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_CTCFKO_CD69nDPWT4-5_vs_CD69nDPCTCFKO1-2.csv',header=T)
wtctcf <- wtctcf[!is.na(wtctcf$padj),]

wtrad21 <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_Rad21KO_DPWT1-2_vs_DPRad21KO1-2.csv',header=T)
wtrad21 <- wtrad21[!is.na(wtrad21$padj),] 

conds <- list(dev=wtwt,dko=wtko,ctcf_ko=wtctcf,rad21_ko=wtrad21)
conds.gr.list <- lapply(conds,difexp_to_granges)
#####
#####resizing from center is bad as is from center of gene?
resize_list <- function(gr){
  x <- resize(gr,11000,fix='center')
  #x <- GenomicRanges::shift(x, -2500)
  return(x)
}

around_genes <- lapply(conds.gr.list,resize_list)
#around_dev_genes <- resize(dev_changes.gr,50000,fix='center')
#around_ko_genes <- resize(ko_changes.gr,50000,fix='center')

ctcf_file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/CTCFCD69negDPWTR1.bw'
rad21_file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/Rad21CD69negDPWTR2.bw'
nipbl_file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/NipblDPWTR1_treat_pileup.bw'
#rad21_CTCFko <- ''

wt.h3k27ac_file <- "~/Documents/Projects/Thymocyte_HiC/chip_seqs/H3K27acCD69negDPWTR1R2_treat.bw"
ko.h3k27ac_file <- "~/Documents/Projects/Thymocyte_HiC/chip_seqs/H3K27acCD69negDPDKOR1R2_treat.bw"

wt.h3k27me3_file <- "~/Documents/Projects/Thymocyte_HiC/chip_seqs/H3K27me3CD69negDPWTR1R2_treat.bw"
ko.h3k27me3_file <- "~/Documents/Projects/Thymocyte_HiC/chip_seqs/H3K27me3CD69negDPDKOR1R2_treat.bw"
wt.h3k4me3_file <- "~/Documents/Projects/Thymocyte_HiC/chip_seqs/H3K4me3CD69negDPWTR1R2_treat.bw"
ko.h3k4me3_file <- "~/Documents/Projects/Thymocyte_HiC/chip_seqs/H3K4me3CD69negDPDKOR1R2.bw"


ctcfko.rad21_file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/Rad21CD69negDPCTCFKOR1.bw'
wt4.rad21_file <- '~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/6_Rad21CD69pos4SPWTR2_mapped_sorted_RemoveDuplicates_MACS_bedGraph/Rad21CD69pos4SPWTR1.bw'
wt4.rad21_file <- '~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/6_Rad21CD69pos4SPWTR2_mapped_sorted_RemoveDuplicates_MACS_bedGraph/Rad21CD69pos4SPWTR1.bw'


marks <- list(CTCF = ctcf_file,RAD21 = rad21_file, ctcfko.RAD21 = ctcfko.rad21_file ,NIPBL = nipbl_file,
              wt.h3k27ac = wt.h3k27ac_file ,ko.h3k27ac = ko.h3k27ac_file ,wt.H3K27me3 = wt.h3k27me3_file,
              ko.H3K27me3 = ko.h3k27me3_file, wt.h3k4me3 = wt.h3k4me3_file, ko.h3k4me3 = ko.h3k4me3_file)

#for each mark, want scoreMatrix list for each condition

# have heatmaps/metaplots for all active genes, demonstrate that is at most genes
all.genes.rad21 <- ScoreMatrixBin(rad21_file,resize(conds.gr.list$dev,2000,fix='center'),bin.num = 101,type="bw",strand.aware = T)
all.genes.ctcf <- ScoreMatrixBin(ctcf_file,resize(conds.gr.list$dev,2000,fix='center'),bin.num = 101,type="bw",strand.aware = T)



#all.genes.ctcf@.Data <- all.genes.ctcf@.Data[order(scaleScoreMatrix(all.genes.ctcf)@.Data[,16] , decreasing = T),]
#all.genes.rad21@.Data <- all.genes.rad21@.Data[order(scaleScoreMatrix(all.genes.rad21)@.Data[,26] , decreasing = T),]
#all.genes.ctcf@.Data <- all.genes.ctcf@.Data[order(scaleScoreMatrix(all.genes.ctcf)@.Data[,26] , decreasing = T),]
#all.genes.rad21@.Data <- all.genes.ctcf@.Data[order(rowSums(rad21.center), decreasing = T),]

heatMatrix(all.genes.ctcf, xcoords = seq(-1,1,length.out = 101), col=brewer.pal(9,"Blues"),
                winsorize = c(1,99), xlab = 'CTCF ChIP-seq signal',order = T) 
heatMatrix(all.genes.rad21, xcoords = seq(-1,1,length.out = 101), col=brewer.pal(9,"Reds"),
           winsorize = c(1,99), xlab = 'RAD21 ChIP-seq signal',order=T)

plotMeta(all.genes.ctcf, xcoords = seq(-1,1,length.out = 101),dispersion = 'se',line.col = rainbow(20)[13],dispersion.col = rainbow(20,alpha=0.5)[13],
           winsorize = c(1,97.5),xlab='Distance from Promoter (kb) (strand orientated)',ylab='DP WT CTCF ChIP-seq Signal')
abline(v=c(0),lty=2)  
title('CTCF signal at promoters',cex.main=1)
plotMeta(all.genes.rad21, xcoords = seq(-1,1,length.out = 101), dispersion='se', line.col = rainbow(20)[1],dispersion.col = rainbow(20,alpha=0.5)[1],
           winsorize = c(1,97.5),xlab='Distance from Promoter (kb) (strand orientated)',ylab='DP WT RAD21 ChIP-seq Signal')
abline(v=c(0),lty=2)
title('RAD21 signal at promoters',cex.main=1)
#quick takeaway -  way more promoters (regardless of if deregulated etc) have rad21 binding than ctcf binding

dev.off()
# barplots of how many genes overlap rad21/ctcf sites 
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

dev.barplot <- rbind(getbars(conds.gr.list$dev[conds.gr.list$dev$padj < 0.05 & conds.gr.list$dev$log2FoldChange > 0], 'Up'),
      getbars(conds.gr.list$dev[conds.gr.list$dev$padj > 0.05], 'Stable'),
      getbars(conds.gr.list$dev[conds.gr.list$dev$padj < 0.05 & conds.gr.list$dev$log2FoldChange < 0], 'Down'))

dev.barplot.gg <- ggplot(dev.barplot, aes(x=d, y=y, fill=x)) + 
  geom_bar(position='fill', stat="identity")  + coord_flip() + xlab('Gene Status') + ylab('Proportion') +
  ggtitle('Distribution of RAD21 and CTCF at promoters of Developmentally Regulated Genes') + scale_fill_jco() +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 13),plot.title = element_text(size = 12),legend.text = element_text(size = 12))

ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/dev.barplot.pdf',
        dev.barplot.gg,width=7,height=2,units = 'in',device = cairo_pdf)

dko.barplot <- rbind(getbars(conds.gr.list$dko[conds.gr.list$dko$padj < 0.05 & conds.gr.list$dko$log2FoldChange > 0], 'Up'),
                     getbars(conds.gr.list$dko[conds.gr.list$dko$padj > 0.05], 'Stable'),
                     getbars(conds.gr.list$dko[conds.gr.list$dko$padj < 0.05 & conds.gr.list$dko$log2FoldChange < 0], 'Down'))

dko.barplot.gg <- ggplot(dko.barplot, aes(x=d, y=y, fill=x)) + 
  geom_bar(position='fill', stat="identity") + coord_flip() + xlab('Gene Status') + ylab('Proportion') +
  ggtitle('Distribution of RAD21 and CTCF at promoters of Genes Deregulated upon ΔCTCF/ΔRAD21') + scale_fill_jco() +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 13),plot.title = element_text(size = 11),legend.text = element_text(size = 12))


ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/dko.barplot.pdf',
        dko.barplot.gg,width=7,height=2,units = 'in',device = cairo_pdf)

ctcfko.barplot <- rbind(getbars(conds.gr.list$ctcf_ko[conds.gr.list$ctcf_ko$padj < 0.05 & conds.gr.list$ctcf_ko$log2FoldChange > 0], 'Up'),
                     getbars(conds.gr.list$ctcf_ko[conds.gr.list$ctcf_ko$padj > 0.05], 'Stable'),
                     getbars(conds.gr.list$ctcf_ko[conds.gr.list$ctcf_ko$padj < 0.05 & conds.gr.list$ctcf_ko$log2FoldChange < 0], 'Down'))

ctcfko.barplot.gg <-  ggplot(ctcfko.barplot, aes(x=d, y=y, fill=x)) + 
  geom_bar(position='fill', stat="identity") + coord_flip() + xlab('Gene Status') + ylab('Proportion') +
  ggtitle('Distribution of RAD21 and CTCF at promoters of Genes Deregulated upon ΔCTCF') + scale_fill_jco() +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 13),plot.title = element_text(size = 12),legend.text = element_text(size = 12))


ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/ctcfko.barplot.pdf',
        ctcfko.barplot.gg,width=7,height=2,units = 'in',device = cairo_pdf)


rad21ko.barplot <- rbind(getbars(conds.gr.list$rad21_ko[conds.gr.list$rad21_ko$padj < 0.05 & conds.gr.list$rad21_ko$log2FoldChange > 0], 'Up'),
                     getbars(conds.gr.list$rad21_ko[conds.gr.list$rad21_ko$padj > 0.05], 'Stable'),
                     getbars(conds.gr.list$rad21_ko[conds.gr.list$rad21_ko$padj < 0.05 & conds.gr.list$rad21_ko$log2FoldChange < 0], 'Down'))

rad21ko.barplot.gg <- ggplot(rad21ko.barplot, aes(x=d, y=y, fill=x)) + 
  geom_bar(position='fill', stat="identity")  + coord_flip()  + xlab(' Gene Status') + ylab('Proportion') +
  ggtitle('Distribution of RAD21 and CTCF at promoters of Genes Deregulated upon ΔRAD21') + scale_fill_jco() +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 13),plot.title = element_text(size = 12),legend.text = element_text(size = 12))


ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/rad21ko.barplot.pdf',
        rad21ko.barplot.gg,width=7,height=2,units = 'in',device = cairo_pdf)



plot(permTest(A=g[g$padj < 0.05 & g$log2FoldChange < 0],B=rad21_bed,randomize.function = resampleRegions,
              evaluate.function = numOverlaps,universe=g))
plot(permTest(A=g[g$padj < 0.05 & g$log2FoldChange > 0],B=rad21_bed,randomize.function = resampleRegions,
              evaluate.function = numOverlaps,universe=g))

plot(permTest(A=g[g$padj < 0.05 & g$log2FoldChange < 0],B=ctcf_bed,randomize.function = resampleRegions,
              evaluate.function = numOverlaps,universe=g))
plot(permTest(A=g[g$padj < 0.05 & g$log2FoldChange > 0],B=ctcf_bed,randomize.function = resampleRegions,
              evaluate.function = numOverlaps,universe=g))

#takeaway - dev reg genes do have different amounts of rad21 binding
# but is mostly number that have rad21 *only* binding (not ctcf as well) - mostly not significant for
#ctcfko dereg genes


#
CairoPDF('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/ctcf_rad21.at.genes.mps.pdf')


fold_change <- 0
get_scorematrix <- function(positions){
  up.sm <- ScoreMatrixBin(mark,positions[positions$padj < 0.05 & positions$log2FoldChange > fold_change],bin.num = 51,type="bw",strand.aware = T)
  down.sm <- ScoreMatrixBin(mark,positions[positions$padj < 0.05 & positions$log2FoldChange < -fold_change],bin.num = 51,type="bw",strand.aware = T)
  noreg.sm <- ScoreMatrixBin(mark,positions[positions$padj > 0.05 ],bin.num = 51,type="bw",strand.aware = T)
  
  sml <- ScoreMatrixList(c( up.sm, down.sm,noreg.sm))
  return(sml)
}

sml.list <- vector("list",length=length(marks))
names(sml.list) <- names(marks)

cond_name <- list('dev'= 'CD69pSP WT','dko'= 'ΔCTCF/ΔRAD21','ctcf_ko'='ΔCTCF','rad21_ko'='ΔRAD21')

for(m in names(marks)){
  mark = marks[[m]]
  mark.sml <- lapply(around_genes, get_scorematrix)
  sml.list[[m]] <- mark.sml #I guess
}

for(m in names(marks)){
  for(cond in names(conds)){
    plotMeta(sml.list[[m]][[cond]],winsorize = c(1, 99),smoothfun=function(x) stats::lowess(x, f = 0.01) ,dispersion = 'se',
             profile.names = c('Upregulated','Downregulated','Stable'),xlab='Distance from Promoter (kb) (strand orientated)',
             ylab=paste0('CD69negDPWT ChIP-seq Signal (',m,')'),xcoords = seq(-5,5,length.out = 51))
             # cex.lab=0.5,cex.axis=0.5) 
    #axis(1,at=c(26),labels=c("Promoter Start"))
    abline(v=c(0),lty=2)
    title(paste0('Differentially Expressed genes in ',cond_name[[cond]])) #,cex.main = 0.5)
  }
} #line.col maybe change brewer, adding titles

dev.off()
#
heatMatrix(sml.list$ko.H3K27me3$dko[[2]],winsorize = c(1,99),clustfun = function(x){kmeans(x,centers=3)$cluster})
heatMatrix(sml.list$wt.H3K27me3$dko[[2]],winsorize = c(1,99.5),order = T,xcoords = c(-15,15))
heatMatrix(sml.list$ko.H3K27me3$dko[[2]],winsorize = c(1,98.5),order = T)

heatMatrix(sml.list$wt.H3K27me3$dko[[1]],winsorize = c(1,99.5),order = T,xcoords = c(-15,15),
           xlab = 'Distance from Promoter (kb)', main = 'H3K27me3 at Upregulated Genes')
heatMatrix(sml.list$wt.H3K27me3$dko[[2]],winsorize = c(1,99.5),order = T,xcoords = c(-15,15),
           xlab = 'Distance from Promoter (kb)', main = 'H3K27me3 at Downregulated Genes')
heatMatrix(sml.list$wt.H3K27me3$dko[[3]],winsorize = c(1,99.8),order = T,xcoords = c(-15,15),
           xlab = 'Distance from Promoter (kb)', main = 'H3K27me3 at Stable Genes')


#plotMeta(sml.list$ctcf$dev,winsorize = c(1, 99),smoothfun=function(x) stats::lowess(x, f = 0.01) ,dispersion = 'se',
#         profile.names = c('Upreg','Downreg','None'),xlab='Promoter Position (strand orientated)',ylab='Average Coverage') #dispersion.col = 'grey',
#axis(1,at=c(26),labels=c("Promoter Start"))
#abline(v=c(26),lty=2)

# sm@.Data[,26] to get mark values at each promoter
# so can then visulaise the distribution. Boxplotsdefault but ecdf would look good also I think

scores_at_peaks <- function(sml_list,mark,cond){  
  up <- as.data.frame(sml_list[[mark]][[cond]][[1]]@.Data[,16])
  colnames(up) <- 'score'
  up$status <- 'Upregulated'
  
  down <- as.data.frame(sml_list[[mark]][[cond]][[2]]@.Data[,16])
  colnames(down) <- 'score'
  down$status <- 'Downregulated'
  
  none <- as.data.frame(sml_list[[mark]][[cond]][[3]]@.Data[,16])
  colnames(none) <- 'score'
  none$status <- 'Stable'
  
  peak_scores <- rbind(up,down,none)
  #peak_scores$status <- factor(peak_scores$status, levels = c('Upregulated','Downregulated','Not Significant'),ordered = TRUE)
  return(peak_scores)
}
#cond <-- eg dev, dko
#want to keep upregulated - red. downregulated - green. non-changing - blue 

my_comparisons <- list( c("Downregulated", "Stable"), c("Downregulated", "Upregulated"), c("Stable", "Upregulated") )

ctcf_ylim = 10
ctcfpos = c(8.4,9.2,10)
rad21_ylim = 3.5
rad21pos = c(2.8,3.1,3.45)
#
peak_scores <- scores_at_peaks(sml.list,'CTCF','dev')
#peak_scores$score[peak_scores$score > quantile(peak_scores$score,0.98)] = quantile(peak_scores$score,0.98)
ctcf_dev_peakboxplot <- ggplot(peak_scores,aes(x=status,y=score)) + geom_boxplot(outlier.shape = NA) + theme_classic() + 
  xlab('CD69nDPWT vs CD69pSPWT') + ylab('CD69negDPWT ChIP-seq Signal (CTCF)') + 
  stat_compare_means(comparisons = my_comparisons,label.y = ctcfpos,tip.length = 0.004) + coord_cartesian(ylim=c(0,ctcf_ylim))  

peak_scores <- scores_at_peaks(sml.list,'RAD21','dev')
#peak_scores$score[peak_scores$score > quantile(peak_scores$score,0.98)] = quantile(peak_scores$score,0.98)
rad21_dev_peakboxplot <- ggplot(peak_scores,aes(x=status,y=score)) + geom_boxplot(outlier.shape = NA) + theme_classic() + 
  xlab('CD69nDPWT vs CD69pSPWT') + ylab('CD69negDPWT ChIP-seq Signal (RAD21)') + 
  stat_compare_means(comparisons = my_comparisons,label.y = rad21pos,tip.length = 0.01) + coord_cartesian(ylim=c(0,rad21_ylim))

#

peak_scores <- scores_at_peaks(sml.list,'CTCF','dko')
#peak_scores$score[peak_scores$score > quantile(peak_scores$score,0.98)] = quantile(peak_scores$score,0.98)
ctcf_dko_peakboxplot <- ggplot(peak_scores,aes(x=status,y=score)) + geom_boxplot(outlier.shape = NA) + theme_classic() + 
  xlab('CD69nDPWT vs CD69nDP ΔCTCF/ΔRAD21') + ylab('CD69negDPWT ChIP-seq Signal (CTCF)') + 
  stat_compare_means(comparisons = my_comparisons,label.y = ctcfpos,tip.length = 0.004) + coord_cartesian(ylim=c(0,ctcf_ylim))  

#

peak_scores <- scores_at_peaks(sml.list,'RAD21','dko')
#peak_scores$score[peak_scores$score > quantile(peak_scores$score,0.98)] = quantile(peak_scores$score,0.98)
rad21_dko_peakboxplot <- ggplot(peak_scores,aes(x=status,y=score)) + geom_boxplot(outlier.shape = NA) + theme_classic() +
  xlab('CD69nDPWT vs CD69nDP ΔCTCF/ΔRAD21') + ylab('CD69negDPWT ChIP-seq Signal (RAD21)') + 
  stat_compare_means(comparisons = my_comparisons,label.y = rad21pos,tip.length = 0.01) + coord_cartesian(ylim=c(0,rad21_ylim))

#
peak_scores <- scores_at_peaks(sml.list,'CTCF','ctcf_ko')
#peak_scores$score[peak_scores$score > quantile(peak_scores$score,0.98)] = quantile(peak_scores$score,0.98)
ctcf_ctcfko_peakboxplot <- ggplot(peak_scores,aes(x=status,y=score)) + geom_boxplot(outlier.shape = NA) + theme_classic() + 
  xlab('CD69nDPWT vs CD69nDP ΔCTCF') + ylab('CD69negDPWT ChIP-seq Signal (CTCF)') + 
  stat_compare_means(comparisons = my_comparisons,label.y = ctcfpos,tip.length = 0.004) + coord_cartesian(ylim=c(0,ctcf_ylim)) 

peak_scores <- scores_at_peaks(sml.list,'RAD21','ctcf_ko')
#peak_scores$score[peak_scores$score > quantile(peak_scores$score,0.98)] = quantile(peak_scores$score,0.98)
rad21_ctcfko_peakboxplot <- ggplot(peak_scores,aes(x=status,y=score)) + geom_boxplot(outlier.shape=NA) + theme_classic() +
  xlab('CD69nDPWT vs CD69nDP ΔCTCF') + ylab('CD69negDPWT ChIP-seq Signal (RAD21)') +
  stat_compare_means(comparisons = my_comparisons,label.y = rad21pos,tip.length = 0.01) + coord_cartesian(ylim=c(0,rad21_ylim))

#
peak_scores <- scores_at_peaks(sml.list,'CTCF','rad21_ko')
#peak_scores$score[peak_scores$score > quantile(peak_scores$score,0.98)] = quantile(peak_scores$score,0.98)
ctcf_rad21ko_peakboxplot <- ggplot(peak_scores,aes(x=status,y=score)) + geom_boxplot(outlier.shape = NA) + theme_classic() + 
  xlab('CD69nDPWT vs CD69nDP ΔRAD21') + ylab('CD69negDPWT ChIP-seq Signal (CTCF)') + 
  stat_compare_means(comparisons = my_comparisons,label.y = ctcfpos,tip.length = 0.004) + coord_cartesian(ylim=c(0,ctcf_ylim)) 

peak_scores <- scores_at_peaks(sml.list,'RAD21','rad21_ko')
#peak_scores$score[peak_scores$score > quantile(peak_scores$score,0.98)] = quantile(peak_scores$score,0.98)
rad21_rad21ko_peakboxplot <- ggplot(peak_scores,aes(x=status,y=score)) + geom_boxplot(outlier.shape = NA) + theme_classic() + 
  xlab('CD69nDPWT vs CD69nDP ΔRAD21') + ylab('CD69negDPWT ChIP-seq Signal (RAD21)') + 
  stat_compare_means(comparisons = my_comparisons,label.y = rad21pos,tip.length = 0.01) + coord_cartesian(ylim=c(0,rad21_ylim)) 

peak_scores <- scores_at_peaks(sml.list,'wt.H3K27me3','dko')
#peak_scores$score[peak_scores$score > quantile(peak_scores$score,0.98)] = quantile(peak_scores$score,0.98)
ggplot(peak_scores,aes(x=status,y=score)) + geom_boxplot(outlier.shape = NA) + theme_classic() + 
  xlab('CD69nDPWT vs CD69nDP ΔCTCF/ΔRAD21') + ylab('CD69negDPWT ChIP-seq Signal (H3K27me3)') +
  stat_compare_means(comparisons = my_comparisons,label.y = c(0.35,0.4,0.45),tip.length = 0.01) + coord_cartesian(ylim=c(0,0.45)) 

peak_scores <- scores_at_peaks(sml.list,'ko.H3K27me3','dko')
#peak_scores$score[peak_scores$score > quantile(peak_scores$score,0.98)] = quantile(peak_scores$score,0.98)
ggplot(peak_scores,aes(x=status,y=score)) + geom_boxplot(outlier.shape = NA) + theme_classic() +
  xlab('CD69nDPWT vs CD69nDP ΔCTCF/ΔRAD21') + ylab('CD69negDP ΔCTCF/ΔRAD21 ChIP-seq Signal (H3K27me3)') +
  stat_compare_means(comparisons = my_comparisons,label.y = c(0.35,0.4,0.45),tip.length = 0.01) + coord_cartesian(ylim=c(0,0.45)) 



f4.mark.peak.boxplots <- plot_grid( rad21_ctcfko_peakboxplot,ctcf_ctcfko_peakboxplot,rad21_rad21ko_peakboxplot,
                                    ctcf_rad21ko_peakboxplot,ncol = 2 )
f4.mark.peak.boxplots2 <- plot_grid( rad21_dev_peakboxplot,ctcf_dev_peakboxplot,rad21_dko_peakboxplot,
                                    ctcf_dko_peakboxplot,ncol = 2 )


ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/ctcf_rad21.at.genes.boxplots.pdf',
        f4.mark.peak.boxplots,width=7,height=10,units = 'in',device = cairo_pdf)
ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/ctcf_rad21.at.genes.boxplots2.pdf',
        f4.mark.peak.boxplots2,width=7,height=10,units = 'in',device = cairo_pdf)

#### quickly get *different in both* (as in up/up or down/down reg in DKO and dev), is more intense?
####z
####

changing_mult_condition <- function(x,y,mark){
  Agenes <- around_genes[[x]]
  Bgenes <- around_genes[[y]]
  #A_up.id <- Agenes[Agenes$padj < 0.05 & Agenes$log2FoldChange >0]$ens_id
  B_up.id <- Bgenes[Bgenes$padj < 0.05 & Bgenes$log2FoldChange >0]$ens_id
  up.both <- Agenes[Agenes$padj < 0.05 & Agenes$log2FoldChange > 0 & Agenes$ens_id %in% B_up.id]
  up.sm <- ScoreMatrixBin(mark,up.both,bin.num = 31,type="bw",strand.aware = T)
  
  B_down.id <- Bgenes[Bgenes$padj < 0.05 & Bgenes$log2FoldChange <0]$ens_id
  down.both <- Agenes[Agenes$padj < 0.05 & Agenes$log2FoldChange < 0 & Agenes$ens_id %in% B_down.id]
  down.sm <- ScoreMatrixBin(mark,down.both,bin.num = 31,type="bw",strand.aware = T)
  
  B_neither.id <- Bgenes[Bgenes$padj > 0.05]$ens_id
  neither.both <- Agenes[Agenes$padj > 0.05 & Agenes$ens_id %in% B_neither.id]
  neither.sm <- ScoreMatrixBin(mark,neither.both,bin.num = 31,type="bw",strand.aware = T)
  
  return(list(c(up.sm,down.sm,neither.sm)))
}

dev_dko_changing.CTCF <- changing_mult_condition('dev','dko',marks$CTCF)
dev_dko_changing.RAD21 <- changing_mult_condition('dev','dko',marks$RAD21)

plotMeta(dev_dko_changing.CTCF[[1]],winsorize = c(1, 99),smoothfun=function(x) stats::lowess(x, f = 0.01) ,dispersion = 'se',
         profile.names = c('Upregulated','Downregulated','Stable'),xlab='Distance from Promoter (kb) (strand orientated)',
         ylab='CD69negDPWT ChIP-seq Signal CTCF', xcoords = seq(-5,5,length.out = 31)) #dispersion.col = 'grey',
abline(v=c(0),lty=2)


plotMeta(dev_dko_changing.RAD21[[1]],winsorize = c(1, 99),smoothfun=function(x) stats::lowess(x, f = 0.01) ,dispersion = 'se',
         profile.names = c('Upregulated','Downregulated','Stable'),xlab='Distance from Promoter (kb) (strand orientated)',
         ylab='CD69negDPWT ChIP-seq Signal RAD21',xcoords = seq(-5,5,length.out = 31)) #dispersion.col = 'grey',
abline(v=c(0),lty=2)



#########3
#########
# heatmaps sorted either by score around center or by chip of interest
# e.g. overlapping rad21 chips

rad21.tmm.wt.file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/bams/normFactor_normedBams/5_Rad21CD69negDPWTR1_mapped_sorted_RemoveDuplicates_TMM.bigwig'
rad21.tmm.ctcfko.file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/bams/normFactor_normedBams/Rad21CD69negDPCTCFKOR1_mapped_sorted_RemoveDuplicates_TMM.bigwig'
rad21.tmm.wt4.file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/bams/normFactor_normedBams/6_Rad21CD69pos4SPWTR2_mapped_sorted_RemoveDuplicates_TMM.bigwig'


ctcf_peaks <- makeGRangesFromDataFrame(read.table('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/1_CTCFCD69negDPWTR1/1_CTCFCD69negDPWTR1_peaks.subpeaks.bed'),seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)
rad21_peaks.wt <- makeGRangesFromDataFrame(read.table('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/Rad21_ChIPseq_in_CD69negDP/WT_CD69negDP/Rad21CD69negDPWTR2_peaks.narrowPeak.bedGraph'),seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)
rad21_peaks.ko <- makeGRangesFromDataFrame(read.table('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/Rad21_ChIPseq_in_CD69negDP/CTCFKO_CD69negDP/Rad21CD69negDPCTCFKOR1_peaks.narrowPeak.bedGraph'),seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)
#need to remove chr x, y m
#maybe filter by strength a bit more
peaks= GRangesList(ctcf = ctcf_peaks , wt.rad21 = rad21_peaks.wt, ko.rad21 = rad21_peaks.ko)

tf.comb = findFeatureComb(peaks, width = 1000)
tf.comb = tf.comb[order(tf.comb$class)]

comb.list <- vector(mode = "list", length = 7) 
for(i in 1:7){
  comb.list[[i]] <- as.integer(rownames(as.data.frame(tf.comb)))[tf.comb$class == i] 
}


bw.files = List(ctcf = ctcf_file,rad21.wt = rad21_file ,rad21.ctcfko = ctcfko.rad21_file)
mypal= ggsci::pal_gsea("default", alpha = 0.7)(9)

make.sml <- function(bw,windows= resize(tf.comb,1250,fix='center'),bin.num=81,filetype='bw'){
  return(ScoreMatrixBin(bw,windows,bin.num = bin.num,type = filetype ))
}

sml = ScoreMatrixList(lapply(bw.files, make.sml)) #ScoreMatrixList(bw.files, tf.comb, bin.num = 20, type = "bw")
#names(sml) = sampleInfo$sampleName[match(names(sml), sampleInfo$fileName)]
#sml.scaled = scaleScoreMatrixList(sml)
#multiHeatMatrix(scaleScoreMatrixList(sml), xcoords = c(-750, 750), col=mypal,matrix.main = c('CTCF','wt.RAD21','ctcfko.RAD21'),
#                winsorize = c(2,97.5),group = comb.list,order = order(rowSums(sml$rad21.wt)))
#multiHeatMatrix(sml, xcoords = c(-750, 750), col=mypal,matrix.main = c('CTCF','wt.RAD21','ctcfko.RAD21'),
#                winsorize = c(2,97.5),group = comb.list,order = order(rowSums(sml$rad21.wt)))
multiHeatMatrix(sml, xcoords = c(-750, 750), col=brewer.pal(9,"Blues"),matrix.main = c('CTCF','WT RAD21','ΔCTCF RAD21'),
                winsorize = c(2,97.5),group = comb.list,order = order(rowSums(sml$rad21.wt)))

#maybe not scaled?
#maybe should have been done with more windows
#I want to get a common scale for WT vs KO, makes the consistent loss a lot clearer
#maybe remove group 5 (no wt rad21 binding, ctcf and ko.rad21 binding) - looks fake? But is sus to remove

#multiHeatMatrix(ScoreMatrixList(list(sml[[1]],sml[[2]],sml[[3]])), xcoords = c(-750, 750), col=brewer.pal(9,"Blues"),matrix.main = c('woof','wt.RAD21','ctcfko.RAD21'),
#                winsorize = c(1,97.5),group = comb.list,order = order(rowSums(sml$rad21.wt)),common.scale = T)

#multiHeatMatrix(ScoreMatrixList(list(sml[[2]],sml[[2]],sml[[3]])), xcoords = c(-750, 750), col=brewer.pal(9,"Blues"),matrix.main = c('woof','wt.RAD21','ctcfko.RAD21'),
#                winsorize = c(1,97.5),group = comb.list,order = order(rowSums(sml$rad21.wt)),common.scale = T)


#cl1 <- function(x) kmeans(x, centers=10)$cluster
#multiHeatMatrix(sml,order=T, xcoords = c(-500, 500), col=brewer.pal(9,"Blues"),
#                matrix.main = c('CTCF','wt.RAD21','ctcfko.RAD21'),winsorize = c(1,95),clustfun = cl1)
#could try other cluster methods.
#quite like this - can see a group of high rad21/ctcf that has minorly weakened ctcfko rad21
#
#or pulling out clusters- meta profile them
#plotMeta(ScoreMatrixList(c(sml$rad21.wt,sml$rad21.ctcfko)),profile.names = c('WT','KO'))
#all in all, only a mild reduction in CTCF at all binding sites
# we can see that sites with very weak CTCF tend to lose RAD21 binding
# can kiiiind of see a group of regions that have no wt rad21 and ctcf, and gain a little ctcf after KO
# and maybe a block of regions with no CTCF, that have a little WT rad21 and gain a little more after KO

#to plot alongside the above heatmap
make_sml.meta <- function(bws,windows= resize(tf.comb,1250,fix='center'),bin.num=41,filetype='bw'){
  x <- vector(mode = "list", length = 7)
  for(i in 1:7){
    win = windows[windows$class == i]
    getmetas <- function(bw){
      return(ScoreMatrixBin(bw,win,bin.num = bin.num))
    }
    x[[i]] <-lapply(bws,getmetas)
  }
  return(x)
  #return(ScoreMatrixBin(bw,windows,bin.num = bin.num,type = filetype ))
  }

x <- make_sml.meta(list(rad21_file,ctcfko.rad21_file))
for(i in 1:7){
  plotMeta(ScoreMatrixList(x[[i]]),winsorize = c(0,97.5),ylim = c(0,3),
           profile.names = c('WT','ΔCTCF'),xlab='Distance from CHiP peak (bp))',
           ylab='RAD21 ChIP-seq Signal', xcoords = seq(-750,750,length.out = 41),dispersion='se',dispersion.col = c('red','blue')) #dispersion.col = 'grey',
  abline(v=c(0),lty=2)
  
}

#ok so super similar, but just for all wt RAD21 binding sites, ordered by wt binding strength I think
#just as a 'show being reduced but not super massively'  

wt.peaks = resize(rad21_peaks.wt[width(rad21_peaks.wt) < 600],2000,fix = 'center')
wt.peaks <- wt.peaks[order(wt.peaks$V5,decreasing = T)]
#wt.peaks <- wt.peaks[order(width(wt.peaks),decreasing = T)]

rad21.binding = ScoreMatrixList(list(ScoreMatrixBin(rad21_file,wt.peaks,bin.num = 121 ),
     ScoreMatrixBin(ctcfko.rad21_file,wt.peaks,bin.num = 121 ),ScoreMatrixBin(wt4.rad21_file,wt.peaks,bin.num = 121 )))

rad21.binding <- intersectScoreMatrixList(rad21.binding)

#rad21.binding.centerA = rad21.binding[[1]]@.Data[,c(19,20,21,22,23)]
#rad21.binding.centerAB = rad21.binding[[2]]@.Data[,c(19,20,21,22,23)]

#rad21.binding[[1]]@.Data <- rad21.binding[[1]]@.Data[order(rowSums(rad21.binding.centerA) , decreasing = T),]
#rad21.binding[[2]]@.Data <- rad21.binding[[2]]@.Data[order(rowSums(rad21.binding.centerA), decreasing = T),]
#rad21.binding[[3]]@.Data <- rad21.binding[[3]]@.Data[order(rowSums(rad21.binding.centerA), decreasing = T),]


multiHeatMatrix(ScoreMatrixList(c(rad21.binding[1],rad21.binding[2])), xcoords = c(-1000, 1000), col=brewer.pal(9,"Blues"),matrix.main = c('WT RAD21','CTCFko RAD21'),
                winsorize = c(1,97.5),order = T,common.scale = T)

multiHeatMatrix(rad21.binding, xcoords = c(-1000, 1000), col=brewer.pal(9,"Blues"),matrix.main = c('DP','ΔCTCF','SP'),
                winsorize = c(1,97.5),order = F,common.scale = F)

plotMeta(ScoreMatrixList(c(rad21.binding[1],rad21.binding[2])),winsorize = c(1, 99),smoothfun=function(x) stats::lowess(x, f = 0.01) ,dispersion = 'se',
         profile.names = c('WT','CTCFKO'),xlab='Distance from Promoter (bp) (strand orientated)',
         ylab='RAD21 ChIP signal',xcoords = seq(-1000,1000,by = 2000/120),
         line.col = rainbow(20)[c(1,13)],dispersion.col = rainbow(20,alpha=0.25)[c(1,13)],overlay=F) #dispersion.col = 'grey',
abline(v=c(0),lty=2)


normWT1 = ScoreMatrixBin(rad21.tmm.wt.file,wt.peaks,bin.num = 121 )
normCTCFko = ScoreMatrixBin(rad21.tmm.ctcfko.file,wt.peaks,bin.num = 121 )
normWT4 = ScoreMatrixBin(rad21.tmm.wt4.file,wt.peaks,bin.num = 121 )

plotMeta(ScoreMatrixList(c(normWT1,normCTCFko,normWT4)),winsorize = c(1, 99),smoothfun=function(x) stats::lowess(x, f = 0.01) ,dispersion = 'se',
         profile.names = c('DP','CTCFKO','SP'),xlab='Distance from Promoter (bp) (strand orientated)',
         ylab='RAD21 ChIP signal',xcoords = seq(-1000,1000,by = 2000/120), overlay =F) #,dev.off()
          #dispersion.col = 'grey',
abline(v=c(0),lty=2)

multiHeatMatrix(ScoreMatrixList(list(normWT1,normCTCFko,normWT4)), xcoords = c(-1000, 1000), col=brewer.pal(9,"Blues"),
                matrix.main = c('DP','ΔCTCF','SP'), winsorize = c(1,97.5),order = F,common.scale = T)
# same as above, but now pull in h3k27me3 or similar?
#should also look at properly getting deferentially bound- and then seeing whether at or near 
#difreg genes?

#maintained rad21 sites- wider (have signal at adjacent bins? Also wonder if the back to back sites |  |
#of adjacent tads would be looking out for each other?- bot hhave to to be lost?
#also can take loop domain with motif locations and use as input
