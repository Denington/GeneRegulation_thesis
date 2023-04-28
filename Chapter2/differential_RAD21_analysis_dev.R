#### relationship of differential rad21 binding (vs CD4+SP and CTCFko)
## on gene expression, change in domain strength/insulation
## and is distance from gained peaks to lost peaks sig closer than expected? 
#gained rad21 sites form new loops?
# # (note that should be able to use calcnormfactors for nicer scaling between rad21 samples)

#bamCoverage -b 6_Rad21CD69pos4SPWTR2_mapped_sorted_RemoveDuplicates.bam -o Rad21CD69pos4SPWTR2.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2472342367 --ignoreForNormalization chrX -bl ~/Documents/Projects/Thymocyte_HiC/Thesis/external_data/mm9_blacklist_andall_chrY.bed

#check whether differential rad21 sites are more likely found in cds and a compartments (with a universe of all sites)

source('~/Documents/Projects/Thymocyte_HiC/Thesis/scripts/general_functions.R')
library(regioneR)
library(genomation)
library(ggplot2)
library(rtracklayer)
library(regioneR)

rad21.sites <- read.table('~/Documents/Projects/Thymocyte_HiC/Differential_chip/WTvsWT4.RAD21.bed',sep = '\t',header = T)


rad21.sites <- makeGRangesFromDataFrame(rad21.sites,keep.extra.columns = T)
out.ranges <- rad21.sites

genes= resize(conds.gr.list$dev,500,fix='center')
plot(permTest(A= genes[genes$padj < 0.05 & genes$log2FoldChange > 0],  B=out.ranges[out.ranges$FDR < 0.05 & out.ranges$direction == 'up'],
              evaluate.function = numOverlaps,randomize.function= resampleRegions, universe = genes, ntimes = 100))
#upregulated genes very minorly gain RAD21 after CTCFKO, and for dev (but not many genes at all!)
#interestingly not something that gets stronger with logfold cutoff for change gene expression
#I had thought that there'd be quite a few upreg genes gaining rad21. Fits with the metaplots though 
#so this holds for all sets of genes interestingly (for ctcfko rad21) (less so for rad21ko), but sometimes in the opposite direction

plot(permTest(A= genes[genes$padj < 0.05 & genes$log2FoldChange < 0],  B=out.ranges[out.ranges$FDR < 0.05 & out.ranges$direction == 'down'],
              evaluate.function = numOverlaps,randomize.function= resampleRegions, universe = genes, ntimes = 100))
#^ a bunch of downregulated genes lose RAD21 it seems (4-fold increase). Still only accounts for 1/10 of downreg genes though
#this maintains/ gets stronger with harder cutoff also
# weirdly see the opposite for developmental genes, where upregulated genes minorly. Should check that isn't a labellig thing on my part..
# prob is. If it is, then also makes sense, otherwise -_- --> might be was looking at ctcfko genes.. which would be interesting
# for dev genes, do see downreg genes overlap downreg rad21 sites

#should note that not that many genes- especially upregulated- associated with differential RAD21 (either for CTCFko or through development)

#plot(permTest(A= genes,  B=out.ranges[out.ranges$FDR < 0.05 & out.ranges$direction == 'down'],
#              evaluate.function = numOverlaps,randomize.function= randomizeRegions, ntimes = 100, genome= 'mm9'))

#plot(permTest(A= genes,  B=out.ranges[out.ranges$FDR < 0.05 & out.ranges$direction == 'up'],
#              evaluate.function = numOverlaps,randomize.function= randomizeRegions, ntimes = 100,genome="mm9"))
#

plot(permTest(B= genes,  A=out.ranges[out.ranges$FDR < 0.05 & out.ranges$direction == 'down'],
              evaluate.function = numOverlaps,randomize.function= resampleRegions, ntimes = 100, universe= out.ranges))

plot(permTest(B= genes,  A=out.ranges[out.ranges$FDR < 0.05 & out.ranges$direction == 'up'],
              evaluate.function = numOverlaps,randomize.function= resampleRegions, ntimes = 100,universe= out.ranges))


#regions that lose RAD21 binding after CTCFKO are less likely to be at promoters
#regions that gain RAD21 binding after CTCFKO are more likely to be at promoters 
#(this is for vs ctcfko rad21. vs CD4+ rad21, tend not to be at genes? Presumably at boundaries)
#should now see if more likely to be loops, boundaries, repeats?
#as would be an obvious- RAD21 moves to *eg* genes, and these then loop to other things
# for dev, changes in rad21 binding tend to be away from genes full stop?


#what about at ctcf sites?
ctcf_bed <- import.bedGraph('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/1_CTCFCD69negDPWTR1/1_CTCFCD69negDPWTR1_peaks.subpeaks.bedgraph')

plot(permTest(B= ctcf_bed,  A=out.ranges[out.ranges$FDR < 0.05 & out.ranges$direction == 'down'],
              evaluate.function = numOverlaps,randomize.function= resampleRegions, ntimes = 100, universe= out.ranges,
              count.once=T))

plot(permTest(B= ctcf_bed,  A=out.ranges[out.ranges$FDR < 0.05 & out.ranges$direction == 'up'],
              evaluate.function = numOverlaps,randomize.function= resampleRegions, ntimes = 100, universe= out.ranges,
              count.once=T))

##########
### now look at Domain Boundaries. Probably want to give some give in position- resize corners to like 20kb or something
#should see if these regions gain loops also
#what about lost at one corner vs both?

ins_folder <- '~/Documents/Projects/Thymocyte_HiC/Thesis/Insulation/'

wt1.Ins <- paste0(ins_folder, 'CD69nDPWT.10kb.insul_score_200000.KRnorm.bw')
dko.Ins <- paste0(ins_folder, 'CD69nDPDKO.10kb.insul_score_200000.KRnorm.bw')
ctcfko.Ins <- paste0(ins_folder, 'CD69nCTCFKO.10kb.insul_score_200000.KRnorm.bw')
wt4.Ins <- paste0(ins_folder, 'CD69pSPWT.10kb.insul_score_200000.KRnorm.bw')

wt1.boundaries <- import.bw('~/Documents/Projects/Thymocyte_HiC/Thesis/Insulation/CD69nDPWT.10kb.insul_pp_200000.KRnorm.bw')
dko.boundaries <- import.bw('~/Documents/Projects/Thymocyte_HiC/Thesis/Insulation/CD69nDPDKO.10kb.insul_pp_200000.KRnorm.bw')
ctcfko.boundaries <- import.bw('~/Documents/Projects/Thymocyte_HiC/Thesis/Insulation/CD69nCTCFKO.10kb.insul_pp_200000.KRnorm.bw')
wt4.boundaries <- import.bw('~/Documents/Projects/Thymocyte_HiC/Thesis/Insulation/CD69pSPWT.10kb.insul_pp_200000.KRnorm.bw')
# filter boundaries to remove super weak boys?
wt1.boundaries <- wt1.boundaries[wt1.boundaries$score > 0.1]
up.boundaries <- wt1.boundaries[unique(queryHits(findOverlaps(wt1.boundaries,rad21.sites[rad21.sites$FDR < 0.05 & rad21.sites$direction == 'up'],maxgap = 30000)))]
down.boundaries <- wt1.boundaries[unique(queryHits(findOverlaps(wt1.boundaries,rad21.sites[rad21.sites$FDR < 0.05 & rad21.sites$direction == 'down'],maxgap = 30000)))]
stab.boundaries <- wt1.boundaries[unique(queryHits(findOverlaps(wt1.boundaries,rad21.sites[rad21.sites$FDR > 0.05],maxgap = 30000)))]
#lots of overlap hmm, might want to get unique and non-unique 
#up.boundaries <- up.boundaries[-queryHits(findOverlaps(up.boundaries,stab.boundaries))]
#down.boundaries <- up.boundaries[-queryHits(findOverlaps(down.boundaries,stab.boundaries))]
#but so many e.g. up boundaries have real close stab boundaries
#oh actually, might be because some of the boundaries are very weak- maybe subset

wtU <- ScoreMatrixBin(wt1.Ins,resize(up.boundaries,250000,fix='center'),bin.num = 25)
ctcfkoU <- ScoreMatrixBin(wt4.Ins,resize(up.boundaries,250000,fix='center'),bin.num = 25)

wtD <- ScoreMatrixBin(wt1.Ins,resize(down.boundaries,250000,fix='center'),bin.num = 25)
ctcfkoD <- ScoreMatrixBin(wt4.Ins,resize(down.boundaries,250000,fix='center'),bin.num = 25)

wtS <- ScoreMatrixBin(wt1.Ins,resize(stab.boundaries,250000,fix='center'),bin.num = 25)
ctcfkoS <- ScoreMatrixBin(wt4.Ins,resize(stab.boundaries,250000,fix='center'),bin.num = 25)

all.sml <- ScoreMatrixList(list(wtD,ctcfkoD,wtU,ctcfkoU,wtS,ctcfkoS))
plotMeta(all.sml,profile.names = c('WTd','KOd','WTu','KOu','WTs','KOs'),dispersion = 'se')

up.diff <- wtU
up.diff@.Data <- wtU@.Data - ctcfkoU@.Data

down.diff <- wtD
down.diff@.Data <- wtD@.Data - ctcfkoD@.Data

stab.diff <- wtS
stab.diff@.Data <- wtS@.Data - ctcfkoS@.Data
plotMeta(ScoreMatrixList(list(up.diff,down.diff,stab.diff)),profile.names = c('Up.diff','Down.diff','Stab.diff'),dispersion = 'se', 
           ylab='CD69-DP - CD69+SP Insulation',xlab = 'Distance from boundary (kb)',xcoords = seq(-125,125,length.out = 25))
abline(v=c(0),lty=2)

#should do the other way round- at boundaries (or I have done already, but not looked at change in rad21?) 
#not what looking for here, but upregulated rad21 tend to be in internal tads? (or at least at subtad boundaries)


difb <- '~/Documents/Projects/Thymocyte_HiC/Thesis/Differential_Insulation_boundaries/200kb_window/'
ctcf_wt_difbounds <- read.table(paste0(difb,'/CTCFvsWT.Insulation.older.bedgraph'),header = F,col.names = c('chrom','start','end','Padj','FoldChange'))
dev_difbounds <- read.table(paste0(difb,'/CD69nDPWTvsCD69pSPWT.Insulation_Boundaries.200kb_window.bedgraph'),header = F,col.names = c('chrom','start','end','Padj','FoldChange'))
dko_wt_difbounds <- read.table(paste0(difb,'/DKOvsWT.Insulation.bedgraph'),header = F,col.names = c('chrom','start','end','Padj','FoldChange'))

ctcf_wt_difbounds <- makeGRangesFromDataFrame(ctcf_wt_difbounds,keep.extra.columns = T) ; dev_difbounds <- makeGRangesFromDataFrame(dev_difbounds,keep.extra.columns = T) ; dko_wt_difbounds <- makeGRangesFromDataFrame(dko_wt_difbounds,keep.extra.columns = T)
ctcf_wt_difbounds$Padj <- ctcf_wt_difbounds$Padj * 5

plot_boundary_vs_rad21 <- function(difbounds){
  downreg.boundaries.rad21fc <- rad21.sites[unique(queryHits(findOverlaps(rad21.sites,difbounds[difbounds$Padj < 0.05 & difbounds$FoldChange < 0],maxgap = 100)))]
  upreg.boundaries.rad21fc <- rad21.sites[unique(queryHits(findOverlaps(rad21.sites,difbounds[difbounds$Padj < 0.05 & difbounds$FoldChange > 0],maxgap = 100)))]
  stable.boundaries.rad21fc <- rad21.sites[unique(queryHits(findOverlaps(rad21.sites,difbounds[difbounds$Padj > 0.5],maxgap = 100)))]
  
  upreg.boundaries.rad21fc$status <- 'Upregulated'
  downreg.boundaries.rad21fc$status <- 'Downregulated'
  stable.boundaries.rad21fc$status <- 'Stable'
  
  change.exp.rad21 <- as.data.frame(c(upreg.boundaries.rad21fc,downreg.boundaries.rad21fc,stable.boundaries.rad21fc))
  comparisons <- list( c("Stable", "Upregulated"), c("Downregulated", "Upregulated")  ,c("Downregulated", "Stable"))
  out <- ggplot(change.exp.rad21,aes(x=status,y=best.logFC , fill=status)) + geom_boxplot(outlier.alpha = 0.1) + theme_classic() + 
    scale_fill_manual(values=c('#0073C2FF','#868686FF','#EFC000FF')) + ggpubr::stat_compare_means(comparisons = comparisons) +
    theme(legend.position = 'None')
  return(out)
}

C3_4_9 <- plot_boundary_vs_rad21(dev_difbounds)
B3_4_10 <- plot_boundary_vs_rad21(ctcf_wt_difbounds)
A3_4_10 <- plot_boundary_vs_rad21(dko_wt_difbounds)
fig_3_4_10 <- plot_grid(A3_4_10,B3_4_10,nrow = 1,labels = c('A','B'))

save_plot('~/Desktop/Fig_3_4_9C.png',C3_4_9,base_asp = 1)
save_plot('~/Desktop/Fig_3_4_10AB.png',fig_3_4_10)

#change in boundary strength correlates with change in RAD21 binding for CTCFko and dev

#@ genes  #genes= resize(conds.gr.list$ctcf_ko,500,fix='center')
genes$class <- 'Stable'
genes[genes$padj < 0.05 & genes$log2FoldChange < 0]$class <- 'Downregulated'
genes[genes$padj < 0.05 & genes$log2FoldChange > 0]$class <- 'Upregulated'

rad21_file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/Rad21CD69negDPWTR2.bw'
ctcfko.rad21_file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/Rad21CD69negDPCTCFKOR1.bw'
wt4.rad21_file <- '~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/6_Rad21CD69pos4SPWTR1_mapped_sorted_RemoveDuplicates_MACS_bedGraph//Rad21CD69pos4SPWTR1.bw'

rad21_file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/bams/normFactor_normedBams/Rad21CD69negDPWTR2_mapped_sorted_RemoveDuplicates_TMM.bigwig'
wt4.rad21_file <- '~/Documents/Projects/Thymocyte_HiC/chip_seqs/bams/normFactor_normedBams/6_Rad21CD69pos4SPWTR2_mapped_sorted_RemoveDuplicates_TMM.bigwig'

bin.num=41
x <- ScoreMatrixBin(rad21_file,resize(genes[genes$class == 'Stable'],2500,fix = 'center'),bin.num = bin.num,type = 'bw' )
y <- ScoreMatrixBin(rad21_file,resize(genes[genes$class == 'Downregulated'],2500,fix = 'center'),bin.num = bin.num,type = 'bw' )
z <- ScoreMatrixBin(rad21_file,resize(genes[genes$class == 'Upregulated'],2500,fix = 'center'),bin.num = bin.num,type = 'bw' )

plotMeta(ScoreMatrixList(list(x,y,z)),winsorize = c(1,99),profile.names = c('Stable','Down','Up'),dispersion='se')
abline(v=20, lty=2)

xx <- ScoreMatrixBin(ctcfko.rad21_file,resize(genes[genes$class == 'Stable'],2500,fix = 'center'),bin.num = bin.num,type = 'bw' )
yy <- ScoreMatrixBin(ctcfko.rad21_file,resize(genes[genes$class == 'Downregulated'],2500,fix = 'center'),bin.num = bin.num,type = 'bw' )
zz <- ScoreMatrixBin(ctcfko.rad21_file,resize(genes[genes$class == 'Upregulated'],2500,fix = 'center'),bin.num = bin.num,type = 'bw' )

plotMeta(ScoreMatrixList(list(xx,yy,zz)),winsorize = c(1,99),profile.names = c('Stable','Downregulated','Upregulated'),dispersion='se')
abline(v=20, lty=2)

x.dif <- x ; x.dif@.Data <- x@.Data - xx@.Data
y.dif <- y ; y.dif@.Data <- y@.Data - yy@.Data
z.dif <- z ; z.dif@.Data <- z@.Data - zz@.Data
plotMeta(ScoreMatrixList(list(x.dif,y.dif,z.dif)),winsorize = c(1,99),profile.names = c('Stable','Downregulated','Upregulated'),dispersion='se')
abline(v=20, lty=2)

#down boundaries genes lose RAD21, but upregulated genes not gaining it after CTCFko
#could also look at logfold changes of sites overlapping tss' rather than pure change in signal

#note that data here is messy, probably want to do some kind of scaling (before splitting genes I think)
#or at least normalise to the stable genes? and get spwt bw

downreg.genes.rad21fc <- rad21.sites[unique(queryHits(findOverlaps(rad21.sites,genes[genes$class == 'Downregulated'],maxgap = 100)))]
upreg.genes.rad21fc <- rad21.sites[unique(queryHits(findOverlaps(rad21.sites,genes[genes$class == 'Upregulated'],maxgap = 100)))]
stable.genes.rad21fc <- rad21.sites[unique(queryHits(findOverlaps(rad21.sites,genes[genes$class == 'Stable'],maxgap = 100)))]

stable.genes.rad21fc <- stable.genes.rad21fc[unique(-queryHits(findOverlaps(stable.genes.rad21fc,downreg.genes.rad21fc,maxgap = 1000)))]
#stable.genes.rad21fc <- stable.genes.rad21fc[unique(-queryHits(findOverlaps(stable.genes.rad21fc,upreg.genes.rad21fc,maxgap = 1000)))]

upreg.genes.rad21fc$status <- 'Upregulated'
downreg.genes.rad21fc$status <- 'Downregulated'
stable.genes.rad21fc$status <- 'Stable'

compare= list(c('Downregulated','Stable'),c('Upregulated','Stable'))
change.exp.rad21 <- as.data.frame(c(upreg.genes.rad21fc,downreg.genes.rad21fc,stable.genes.rad21fc))
ggplot(change.exp.rad21,aes(x=status,y=best.logFC , fill=status)) + geom_boxplot(outlier.alpha = 0.01) + theme_classic() + 
  scale_fill_manual(values=c('#0073C2FF','#868686FF','#EFC000FF')) + ggpubr::stat_compare_means(comparisons = compare) +
  ylab('log Fold Change in RAD21 signal') + xlab('Gene Status DP vs. SP')

#the above boxplots also look good for dev changes. But bams for wt4 look eurgh, and not enough time to clean myself
############### now stepping back a bit- rad21 over promoters, ctcf binding sites and domains

rad21_bed <- import.bedGraph('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/Rad21_ChIPseq_in_CD69negDP/WT_CD69negDP/Rad21CD69negDPWTR1_mapped_sorted_RemoveDuplicates_peaks.subpeaks.bed')

ctcf_bed <- ctcf_bed[order(ctcf_bed$score,decreasing=T)]

all.genes.wt.rad21 <- ScoreMatrixBin(rad21_file,resize(conds.gr.list$dev,2000,fix='center'),bin.num = 101,type="bw",strand.aware = T)
all.genes.wt4.rad21 <- ScoreMatrixBin(wt4.rad21_file,resize(conds.gr.list$dev,2000,fix='center'),bin.num = 101,type="bw",strand.aware = T)

ctcfsites.wt.rad21 <- ScoreMatrixBin(rad21_file,resize(ctcf_bed,2000,fix='center'),bin.num = 101,type="bw",strand.aware = T)
ctcfsites.wt4.rad21 <- ScoreMatrixBin(wt4.rad21_file,resize(ctcf_bed,2000,fix='center'),bin.num = 101,type="bw",strand.aware = T)

#at active genes
multiHeatMatrix(ScoreMatrixList(list(all.genes.wt.rad21,all.genes.wt4.rad21)),winsorize = c(1,99),
                col=blues9,matrix.main = c('DP','SP'),order = T,
                xlab = 'Distance from promoter (bp)',xcoords = c(-1000,1000))
plotMeta(ScoreMatrixList(list(all.genes.wt.rad21,all.genes.wt4.rad21)),winsorize = c(1,99),profile.names = c('DP','SP'),
         dispersion='se',line.col = rainbow(20)[c(1,13)],dispersion.col = rainbow(20,alpha=0.25)[c(1,13)],xcoords = c(-1000,1000),
         xlab = 'Distance from promoter (bp)',ylab='RAD21 signal')
abline(v=0, lty=2)

#at all ctcf binding sites
multiHeatMatrix(ScoreMatrixList(list(ctcfsites.wt.rad21,ctcfsites.wt4.rad21)),winsorize = c(1,97.5),
                col=blues9,matrix.main = c('DP','SP'),order = T,xcoords = c(-1000,1000),common.scale = T)
plotMeta(ScoreMatrixList(list(ctcfsites.wt.rad21,ctcfsites.wt4.rad21)),winsorize = c(1,99),profile.names = c('DP','SP'),
         dispersion='se',line.col = rainbow(20)[c(1,13)],dispersion.col = rainbow(20,alpha=0.25)[c(1,13)],
         xcoords = c(-1000,1000),xlab='Distance from CTCF site (bp)')
abline(v=0, lty=2)

#####################
#####################


getlfs <- function(genes,nam){
  genes$class <- 'Stable'
  genes[genes$padj < 0.05 & genes$log2FoldChange < 0]$class <- 'Downregulated'
  genes[genes$padj < 0.05 & genes$log2FoldChange > 0]$class <- 'Upregulated'
  
  downreg.genes.rad21fc <- rad21.sites[unique(queryHits(findOverlaps(rad21.sites,genes[genes$class == 'Downregulated'],maxgap = 100)))]
  upreg.genes.rad21fc <- rad21.sites[unique(queryHits(findOverlaps(rad21.sites,genes[genes$class == 'Upregulated'],maxgap = 100)))]
  stable.genes.rad21fc <- rad21.sites[unique(queryHits(findOverlaps(rad21.sites,genes[genes$class == 'Stable'],maxgap = 100)))]
  
  stable.genes.rad21fc <- stable.genes.rad21fc[unique(-queryHits(findOverlaps(stable.genes.rad21fc,downreg.genes.rad21fc,maxgap = 1000)))]
  #stable.genes.rad21fc <- stable.genes.rad21fc[unique(-queryHits(findOverlaps(stable.genes.rad21fc,upreg.genes.rad21fc,maxgap = 1000)))]
  
  upreg.genes.rad21fc$status <- 'Upregulated'
  downreg.genes.rad21fc$status <- 'Downregulated'
  stable.genes.rad21fc$status <- 'Stable'
  
  compare= list(c('Downregulated','Stable'),c('Upregulated','Stable'))
  change.exp.rad21 <- as.data.frame(c(upreg.genes.rad21fc,downreg.genes.rad21fc,stable.genes.rad21fc))
  ggplot(change.exp.rad21,aes(x=status,y=best.logFC , fill=status)) + geom_boxplot(outlier.alpha = 0.01) + theme_classic() + 
    scale_fill_manual(values=c('#0073C2FF','#868686FF','#EFC000FF')) + ggpubr::stat_compare_means(comparisons = compare) +
    ylab('log Fold Change in RAD21 signal through Development') + xlab(paste0('Gene Status ',nam))
}

A <- getlfs(resize(conds.gr.list$ctcf_ko,500,fix='center'),'post-ΔCTCF')
B <- getlfs(resize(conds.gr.list$dko,500,fix='center'),'post-ΔCTCF/ΔRAD21')
getlfs(resize(conds.gr.list$rad21_ko,500,fix='center'),'post-ΔRAD21')

