#####
##### A quick couple of extra loop questions
##### that I don't want to further burdle up the main scripts with
##### a) putting the loopdomain/across/within barcharts into a function. so can cross compare more easily- lost vs gained vs stable
##### b) base expression and change in expression of genes with loops, no loops, lost loops gained loops
##### c) and same, but for promoter-promoter loops, and for number of loops gene has - these all just look ick :/
##### d) rewriting some of the permutation stuff just as fisher tests
##### e) size of loops to dereg genes vs stab
##### f) numbers for rad21/diff rad21 overlapping gained/lost loops/ pp loops

library(rtracklayer)
library(GenomicRanges)
library(InteractionSet)
library(ggsci)
library(ggplot2)
options(scipen=3)
#### Import gene tracks, saved as previously processed, get pvals etc

source('~/Documents/Projects/Thymocyte_HiC/Thesis/scripts/general_functions.R')
load('~/Documents/Projects/Thymocyte_HiC/Thesis/scripts/mm9_gene_expression.rdata.RData')
library(GenomicInteractions)

####
####

setwd('~/Documents/Projects/Thymocyte_HiC/Thesis/loops/')

## import loops. Careful as juicer chrom output is 1-X, rather than chr1-chrX etc

load.loops <- function(loopfile){
  x <- import(loopfile, format = "bedpe")
  first(x) <- renameSeqlevels(first(x), paste0('chr',seqlevels(first(x))))
  second(x) <- renameSeqlevels(second(x), paste0('chr',seqlevels(second(x))))
  outloops <-  makeGInteractionsFromGRangesPairs(x)
  return(outloops)
}

cd69nDPDKO.loops <- load.loops('KO_loops/KO_merged_loops_with_motifs.bedpe') 

setwd('loops_Nov18')
cd69nDPWT.loops <- load.loops('CD69nDPWT/merged_loops.bedpe') 
#cd69nDPDKO.loops <- load.loops('CD69nDPDKO/merged_loops.bedpe') 
cd69nDPCTCFKO.loops <- load.loops('CD69nDPCTCFKO/merged_loops.bedpe') 
cd69pSPWT.loops <- load.loops('CD69posCD4SPWT/merged_loops.bedpe')
#cd69pSPDKO.loops <- load.loops('CD69posCD4SPDKO/merged_loops.bedpe') 

####
setwd('../../../Thesis/loops/differential_loops/lenient')

ctcfko_ctcfko.only.loops = load.loops('hiccups_dif25_wt1_ctcf/differential_loops2.bedpe') 
ctcfko_wt.only.loops = load.loops('hiccups_dif25_wt1_ctcf/differential_loops1.bedpe') 
ctcfko_wt.both.loops = cd69nDPWT.loops[subjectHits(findOverlaps(cd69nDPCTCFKO.loops,cd69nDPWT.loops))]

nonmutual.CTCFKO <- cd69nDPCTCFKO.loops[-subjectHits(findOverlaps(cd69nDPWT.loops,cd69nDPCTCFKO.loops))]
nonmutual.WT.CTCFKO <- cd69nDPWT.loops[-queryHits(findOverlaps(cd69nDPWT.loops,cd69nDPCTCFKO.loops))]

#
dev_wt4.only.loops = load.loops('hiccups_dif25_oldWT1_wt4/differential_loops2.bedpe') 
dev_wt1.only.loops = load.loops('hiccups_dif25_oldWT1_wt4/differential_loops1.bedpe') 
dev_both.loops = cd69nDPWT.loops[subjectHits(findOverlaps(cd69pSPWT.loops,cd69nDPWT.loops))]

nonmutual.WT4 <- cd69pSPWT.loops[-subjectHits(findOverlaps(cd69nDPWT.loops,cd69pSPWT.loops))]
nonmutual.WT1 <- cd69nDPWT.loops[-queryHits(findOverlaps(cd69nDPWT.loops,cd69pSPWT.loops))]
#
dko.only.loops = load.loops('hiccups_dif25_oldWT_ko1/differential_loops2.bedpe') 
dko_wt.only.loops = load.loops('hiccups_dif25_oldWT_ko1/differential_loops1.bedpe') 
dko_both.loops = cd69nDPWT.loops[subjectHits(findOverlaps(cd69nDPDKO.loops,cd69nDPWT.loops))]


### load enhancers
cd69nDPWT.enhancers <- makeGRangesFromDataFrame(read.table('~/Documents/Projects/Thymocyte_HiC/Ya_enhancers_june18/WT_CD69nDP_enhancers.bed'),seqnames.field = 'V1',start.field = 'V2',end.field = 'V3')
cd69pSPWT.enhancers <- makeGRangesFromDataFrame(read.table('~/Documents/Projects/Thymocyte_HiC/Ya_enhancers_june18/WT_CD4SP_enhancers.bed'),seqnames.field = 'V1',start.field = 'V2',end.field = 'V3')

cd69nDPWT.superenhancers <- makeGRangesFromDataFrame(read.table('~/Documents/Projects/Thymocyte_HiC/Ya_enhancers_june18/WT_CD69nDP_superenhancers.bed'),seqnames.field = 'V1',start.field = 'V2',end.field = 'V3')
cd69pSPWT.superenhancers <- makeGRangesFromDataFrame(read.table('~/Documents/Projects/Thymocyte_HiC/Ya_enhancers_june18/WT_CD4SP_superenhancers.bed'),seqnames.field = 'V1',start.field = 'V2',end.field = 'V3')




#### (a) ###################################################################################################################
#putting the loopdomain/across/within barcharts into a function. so can cross compare more easily- lost vs gained vs stable#
############################################################################################################################

#cd69nDPDKO.loops <- load.loops('~/Documents/Projects/Thymocyte_HiC/Thesis/loops/loops_Nov18/CD69nDPDKO/merged_loops.bedpe') 

wt1.cds.gr <- read_bedpe('~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',F)
ctcfko.cds.gr <- read_bedpe('~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/CTCFko.merged.ContactDomains_keepSmall.bedpe',F)
ko1.cds.gr <- read_bedpe('~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/KO1.merged.ContactDomains_keepSmall.bedpe',F)
wt4.cds.gr <- read_bedpe('~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT4.merged.ContactDomains_keepSmall.bedpe',F)

#eurgh I might hve to set header to false from now on
#this is as a bed 
#make a loops, with each anchor being the 
dog <- function(a){return(min(a,50000))} #lol. should be 50kb 
get.loops.at.cd.corners <- function(cds,loops){
  cds.asloops <- load.loops(cds)
  cds.aspairs <- pairs(cds.asloops)
  #min doesn't work here..
  
  first(cds.aspairs) <- resize(first(cds.aspairs),width=sapply(0.5*width(first(cds.aspairs)),dog),fix ='start')
  second(cds.aspairs) <- resize(second(cds.aspairs),width=sapply(0.5*width(second(cds.aspairs)),dog),fix ='end')
  #
  z <- width(first(cds.aspairs))/2
  start(first(cds.aspairs)) <- start(first(cds.aspairs)) - z
  end(first(cds.aspairs)) <- end(first(cds.aspairs)) - z
  
  z <- width(second(cds.aspairs))/2
  start(second(cds.aspairs)) <- start(second(cds.aspairs)) + z
  end(second(cds.aspairs)) <- end(second(cds.aspairs)) + z
  
  cds.asloops <- makeGInteractionsFromGRangesPairs(cds.aspairs)
  
  loops.at.corners <- loops[queryHits(findOverlaps(loops,cds.asloops))]
  loops.away.corners <- loops[-queryHits(findOverlaps(loops,cds.asloops))]
  loopdomains <- cds.asloops[subjectHits(findOverlaps(loops,cds.asloops))]
  nonloopdomains <- cds.asloops[-subjectHits(findOverlaps(loops,cds.asloops))]
  return(list(unique(loops.at.corners),unique(loops.away.corners),unique(loopdomains),unique(nonloopdomains)))
}

loops.within.domains <- function(loops,domains){
  #assumes already filtered out loops at domain corners
  loops.df <- as.data.frame(loops)[,c('seqnames1','start1','end2')] 
  loops.gr <- makeGRangesFromDataFrame(loops.df,seqnames.field = 'seqnames1',start.field = 'start1',end.field = 'end2')
  loops.within.domains <- loops[queryHits(findOverlaps(loops.gr,domains,type = 'within'))] # i think the order should be maintained
  loops.notwithin.domains <- loops[-queryHits(findOverlaps(loops.gr,domains,type = 'within'))] # i think the order should be maintained
  return(list(unique(loops.within.domains), unique(loops.notwithin.domains)))
}  #should note that not within means not totally within

######### loops at contact domain corners (loops that make up loop domains)

types.of.loops <- function(loops,domains.txt,domains.gr){
  loops.loopdomains <- get.loops.at.cd.corners(domains.txt,loops)

  loop.ld <- length(loops.loopdomains[[1]])/length(loops)
  #
  loops.within.domains <- loops.within.domains(loops.loopdomains[[2]],domains.gr)
  #i have to be consistent here- if i use wt1.cds. gr here, I need for kos, i need to use it everywhere
  loop.wn <- length(loops.within.domains[[1]])/length(loops)
  #
  loops.across.domains <- unique(loops.within.domains[[2]][queryHits(findOverlaps(loops.within.domains[[2]],domains.gr))])

  loop.across <- length(loops.across.domains)/length(loops)
  loop.nocd <- 1 - loop.across - loop.ld - loop.wn
  return(c(loop.ld,loop.wn,loop.across,loop.nocd))
}

#####  
mutual.ctcfko.types <- types.of.loops(ctcfko_wt.both.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)
wt.only.types <- types.of.loops(ctcfko_wt.only.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)
ctcfko.only.types <- types.of.loops(ctcfko_ctcfko.only.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)

  
loop.type.distributions <- data.frame(Proportion = c(mutual.ctcfko.types,wt.only.types,ctcfko.only.types),
                                      condition=c(rep('Maintained loops',4),rep('Lost loops',4),rep('ΔCTCF gained loops',4)),
                                      type.loop = rep(c('Loop Domain','Within','Across','Non-associated'),3))

ggplot(loop.type.distributions,aes(x=condition,y=Proportion,fill=type.loop)) + geom_bar(stat='identity',position='fill') + 
  theme_classic() + ggsci::scale_fill_jco() + coord_flip() 

#####################

mutual.dev.types <- types.of.loops(dev_both.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)
wt1.only.types <- types.of.loops(dev_wt1.only.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)
wt4.only.types <- types.of.loops(dev_wt4.only.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT4.merged.ContactDomains_keepSmall.bedpe',wt4.cds.gr)


loop.type.distributions <- data.frame(Proportion = c(mutual.dev.types,wt1.only.types,wt4.only.types),
                                      condition=c(rep('Mutual Loops',4),rep('DP only',4),rep('SP only',4)),
                                      type.loop = rep(c('Loop Domain','Within','Across','Non-associated'),3))
                                      
ggplot(loop.type.distributions,aes(x=condition,y=Proportion,fill=type.loop)) + geom_bar(stat='identity',position='fill') + 
                                    theme_classic() + ggsci::scale_fill_jco() + coord_flip() 
                                      
#could also do ko

mutual.dko.types <- types.of.loops(dko_both.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)
dko.lost.types <- types.of.loops(dko_wt.only.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)


loop.type.distributions <- data.frame(Proportion = c(mutual.dko.types,dko.lost.types),
                                      condition=c(rep('Mutual Loops',4),rep('ΔCTCF/ΔRAD21 lost',4)),
                                      type.loop = rep(c('Loop Domain','Within','Across','Non-associated'),2))

ggplot(loop.type.distributions,aes(x=condition,y=Proportion,fill=type.loop)) + geom_bar(stat='identity',position='fill') + 
  theme_classic() + ggsci::scale_fill_jco() + coord_flip() 


#### (b) ###################################################################################################################
#base expression and change in expression of genes with loops, no loops, lost loops gained loops############################
############################################################################################################################

conds.gr.list$ctcf_ko <- resize(conds.gr.list$ctcf_ko, 5000, fix='center')
#pp
ctcfko.mutual.genes <- as.data.frame(conds.gr.list$ctcf_ko[unique(linkOverlaps(ctcfko_wt.both.loops,conds.gr.list$ctcf_ko)$subject1)])
ctcfko.only.genes <- as.data.frame(conds.gr.list$ctcf_ko[unique(linkOverlaps(ctcfko_ctcfko.only.loops,conds.gr.list$ctcf_ko)$subject1)])
ctcfko.WTonly.genes <- as.data.frame(conds.gr.list$ctcf_ko[unique(linkOverlaps(ctcfko_wt.only.loops,conds.gr.list$ctcf_ko)$subject1)])
noloopgenes <- conds.gr.list$ctcf_ko[-unique(linkOverlaps(cd69nDPWT.loops,conds.gr.list$ctcf_ko)$subject1)]
noloopgenes <- as.data.frame(noloopgenes[-unique(linkOverlaps(cd69nDPCTCFKO.loops,noloopgenes)$subject1)])

#p-anything
ctcfko.mutual.ep <- as.data.frame(conds.gr.list$ctcf_ko[unique(subjectHits(findOverlaps(ctcfko_wt.both.loops,conds.gr.list$ctcf_ko)))])
ctcfko.only.ep <- as.data.frame(conds.gr.list$ctcf_ko[unique(subjectHits(findOverlaps(ctcfko_ctcfko.only.loops,conds.gr.list$ctcf_ko)))])
ctcfko.WTonly.ep <- as.data.frame(conds.gr.list$ctcf_ko[unique(subjectHits(findOverlaps(ctcfko_wt.only.loops,conds.gr.list$ctcf_ko)))])
noloop.ep <- conds.gr.list$ctcf_ko[-unique(subjectHits(findOverlaps(cd69nDPWT.loops,conds.gr.list$ctcf_ko)))]
noloop.ep <- as.data.frame(noloop.ep[-unique(subjectHits(findOverlaps(cd69nDPCTCFKO.loops,noloop.ep)))])


ctcfko.mutual.genes$status <- 'Mutual Loops'
ctcfko.only.genes$status <- 'ΔCTCF Gained'
ctcfko.WTonly.genes$status <- 'Lost'
noloopgenes$status <- 'No Loops'

ctcfko.mutual.ep$status <- 'Mutual Loops'
ctcfko.only.ep$status <- 'ΔCTCF Gained'
ctcfko.WTonly.ep$status <- 'ΔCTCF Lost'
noloop.ep$status <- 'No Loops'

#test <- rbind(ctcfko.mutual.genes,ctcfko.only.genes,ctcfko.WTonly.genes, noloopgenes)
#ggplot(test,aes(x=status, y=baseMean,fill=status)) + geom_boxplot() + theme_classic() +  ggsci::scale_fill_jco() + coord_cartesian(ylim=c(0,5000))
#ggplot(test,aes(x=status, y=log2FoldChange,fill=status)) + geom_boxplot() + theme_classic() +  ggsci::scale_fill_jco() + coord_cartesian(ylim=c(-1.5,1.5))

comparisons <- list(c('No Loops','Mutual Loops'),c('ΔCTCF Lost','Mutual Loops'),c('ΔCTCF Lost','ΔCTCF Gained'),
                    c('Mutual Loops','ΔCTCF Gained'))

test <- rbind(ctcfko.mutual.ep,ctcfko.only.ep,ctcfko.WTonly.ep, noloop.ep)
test$status <- factor(test$status,levels = c('No Loops','ΔCTCF Lost','Mutual Loops','ΔCTCF Gained'))
test$manhat <- -log10(test$padj)
test$postexp <- test$baseMean * 2^test$log2FoldChange
g1 <- ggplot(test,aes(x=status, y=baseMean,fill=status)) + geom_boxplot(outlier.alpha = 0.05) + theme_classic() +  ggsci::scale_fill_jco() +
  ylab('DP WT gene expression') + xlab('Loop status')  + scale_y_continuous(trans='log10') + ggpubr::stat_compare_means(comparisons = comparisons) #+ coord_cartesian(ylim=c(0,4000))
g2 <- ggplot(test,aes(x=status, y=postexp,fill=status)) + geom_boxplot(outlier.alpha = 0.05) + theme_classic() +  ggsci::scale_fill_jco() + 
  ylab('ΔCTCF gene expression') + xlab('Loop status')  + scale_y_continuous(trans='log10') + ggpubr::stat_compare_means(comparisons = comparisons) #+ coord_cartesian(ylim=c(0,4000))
#ggplot(test,aes(x=status, y=log2FoldChange,fill=status)) + geom_boxplot() + theme_classic() +  ggsci::scale_fill_jco() + coord_cartesian(ylim=c(-1.5,1.5))
g3 <- ggplot(test,aes(x=status, y=log2FoldChange,fill=status)) + geom_violin() + theme_classic() +  ggsci::scale_fill_jco() + 
  coord_cartesian(ylim=c(-3,3)) + stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",colour = "black") +
  ylab('Log2FoldChange gene expression WT vs. ΔCTCF')  + xlab('Loop status')                          

x1<-cowplot::plot_grid(g1,g2,g3,ncol = 2)
save_plot(filename = '~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/2.13A1.png',plot = x1,base_height = 7)


#ggplot(test,aes(x=log2FoldChange,color=status)) + stat_ecdf(geom='line') + theme_classic() +  ggsci::scale_color_jco() + coord_cartesian(xlim=c(-3,5))

#it looks like there is nothing real or interesing here unfortunately. if geom_violin it, can see a v mild difference
# or.. there's a slight difference, as loop loss leads to a touch greater loss in gene expression
# but gained loops don't seem to so strongly lead to increase in expression
conds.gr.list$dev <- resize(conds.gr.list$dev, 5000, fix='center')

dev.mutual.genes <- as.data.frame(conds.gr.list$dev[unique(linkOverlaps(dev_both.loops,conds.gr.list$dev)$subject1)])
dev_wt4.only.genes <- as.data.frame(conds.gr.list$dev[unique(linkOverlaps(dev_wt4.only.loops,conds.gr.list$dev)$subject1)])
dev_wt1.only.genes <- as.data.frame(conds.gr.list$dev[unique(linkOverlaps(dev_wt1.only.loops,conds.gr.list$dev)$subject1)])
noloopgenes.dev <- conds.gr.list$dev[-unique(linkOverlaps(cd69nDPWT.loops,conds.gr.list$dev)$subject1)]
noloopgenes.dev <- as.data.frame(noloopgenes[-unique(linkOverlaps(cd69pSPWT.loops,noloopgenes.dev)$subject1),])

#p-anything
dev.mutual.ep <- as.data.frame(conds.gr.list$dev[unique(subjectHits(findOverlaps(dev_both.loops,conds.gr.list$dev)))])
dev_wt4.only.ep <- as.data.frame(conds.gr.list$dev[unique(subjectHits(findOverlaps(dev_wt4.only.loops,conds.gr.list$dev)))])
dev_wt1.only.ep <- as.data.frame(conds.gr.list$dev[unique(subjectHits(findOverlaps(dev_wt1.only.loops,conds.gr.list$dev)))])
noloop.ep.dev <- conds.gr.list$dev[-unique(subjectHits(findOverlaps(cd69nDPWT.loops,conds.gr.list$dev)))]
noloop.ep.dev <- as.data.frame(noloop.ep.dev[-unique(subjectHits(findOverlaps(cd69pSPWT.loops,noloop.ep.dev))),])
noloopgenes.dev$width <- NaN

dev.mutual.genes$status <- 'Mutual Loops'
dev_wt4.only.genes$status <- 'SP only'
dev_wt1.only.genes$status <- 'DP only'
noloopgenes.dev$status <- 'No Loops'

dev.mutual.ep$status <- 'Mutual Loops'
dev_wt4.only.ep$status <- 'SP only'
dev_wt1.only.ep$status <- 'DP only'
noloop.ep.dev$status <- 'No Loops'

#test <- rbind(dev.mutual.genes,dev_wt4.only.genes,dev_wt1.only.genes, noloopgenes.dev)
#ggplot(test,aes(x=status, y=baseMean,fill=status)) + geom_boxplot() + theme_classic() +  ggsci::scale_fill_jco() + coord_cartesian(ylim=c(0,5000))
#ggplot(test,aes(x=status, y=log2FoldChange,fill=status)) + geom_boxplot() + theme_classic() +  ggsci::scale_fill_jco() + coord_cartesian(ylim=c(-1.5,1.5))

test <- rbind(dev.mutual.ep,dev_wt4.only.ep,dev_wt1.only.ep, noloop.ep.dev)
test$status <- factor(test$status,levels = c('No Loops','DP only','Mutual Loops','SP only'))
test$postexp <- test$baseMean * 2^test$log2FoldChange
test$manhat <- -log10(test$padj)

comparisons <- list(c('No Loops','Mutual Loops'),c('DP only','Mutual Loops'),c('DP only','SP only'),
                    c('Mutual Loops','SP only'))

g1 <- ggplot(test,aes(x=status, y=baseMean,fill=status)) + geom_boxplot(outlier.alpha = 0.05) + theme_classic() +  ggsci::scale_fill_jco() +
 ylab('DP gene expression') + xlab('Loop status')  + scale_y_continuous(trans='log10') + ggpubr::stat_compare_means(comparisons = comparisons) #+ coord_cartesian(ylim=c(0,4000))
g2 <- ggplot(test,aes(x=status, y=postexp,fill=status)) + geom_boxplot(outlier.alpha = 0.05) + theme_classic() +  ggsci::scale_fill_jco() + 
  ylab('SP gene expression')  + xlab('Loop status')   + scale_y_continuous(trans='log10') + ggpubr::stat_compare_means(comparisons = comparisons) #+ coord_cartesian(ylim=c(0,4000))
#ggplot(test,aes(x=status, y=log2FoldChange,fill=status)) + geom_boxplot() + theme_classic() +  ggsci::scale_fill_jco() + coord_cartesian(ylim=c(-1.5,2.5))
g3 <- ggplot(test,aes(x=status, y=log2FoldChange,fill=status)) + geom_violin() + theme_classic() +  ggsci::scale_fill_jco() +
  coord_cartesian(ylim=c(-5,5)) + stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",colour = "black") + 
  ylab('Log2FoldChange gene expression DP vs. SP')  + xlab('Loop status')                          

x1<-cowplot::plot_grid(g1,g2,g3,ncol = 2)
save_plot(filename = '~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/2.13A2.png',plot = x1,base_height = 7)

  
#ggplot(test,aes(x=log2FoldChange,color=status)) + stat_ecdf(geom='line') + theme_classic() +  ggsci::scale_color_jco() 
######################################################################################################################################
######################################################################################################################################
conds.gr.list$dko <- resize(conds.gr.list$dko, 5000, fix='center')


dko.mutual.genes <- as.data.frame(conds.gr.list$dko[unique(linkOverlaps(dko_both.loops,conds.gr.list$dko)$subject1)])
#dko.only.genes <- as.data.frame(conds.gr.list$dko[unique(linkOverlaps(dev_wt4.only.loops,conds.gr.list$dko)$subject1)])
dko_wt.only.genes <- as.data.frame(conds.gr.list$dko[unique(linkOverlaps(dko_wt.only.loops,conds.gr.list$dko)$subject1)])
noloopgenes.dko <- conds.gr.list$dko[-unique(linkOverlaps(cd69nDPWT.loops,conds.gr.list$dko)$subject1)]
noloopgenes.dko <- as.data.frame(noloopgenes[-unique(linkOverlaps(cd69nDPDKO.loops,noloopgenes.dko)$subject1),])

#p-anything
dko.mutual.ep <- as.data.frame(conds.gr.list$dko[unique(subjectHits(findOverlaps(dko_both.loops,conds.gr.list$dko)))])
#dko.only.ep <- as.data.frame(conds.gr.list$dko[unique(subjectHits(findOverlaps(dev_wt4.only.loops,conds.gr.list$dko)))])
dko_wt.only.ep <- as.data.frame(conds.gr.list$dko[unique(subjectHits(findOverlaps(dko_wt.only.loops,conds.gr.list$dko)))])
noloop.ep.dko <- conds.gr.list$dko[-unique(subjectHits(findOverlaps(cd69nDPWT.loops,conds.gr.list$dko)))]
noloop.ep.dko <- as.data.frame(noloop.ep.dko[-unique(subjectHits(findOverlaps(cd69nDPDKO.loops,noloop.ep.dko))),])
 

dko.mutual.genes$status <- 'Mutual Loops'
#dko.only.genes$status <- 'Dev Gained'
dko_wt.only.genes$status <- 'Lost'
noloopgenes.dko$status <- 'No Loops'

dko.mutual.ep$status <- 'Mutual Loops'
#dko.only.ep$status <- 'Dev Gained'
dko_wt.only.ep$status <- 'Lost'
noloop.ep.dko$status <- 'No Loops'

#test <- rbind(dko.mutual.genes,dko_wt.only.genes, noloopgenes.dko)
#ggplot(test,aes(x=status, y=baseMean,fill=status)) + geom_boxplot() + theme_classic() +  ggsci::scale_fill_jco() + coord_cartesian(ylim=c(0,5000))
#ggplot(test,aes(x=status, y=log2FoldChange,fill=status)) + geom_boxplot() + theme_classic() +  ggsci::scale_fill_jco() + coord_cartesian(ylim=c(-1.5,1.5))

#test <- rbind(dko.mutual.ep,dko_wt.only.ep, noloop.ep.dko)
test <- rbind(dko_wt.only.ep, noloop.ep.dko)
test$status <- factor(test$status,levels = c('No Loops','Lost','Mutual Loops'))

test$manhat <- -log10(test$padj)
test$postexp <- test$baseMean * 2^test$log2FoldChange

g1 <- ggplot(test,aes(x=status, y=baseMean,fill=status)) + geom_boxplot(outlier.alpha = 0.05) + theme_classic() +  ggsci::scale_fill_jco() + 
  ylab('WT gene expression') + xlab('Loop status')  + scale_y_continuous(trans='log10') + ggpubr::stat_compare_means(comparisons = comparisons) #+ coord_cartesian(ylim=c(0,4000))
g2 <- ggplot(test,aes(x=status, y=postexp,fill=status)) + geom_boxplot(outlier.alpha = 0.05) + theme_classic() +  ggsci::scale_fill_jco() + 
  ylab('ΔCTCF/ΔRAD21 gene expression') + xlab('Loop status')  + scale_y_continuous(trans='log10') + ggpubr::stat_compare_means(comparisons = comparisons) #+ coord_cartesian(ylim=c(0,4000))

g3 <- ggplot(test,aes(x=status, y=log2FoldChange,fill=status)) + geom_violin() + theme_classic() +  ggsci::scale_fill_jco() +
  coord_cartesian(ylim=c(-3,3)) + stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",colour = "black") + 
  ylab('Log2FoldChange gene expression WT vs. ΔCTCF/ΔRAD21')  + xlab('Loop status')                          

#ggplot(test,aes(x=log2FoldChange,color=status)) + stat_ecdf(geom='line') + theme_classic() +  ggsci::scale_color_jco() + coord_cartesian(xlim=c(-5,5))

x1<-cowplot::plot_grid(g1,g2,g3,ncol = 2)
save_plot(filename = '~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/2.13A3.png',plot = x1,base_height = 7)

######

summary(pairdist(ctcfko_ctcfko.only.loops))
summary(pairdist(ctcfko_wt.only.loops))
summary(pairdist(ctcfko_wt.both.loops))

#### (a) ###################################################################################################################
#fishers tests#
############################################################################################################################


f.t.loops <- function(genes,loopsA,loopsB){
  downA <- length(subjectHits(findOverlaps(loopsA,genes[genes$padj < 0.05 & genes$log2FoldChange < 0])))
  upA <- length(subjectHits(findOverlaps(loopsA,genes[genes$padj < 0.05 & genes$log2FoldChange > 0])))
  sA <- length(subjectHits(findOverlaps(loopsA,genes[genes$padj > 0.05])))

  downB <- length(subjectHits(findOverlaps(loopsB,genes[genes$padj < 0.05 & genes$log2FoldChange < 0])))
  upB <- length(subjectHits(findOverlaps(loopsB,genes[genes$padj < 0.05 & genes$log2FoldChange > 0])))
  sB <- length(subjectHits(findOverlaps(loopsB,genes[genes$padj > 0.05])))
  
  print(fisher.test(matrix(c(downA,upA,downB,upB),nrow=2)))
  print(fisher.test(matrix(c(downA,sA,downB,sB),nrow=2)))
  print(fisher.test(matrix(c(upA,sA,upB,sB),nrow=2)))
 
   paste0('nDownA:',downA,' nDownB:',downB,' nUpA:',upA,' nupB:',upB,' nStabA:',sA,' nStabB:',sB)
  
  }

f.t.loops(conds.gr.list$ctcf_ko,nonmutual.WT.CTCFKO,nonmutual.CTCFKO)
f.t.loops(conds.gr.list$dev,nonmutual.WT1,nonmutual.WT4)
f.t.loops(conds.gr.list$dko,cd69nDPWT.loops,cd69nDPDKO.loops)

f.t.loops(conds.gr.list$ctcf_ko,ctcfko_wt.only.loops,ctcfko_ctcfko.only.loops)
f.t.loops(conds.gr.list$dev,dev_wt1.only.loops,dev_wt4.only.loops)

f.t.loops(conds.gr.list$dko,dko_wt.only.loops,cd69nDPWT.loops[-queryHits(findOverlaps(cd69nDPWT.loops,dko_wt.only.loops))])

# dev deregulated genes can be linked to differences in loops
# but others can't. in case of dko- too few loops (would need to be a slightly different question)
# for ctcfko, just seems to not be that linked? Is near significant, I think if there were more loops it would be sig.
# if increase gene size etc, actually becomes less significant so hmmm
# and up vs down (rather than up vs stable and down vs stable) is about significant
#could also use this same function to compare different types of tads from genes_vs_cds
#using fithic significant interactions, all is very stat significant but the effect size is smaller
#in both directions also

f.t.inout.loops <- function(genes,loops){
  downL <- length(genes[subjectHits(findOverlaps(loops,genes[genes$padj < 0.05 & genes$log2FoldChange < 0]))])
  upL <- length(genes[subjectHits(findOverlaps(loops,genes[genes$padj < 0.05 & genes$log2FoldChange > 0]))])
  sL <- length(genes[subjectHits(findOverlaps(loops,genes[genes$padj > 0.05]))])
  
  downO.5 <- genes[-subjectHits(findOverlaps(loops,genes[genes$padj < 0.05 & genes$log2FoldChange < 0]))]
  downO <- length(downO.5[downO.5$padj < 0.05 & downO.5$log2FoldChange < 0])
  #downO <- length(genes[-subjectHits(findOverlaps(loops,genes[genes$padj < 0.05 & genes$log2FoldChange < 0]))])
  upO.5 <- genes[-subjectHits(findOverlaps(loops,genes[genes$padj < 0.05 & genes$log2FoldChange > 0]))]
  upO <- length(upO.5[upO.5$padj < 0.05 & upO.5$log2FoldChange > 0])
  #upO <- length(genes[-subjectHits(findOverlaps(loops,genes[genes$padj < 0.05 & genes$log2FoldChange > 0]))])
  sO.5 <- genes[-subjectHits(findOverlaps(loops,genes[genes$padj > 0.05]))]
  sO <- length(sO.5[sO.5$padj > 0.05])
  
  print(fisher.test(matrix(c(downL,upL,downO,upO),nrow=2)))
  print(fisher.test(matrix(c(downL,sL,downO,sO),nrow=2)))
  print(fisher.test(matrix(c(upL,sL,upO,sO),nrow=2)))
  
  paste0('nDownIn:',downL,' nDownOut:',downO,' nUpIn:',upL,' nupOut:',upO,' nStabIn:',sL,' nStabOut:',sO)
  
}
#check how changes if unique is required for genes
f.t.inout.loops(conds.gr.list$ctcf_ko,cd69nDPWT.loops)
f.t.inout.loops(conds.gr.list$ctcf_ko,cd69nDPCTCFKO.loops)

f.t.inout.loops(conds.gr.list$dev,cd69nDPWT.loops)
f.t.inout.loops(conds.gr.list$dev,cd69pSPWT.loops)

f.t.inout.loops(conds.gr.list$dko,cd69nDPWT.loops)
f.t.inout.loops(conds.gr.list$dko,cd69nDPDKO.loops)

#in general, differentially expressed genes are linked to loops
# interestingly wt1 loops at dev reg genes only significant for up v down (but wt4 loops are significant)
# down vs stab, up vs stab are weak/ not sig (p=~0.1) for underepresentation
# im guessing because a bunch of loops are maintained, 
# and dko upreg genes aren't linked to wt loops, which also makes sense
#for ctcfko- we see upreg and downreg linked with wt and ctcfko loops, in both cases overrepresented
#with up more overrepresented than down- even for wt loops!
#i think this can then feed into the 'changes in loops strength at differentially regulated genes' stuff (ie the old boxplots)

f.t.loops.ep <- function(genes,loopsA,loopsB,mark){
  downA <- length(linkOverlaps(loopsA,genes[genes$padj < 0.05 & genes$log2FoldChange < 0], mark)$subject1)
  upA <- length(linkOverlaps(loopsA,genes[genes$padj < 0.05 & genes$log2FoldChange > 0],mark)$subject1)
  sA <- length(linkOverlaps(loopsA,genes[genes$padj > 0.05],mark)$subject1)
  
  downB <- length(linkOverlaps(loopsB,genes[genes$padj < 0.05 & genes$log2FoldChange < 0],mark)$subject1)
  upB <- length(linkOverlaps(loopsB,genes[genes$padj < 0.05 & genes$log2FoldChange > 0],mark)$subject1)
  sB <- length(linkOverlaps(loopsB,genes[genes$padj > 0.05],mark)$subject1)
  
  print(fisher.test(matrix(c(downA,upA,downB,upB),nrow=2)))
  print(fisher.test(matrix(c(downA,sA,downB,sB),nrow=2)))
  print(fisher.test(matrix(c(upA,sA,upB,sB),nrow=2)))
  paste0('nDownA:',downA,' nDownB:',downB,' nUpA:',upA,' nupB:',upB,' nStabA:',sA,' nStabB:',sB)
  
  
}

f.t.loops.ep(conds.gr.list$ctcf_ko,nonmutual.WT.CTCFKO,nonmutual.CTCFKO,cd69nDPWT.enhancers)
f.t.loops.ep(conds.gr.list$dev,nonmutual.WT1,nonmutual.WT4,cd69nDPWT.enhancers)
f.t.loops.ep(conds.gr.list$dev,nonmutual.WT1,nonmutual.WT4,cd69pSPWT.enhancers)
f.t.loops.ep(conds.gr.list$dko,cd69nDPWT.loops,cd69nDPDKO.loops,cd69nDPWT.enhancers)

f.t.loops.ep(conds.gr.list$ctcf_ko,ctcfko_wt.only.loops,ctcfko_ctcfko.only.loops,cd69nDPWT.enhancers)
f.t.loops.ep(conds.gr.list$dev,dev_wt1.only.loops,dev_wt4.only.loops,cd69nDPWT.enhancers)
f.t.loops.ep(conds.gr.list$dev,dev_wt1.only.loops,dev_wt4.only.loops,cd69pSPWT.enhancers)

#if specifically look at ep, can see a weak ctcfko gene exp loop relationship- if don't have unique. Which is odd to me
# does this mean lots of genes have multiple loops?- will be because of adjacent enhancers, hmmmm
#yep so loop w genes on one end, multiple other enhancers on other anchor. Kind of makes some sense
#but requiring 'unique' does remove the impact
#i think maybe don't discuss this
# for up as well as downregulated genes


#### (x) ###################################################################################################################
#same as above, but using fithic contacts
############################################################################################################################


wt1.sigc.l.only <-  wt.sigc.l[-queryHits(findOverlaps(wt.sigc.l,wt4.sigc.l))]
wt4.sigc.l.only <- wt4.sigc.l[-subjectHits(findOverlaps(wt.sigc.l,wt4.sigc.l))]

wt.sigc.l.only <-  wt.sigc.l[-queryHits(findOverlaps(wt.sigc.l,ctcfko.sigc.l))]
ctcfko.sigc.l.only <- ctcfko.sigc.l[-subjectHits(findOverlaps(wt.sigc.l,ctcfko.sigc.l))]

f.t.loops(conds.gr.list$ctcf_ko,wt.sigc.l.only,ctcfko.sigc.l.only)
f.t.loops(conds.gr.list$dev,wt1.sigc.l.only,wt4.sigc.l.only)

#my takes here: there are differences
#but they are not huge
#though will be slightly larger if using actual differential
#makes me want to have another go at hiccups dif
#as crap as it is eurgh
#though.... can I cheat and choose my own intermediates?

#### (x) ###################################################################################################################
#parsing hiccompare#
#ignore this, hiccompre is actually pretty crap- throws out super long distance things
############################################################################################################################

setwd('~/Documents/Projects/Thymocyte_HiC/HiC_compare/')
'WT1vsCTCF.'

hiccompare.to.loops <- function(interactions){
   anchors1u <- makeGRangesFromDataFrame(interactions[interactions$Z > 0,],seqnames.field = 'chr1',start.field = 'start1',end.field = 'end1',keep.extra.columns = F,ignore.strand = T)
   end(anchors1u) <- start(anchors1u) + resolution
   anchors2u <- makeGRangesFromDataFrame(interactions[interactions$Z > 0,],seqnames.field = 'chr2',start.field = 'start2',end.field = 'end2',keep.extra.columns = T,ignore.strand = T)
   end(anchors2u) <- start(anchors2u) + resolution
   
   anchors1d <- makeGRangesFromDataFrame(interactions[interactions$Z < 0,],seqnames.field = 'chr1',start.field = 'start1',end.field = 'end1',keep.extra.columns = F,ignore.strand = T)
   end(anchors1d) <- start(anchors1d) + resolution
   anchors2d <- makeGRangesFromDataFrame(interactions[interactions$Z < 0,],seqnames.field = 'chr2',start.field = 'start2',end.field = 'end2',keep.extra.columns = T,ignore.strand = T)
   end(anchors2d) <- start(anchors2d) + resolution
      
   return(list(GInteractions(anchors1u,anchors2u),GInteractions(anchors1d,anchors2d)))
 } 
 
parse_hiccompare <- function(comp,chroms=seq(1:19),qval=0.05,mindist=50000,maxdist=2000000){
  outup = GInteractions()
  outdown = GInteractions()
  for(i in chroms){
    x <- data.table::fread(paste0(comp,'.',as.character(i),'.hic_compare.10kb.tsv'))
    x <- x[x$distance >= mindist,]
    x <- x[x$distance <= maxdist,]
    x <- x[x$p.adj <= qval,]
    ud <- hiccompare.to.loops(x)
    outup <- c(outup,ud[[1]])
    outdown <- c(outdown,ud[[2]])
  }
  
  return(list(outup,outdown))
}

test<-parse_hiccompare('WT1vsCTCF')
test2<-parse_hiccompare('WT1vsWT4')

f.t.loops(conds.gr.list$ctcf_ko,test[[1]],test[[2]])



##As detailed in (Rao and Huntley et al., 2014),
##differential loops are only called if a loop is annotated for one file, no overlapping loop was identified in the
##other file, and the peak pixel displays less than 1.3-fold enrichment over all local neighborhoods in the other file
##^ is that it? can i just mod coolpup to get this for me
##

#####
# differential rad21 and loops. w/in and across domains.
####

rad21.ctcfko.sites <- read.table('~/Documents/Projects/Thymocyte_HiC/Differential_chip/WTvsCTCFko.RAD21.bed',sep = '\t',col.names = 
                            strsplit('chrom	start	end	width	strand	num.tests	num.up.logFC	num.down.logFC	PValue	FDR	direction	rep.test	rep.logFC	best.pos	best.logFC
','\t')[[1]])
rad21.ctcfko.sites <- makeGRangesFromDataFrame(rad21.ctcfko.sites,keep.extra.columns = T)


rad21.dev.sites <- read.table('~/Documents/Projects/Thymocyte_HiC/Differential_chip/WTvsWT4.RAD21.bed',sep = '\t',skip = 1,col.names = 
                                   strsplit('chrom	start	end	width	strand	num.tests	num.up.logF	num.down.logFC	PValue	FDR	direction	rep.test	rep.logFC	best.pos	best.logFC
','\t')[[1]])
rad21.dev.sites <- makeGRangesFromDataFrame(rad21.dev.sites,keep.extra.columns = T)

looptypes <- function(loops,domains.txt,domains.gr){
  loops.loopdomains <- get.loops.at.cd.corners(domains.txt,loops)
  
  loops.within.domains <- loops.within.domains(loops.loopdomains[[2]],domains.gr)
  #i have to be consistent here- if i use wt1.cds. gr here, I need for kos, i need to use it everywhere
  #
  loops.across.domains <- unique(loops.within.domains[[2]][queryHits(findOverlaps(loops.within.domains[[2]],domains.gr))])
  otherloops <- unique(loops.within.domains[[2]][-queryHits(findOverlaps(loops.within.domains[[2]],domains.gr))])
  
  return(list(loops.loopdomains[[1]],loops.within.domains[[1]],loops.across.domains,otherloops))
}

mutual.ctcfko.types2 <- looptypes(ctcfko_wt.both.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)
wt.only.types2 <- looptypes(ctcfko_wt.only.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)
ctcfko.only.types2 <- looptypes(ctcfko_ctcfko.only.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)

#so i guess make some barplots with these?
#withins
win.ctcfko.mutL.up <- length(mutual.ctcfko.types2[[2]][unique(queryHits(findOverlaps(mutual.ctcfko.types2[[2]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. > 0,])))])
win.ctcfko.mutL.down <- length(mutual.ctcfko.types2[[2]][unique(queryHits(findOverlaps(mutual.ctcfko.types2[[2]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. < 0,])))])
win.ctcfko.mutL.stab <- length(mutual.ctcfko.types2[[2]][unique(queryHits(findOverlaps(mutual.ctcfko.types2[[2]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR > 0.05 ,])))])

win.ctcfko.gainL.up <- length(ctcfko.only.types2[[2]][unique(queryHits(findOverlaps(ctcfko.only.types2[[2]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. > 0,])))])
win.ctcfko.gainL.down <- length(ctcfko.only.types2[[2]][unique(queryHits(findOverlaps(ctcfko.only.types2[[2]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. < 0,])))])
win.ctcfko.gainL.stab <- length(ctcfko.only.types2[[2]][unique(queryHits(findOverlaps(ctcfko.only.types2[[2]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR > 0.05 ,])))])

win.ctcfko.lossL.up <- length(wt.only.types2[[2]][unique(queryHits(findOverlaps(wt.only.types2[[2]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. > 0,])))])
win.ctcfko.lossL.down <- length(wt.only.types2[[2]][unique(queryHits(findOverlaps(wt.only.types2[[2]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. < 0,])))])
win.ctcfko.lossL.stab <- length(wt.only.types2[[2]][unique(queryHits(findOverlaps(wt.only.types2[[2]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR > 0.05 ,])))])
#


#across
across.ctcfko.mutL.up <- length(mutual.ctcfko.types2[[3]][unique(queryHits(findOverlaps(mutual.ctcfko.types2[[3]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. > 0,])))])
across.ctcfko.mutL.down <-length(mutual.ctcfko.types2[[3]][unique(queryHits(findOverlaps(mutual.ctcfko.types2[[3]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. < 0,])))])
across.ctcfko.mutL.stab <-length(mutual.ctcfko.types2[[3]][unique(queryHits(findOverlaps(mutual.ctcfko.types2[[3]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR > 0.05 ,])))])

across.ctcfko.gainL.up <- length(ctcfko.only.types2[[3]][unique(queryHits(findOverlaps(ctcfko.only.types2[[3]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. > 0,])))])
across.ctcfko.gainL.down <- length(ctcfko.only.types2[[3]][unique(queryHits(findOverlaps(ctcfko.only.types2[[3]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. < 0,])))])
across.ctcfko.gainL.stab <- length(ctcfko.only.types2[[3]][unique(queryHits(findOverlaps(ctcfko.only.types2[[3]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR > 0.05 ,])))])

across.ctcfko.lossL.up <- length(wt.only.types2[[3]][unique(queryHits(findOverlaps(wt.only.types2[[3]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. > 0,])))])
across.ctcfko.lossL.down <- length(wt.only.types2[[3]][unique(queryHits(findOverlaps(wt.only.types2[[3]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR < 0.05 & rad21.ctcfko.sites$best.logFC. < 0,])))])
across.ctcfko.lossL.stab <- length(wt.only.types2[[3]][unique(queryHits(findOverlaps(wt.only.types2[[3]],
                                  rad21.ctcfko.sites[rad21.ctcfko.sites$FDR > 0.05 ,])))])

#and for dev

mutual.dev.types2 <- looptypes(dev_both.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)
wt1.only.types2 <- looptypes(dev_wt1.only.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)
wt4.only.types2 <- looptypes(dev_wt4.only.loops,'~/Documents/Projects/Thymocyte_HiC/Ya_small_domains/WT1.merged.ContactDomains_keepSmall.bedpe',wt1.cds.gr)
#
win.dev.mutL.up <- length(mutual.dev.types2[[2]][unique(queryHits(findOverlaps(mutual.dev.types2[[2]],
                        rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. > 0,])))])
win.dev.mutL.down <- length(mutual.dev.types2[[2]][unique(queryHits(findOverlaps(mutual.dev.types2[[2]],
                        rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. < 0,])))])
win.dev.mutL.stab <- length(mutual.dev.types2[[2]][unique(queryHits(findOverlaps(mutual.dev.types2[[2]],
                        rad21.dev.sites[mutual.dev.types2$FDR > 0.05 ,])))])

win.dev.gainL.up <- length(wt4.only.types2[[2]][unique(queryHits(findOverlaps(wt4.only.types2[[2]],
                          rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. > 0,])))])
win.dev.gainL.down <- length(wt4.only.types2[[2]][unique(queryHits(findOverlaps(wt4.only.types2[[2]],
                          rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. < 0,])))])
win.dev.gainL.stab <- length(wt4.only.types2[[2]][unique(queryHits(findOverlaps(wt4.only.types2[[2]],
                          rad21.dev.sites[rad21.dev.sites$FDR > 0.05 ,])))])

win.dev.lossL.up <- length(wt1.only.types2[[2]][unique(queryHits(findOverlaps(wt1.only.types2[[2]],
                         rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. > 0,])))])
win.dev.lossL.down <- length(wt1.only.types2[[2]][unique(queryHits(findOverlaps(wt1.only.types2[[2]],
                         rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. < 0,])))])
win.dev.lossL.stab <- length(wt1.only.types2[[2]][unique(queryHits(findOverlaps(wt1.only.types2[[2]],
                        rad21.dev.sites[rad21.dev.sites$FDR > 0.05 ,])))])
#


#across
across.dev.mutL.up <- length(mutual.dev.types2[[3]][unique(queryHits(findOverlaps(mutual.dev.types2[[3]],
                          rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. > 0,])))])
across.dev.mutL.down <-length(mutual.dev.types2[[3]][unique(queryHits(findOverlaps(mutual.dev.types2[[3]],
                          rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. < 0,])))])
across.dev.mutL.stab <-length(mutual.dev.types2[[3]][unique(queryHits(findOverlaps(mutual.dev.types2[[3]],
                          rad21.dev.sites[rad21.dev.sites$FDR > 0.05 ,])))])

across.dev.gainL.up <- length(wt4.only.types2[[3]][unique(queryHits(findOverlaps(wt4.only.types2[[3]],
                           rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. > 0,])))])
across.dev.gainL.down <- length(wt4.only.types2[[3]][unique(queryHits(findOverlaps(wt4.only.types2[[3]],
                           rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. < 0,])))])
across.dev.gainL.stab <- length(wt4.only.types2[[3]][unique(queryHits(findOverlaps(wt4.only.types2[[3]],
                            rad21.dev.sites[rad21.dev.sites$FDR > 0.05 ,])))])

across.dev.lossL.up <- length(wt1.only.types2[[3]][unique(queryHits(findOverlaps(wt1.only.types2[[3]],
                            rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. > 0,])))])
across.dev.lossL.down <- length(wt1.only.types2[[3]][unique(queryHits(findOverlaps(wt1.only.types2[[3]],
                            rad21.dev.sites[rad21.dev.sites$FDR < 0.05 & rad21.dev.sites$best.logFC. < 0,])))])
across.dev.lossL.stab <- length(wt1.only.types2[[3]][unique(queryHits(findOverlaps(wt1.only.types2[[3]],
                            rad21.dev.sites[rad21.dev.sites$FDR > 0.05 ,])))])






#################hand stiched bar charts

ctcko.lpfate <- data.frame(condition=c('WT','WT','WT','WT','ΔCTCF','ΔCTCF','ΔCTCF','ΔCTCF'),
           type.loop=c('Promoter-Promoter','Promoter-Other','Loop-domain','Other'),
           Freq=c(18,107,94,237,20,45,2,16))

dev.lpfate <- data.frame(condition=c('DP','DP','DP','DP','SP','SP','SP','SP'),
                           type.loop=c('Promoter-Promoter','Promoter-Other','Loop-domain','Other'),
                           Freq=c(38,109,34,38,3,20,7,36))


ggplot(ctcko.lpfate,aes(x=condition,y=Freq,fill=type.loop)) + geom_bar(stat='identity',position='fill') + 
  theme_classic() + ggsci::scale_fill_jco() + coord_flip() 

ggplot(dev.lpfate,aes(x=condition,y=Freq,fill=type.loop)) + geom_bar(stat='identity',position='fill') + 
  theme_classic() + ggsci::scale_fill_jco() + coord_flip() 


#
ggplot(dev.lpfate[dev.lpfate$type.loop != 'Other',],aes(x=condition,y=Freq,fill=type.loop)) + geom_bar(stat='identity',position='fill') + 
  theme_classic() + ggsci::scale_fill_jco() + coord_flip() 
ggplot(ctcko.lpfate[ctcko.lpfate$type.loop != 'Other',],aes(x=condition,y=Freq,fill=type.loop)) + geom_bar(stat='identity',position='fill') +
  theme_classic() + ggsci::scale_fill_jco() + coord_flip() 
