'''
Just using Yas analysed RNA seq data
should check ,liz filters for gene expression
because I want to know why our dev vs rad21 ko logfold vs logfold plots look so different

should add a "dev change on CTCF vs RAD21 ko graph" 
should colour just 

note that in particular, upregulated in dev genes *especially* tend not to be downregulated in the KOs

do I have a "gene expression at differential boundaries" thing anywhere?
sooo. its laziness, but just run this 3 times with/ without / 

filtered is for whether the barplots and fishers tests look at
all genes / genes without RAD21/CTCF binding nearby / genes *with* proximal CTCF or RAD21)
if filtered is on, the .X refers to genes that do not have overlapping rad21/ctcf
other do
'''

source('~/Documents/Projects/Thymocyte_HiC/Thesis/scripts/general_functions.R')

library(ggplot2)
library(RColorBrewer)
library(regioneR)
library(cowplot)
library(ggsci)
library(ggpubr)

setwd('~/Documents/Projects/Thymocyte_HiC/DEseq2_in_vivo_Development_YaRNAseq/')

theme_set(theme_classic(base_size = 6))

filtered = T

#wtko <- read.csv('01_CD69nDPWT1-3_vs_CD69nDPDKO1-3.csv',header=T)
#wtwt <- read.csv('04_CD69nDPWT1-3_vs_CD69pCD4SPWT1-2.csv',header=T)
wtko <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_DKO_CD69nDPWT1-3_vs_CD69nDPDKO1-3.csv',header=T)
wtko <- wtko[!is.na(wtko$padj),]
#wtko$padj[is.na(wtko$padj)] <- 1 #none of the pvals that don't get converted to padj (why???) are strong enough to be actual dif genes

wtwt <- read.csv('04_CD69nDPWT1-3_vs_CD69pCD4SPWT1-2.csv',header=T)
wtwt <- wtwt[!is.na(wtwt$padj),]  #genes that are na here are na as their basemean is super low
#wtwt$padj[is.na(wtwt$padj)] <- 1 #none of the pvals that don't get converted to padj (why???) are strong enough to be actual dif genes


wtctcf  <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_CTCFKO_CD69nDPWT4-5_vs_CD69nDPCTCFKO1-2.csv',header=T)
wtctcf <- wtctcf[!is.na(wtctcf$padj),]
#wtctcf$padj[is.na(wtctcf$padj)] <- 1 #none of the pvals that don't get converted to padj (why???) are strong enough to be actual dif genes


wtrad21 <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_Rad21KO_DPWT1-2_vs_DPRad21KO1-2.csv',header=T)
wtrad21 <- wtrad21[!is.na(wtrad21$padj),]
#wtrad21$padj[is.na(wtrad21$padj)] <- 1 #none of the pvals that don't get converted to padj (why???) are strong enough to be actual dif genes

#remove duplicates
wtko <- wtko[order(wtko$ensembl_gene_id, abs(wtko$pvalue) ), ]
wtko <-  wtko[ !duplicated(wtko$ensembl_gene_id), ] 

wtwt <- wtwt[order(wtwt$ensembl_gene_id, abs(wtwt$pvalue) ), ]
wtwt <-  wtwt[ !duplicated(wtwt$ensembl_gene_id), ] 

wtctcf <- wtctcf[order(wtctcf$ensembl_gene_id, abs(wtctcf$pvalue) ), ]
wtctcf <-  wtctcf[ !duplicated(wtctcf$ensembl_gene_id), ] 

wtrad21 <- wtrad21[order(wtrad21$ensembl_gene_id, abs(wtrad21$pvalue) ), ]
wtrad21 <-  wtrad21[ !duplicated(wtrad21$ensembl_gene_id), ] 
#merge


#fold_changes <- merge(wtko,wtwt,by='ensembl_gene_id')
#fold_changes[is.na(fold_changes)] <- 1 #remove the padj which go to NA, upset stuff
#fold_changes$status <- 'unchanged'
#fold_changes[fold_changes$padj.x <= 0.05,]$status <- 'Cohesin/CTCF KO'
#fold_changes[fold_changes$padj.y <= 0.05,]$status <- 'Development'
#fold_changes[fold_changes$padj.y <= 0.05 & fold_changes$padj.x <= 0.05,]$status <- 'Both'
# looks better if theres a cutoff on the fold changes logfold

merge_diffgene_exps <- function(diffgene1,diffgene2,pvalcut1,pvalcut2,logfoldcut1,logfoldcut2,cond1,cond2){
  merged_genes <- merge(diffgene1,diffgene2,by='ensembl_gene_id')
  #fold_changes[is.na(fold_changes)] <- 1 #remove the padj which go to NA, upset stuff
  merged_genes$status <- 'unchanged'
  merged_genes[merged_genes$padj.x <= pvalcut1 & abs(merged_genes$log2FoldChange.x) >logfoldcut1,]$status <- cond1
  merged_genes[merged_genes$padj.y <= pvalcut2 & abs(merged_genes$log2FoldChange.y) >logfoldcut2,]$status <- cond2
  merged_genes[merged_genes$padj.y <= pvalcut2 & merged_genes$padj.x <= pvalcut1 & abs(merged_genes$log2FoldChange.x) >logfoldcut1 & abs(merged_genes$log2FoldChange.y) >logfoldcut2,]$status <- 'Both'
  return(merged_genes)
}

dev_vs_dko <- merge_diffgene_exps(wtwt,wtko,0.05,0.05,1.5,1.5,'Development','ΔCTCF/ΔRAD21')
dev_vs_ctcf <- merge_diffgene_exps(wtwt,wtctcf,0.05,0.05,1.5,1.5,'Development','ΔCTCF')
dev_vs_rad21 <- merge_diffgene_exps(wtwt,wtrad21,0.05,0.05,1.5,1.1,'Development','ΔRAD21')
ctcf_vs_dko <- merge_diffgene_exps(wtctcf,wtko,0.05,0.05,1.5,1.5,'ΔCTCF','ΔCTCF/ΔRAD21')
ctcf_vs_rad21 <- merge_diffgene_exps(wtctcf,wtrad21,0.05,0.05,1.5,1.1,'ΔCTCF','ΔRAD21')
rad21_vs_dko <- merge_diffgene_exps(wtrad21,wtko,0.05,0.05,1.1,1.5,'ΔRAD21','ΔCTCF/ΔRAD21')

fc.data <- dev_vs_ctcf
f3.dev_vs_ctcf.scatter <- ggplot() + geom_point(data = fc.data[fc.data$status== 'unchanged',], aes(x=log2FoldChange.x,y=log2FoldChange.y),size=0.25,alpha=0.5) + 
  geom_point(data=fc.data[fc.data$status != 'unchanged',] , aes(x=log2FoldChange.x,y=log2FoldChange.y,color=status),size=0.5,alpha=0.7) + 
  scale_colour_manual(values = c('red3','blue3',"dodgerblue")) + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + 
  geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) + 
  xlab('Log2 Fold Change CD69-DP vs CD69+SP') + ylab('Log2 Fold Change CD69-DP vs ΔCTCF') +
  theme(legend.position = 'bottom',legend.spacing.x = unit(0,'cm'))

fc.data <- dev_vs_rad21
f3.dev_vs_rad21.scatter <- ggplot() + geom_point(data = fc.data[fc.data$status== 'unchanged',], aes(x=log2FoldChange.x,y=log2FoldChange.y),size=0.25,alpha=0.5) +  geom_point(data=fc.data[fc.data$status != 'unchanged',] , aes(x=log2FoldChange.x,y=log2FoldChange.y,color=status),size=0.5,alpha=0.7) + 
  scale_colour_manual(values = c('red3','blue3',"dodgerblue")) + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + 
  geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) + 
  xlab('Log2 Fold Change CD69-DP vs CD69+SP') + ylab('Log2 Fold Change CD69-DP vs ΔRAD21') +
  theme(legend.position = 'bottom',legend.spacing.x = unit(0,'cm'))


fc.data <- ctcf_vs_dko
f2.ctcf_vs_dko.scatter <- ggplot() + geom_point(data = fc.data[fc.data$status== 'unchanged',], aes(x=log2FoldChange.x,y=log2FoldChange.y),size=0.25,alpha=0.5) +  geom_point(data=fc.data[fc.data$status != 'unchanged',] , aes(x=log2FoldChange.x,y=log2FoldChange.y,color=status),size=0.5,alpha=0.7) + scale_colour_manual(values = c('red3','blue3',"dodgerblue")) + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) + 
  xlab('Log2 Fold Change CD69-DP vs ΔCTCF') + ylab('Log2 Fold Change CD69-DP vs ΔCTCF/ΔRAD21') +
  theme(legend.position = 'bottom',legend.spacing.x = unit(0,'cm'))

fc.data <- ctcf_vs_rad21
f2.ctcf_vs_rad21.scatter <- ggplot() + geom_point(data = fc.data[fc.data$status== 'unchanged',], aes(x=log2FoldChange.x,y=log2FoldChange.y),size=0.25,alpha=0.5) +  geom_point(data=fc.data[fc.data$status != 'unchanged',] , aes(x=log2FoldChange.x,y=log2FoldChange.y,color=status),size=0.5,alpha=0.7) + scale_colour_manual(values = c('red3','blue3',"dodgerblue")) + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) + 
  xlab('Log2 Fold Change CD69-DP vs ΔCTCF') + ylab('Log2 Fold Change CD69-DP vs ΔRAD21') +
  theme(legend.position = 'bottom',legend.spacing.x = unit(0,'cm'))

fc.data <- rad21_vs_dko
f2.rad21_vs_dko.scatter <- ggplot() +  geom_point(data = fc.data[fc.data$status== 'unchanged',], aes(x=log2FoldChange.x,y=log2FoldChange.y),size=0.25,alpha=0.5) +  geom_point(data=fc.data[fc.data$status != 'unchanged',] , aes(x=log2FoldChange.x,y=log2FoldChange.y,color=status),size=0.5,alpha=0.7)  + scale_colour_manual(values = c('red3','blue3',"dodgerblue")) + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) + 
  xlab('Log2 Fold Change CD69-DP vs ΔRAD21') + ylab('Log2 Fold Change CD69-DP vs ΔCTCF/ΔRAD21') +
  theme(legend.position = 'bottom',legend.spacing.x = unit(0,'cm'))

f2.scatter.legend <- get_legend(ggplot() +  geom_point(data = fc.data[fc.data$status== 'unchanged',], aes(x=log2FoldChange.x,y=log2FoldChange.y),size=0.25,alpha=0.5) +  geom_point(data=fc.data[fc.data$status != 'unchanged',] , aes(x=log2FoldChange.x,y=log2FoldChange.y,color=status),size=0.5,alpha=0.7) + scale_colour_manual(values = c('red3','blue3',"dodgerblue")) + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) + 
                                  xlab('Log2 Fold Change CD69-DP vs CD69+SP') + ylab('Log2 Fold Change CD69-DP vs ΔCTCF/ΔRAD21'))

#f2.scatters <- plot_grid(f2.ctcf_vs_rad21.scatter,f2.ctcf_vs_dko.scatter,f2.rad21_vs_dko.scatter,f2.scatter.legend,
#                         nrow = 2,ncol=3,labels = c('D','E','F'),label_size = 9,hjust = 0)

f2.scatters <- plot_grid(f2.ctcf_vs_rad21.scatter,f2.ctcf_vs_dko.scatter,f2.rad21_vs_dko.scatter,
                         nrow = 1,labels = c('D','E','F'),label_size = 9,hjust = 0)

ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/kovsko.scatters.pdf',
        f2.scatters,width=7,height=3,units = 'in',device = cairo_pdf)

#height*0.5 should equal ~1/3 with

#ggplot(fold_changes,aes(x=log2FoldChange.x,y=log2FoldChange.y,color=status)) + geom_point(size=0.25,alpha=0.5) + theme_classic() + scale_colour_manual(values = c('red4','dodgerblue','blue3','grey'))
##+ xlim(c(-10,10)) + ylim(c(-10,10))
#+ geom_hline(yintercept = 0,linetype='dashed',size=0.25) + geom_vline(xintercept = 0,linetype='dashed',size=0.25)
##actually lets just do black red blue, make dif ex dots larger (like plotMA etc)
##size looks better larger. Should plot differnetial point larger and on top


fc.data <- dev_vs_dko
f4.dev_vs_dko.scatter <- ggplot() + geom_point(data = fc.data[fc.data$status== 'unchanged',], aes(x=log2FoldChange.x,y=log2FoldChange.y),size=0.25,alpha=0.5) +  geom_point(data=fc.data[fc.data$status != 'unchanged',] , aes(x=log2FoldChange.x,y=log2FoldChange.y,color=status),size=0.5,alpha=0.7) + 
  scale_colour_manual(values = c('red3','blue3',"dodgerblue")) + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + 
  geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) + 
  xlab('Log2 Fold Change CD69-DP vs CD69+SP') + ylab('Log2 Fold Change CD69-DP vs CD69-DP ΔCTCF/ΔRAD21') +
  theme(legend.position = 'bottom',legend.spacing.x = unit(0,'cm'))

f3.scatters <- plot_grid(f3.dev_vs_ctcf.scatter,f3.dev_vs_rad21.scatter,f4.dev_vs_dko.scatter,
                         nrow = 1,labels = c('A','B','C'),label_size = 9,hjust = 0)

ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/devvsko.scatters.pdf',
        f3.scatters,width=7,height=3,units = 'in',device = cairo_pdf)


#we see here that gene upreg in dev are very rarely downregulated by dko
# we also see that the dko correlates much better to the rad21ko than the ctcfko
# and that the similar change in direction vs dev comes more from rad21 than ctcf ko

#### same as above, but removing genes overlapping rad21 or ctcf binding sites
#requires fc.data to be using dev_vs_dko, which I'm fine with
fc.data$chrom <- paste0('chr',fc.data$chromosome_name)
fc.data$strand2 <- '*'
fc.data$strand2[fc.data$strand == -1] <- '-'
fc.data$strand2[fc.data$strand == 1] <- '+'
fc.data$strand <- fc.data$strand2
fc.data.gr <- makeGRangesFromDataFrame(fc.data,seqnames.field = 'chrom',start.field = 'start_position',
                                       strand.field = 'strand',end.field = 'end_position',keep.extra.columns = T)
fc.data.gr <- promoters(fc.data.gr)
ctcf.sites.gr <- rtracklayer::import.bed('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/1_CTCFCD69negDPWTR1/1_CTCFCD69negDPWTR1_peaks.subpeaks.bed')
rad21.sites.gr <- rtracklayer::import.bed('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/Rad21_ChIPseq_in_CD69negDP/WT_CD69negDP/Rad21CD69negDPWTR1_mapped_sorted_RemoveDuplicates_peaks.subpeaks.bed')

fc.data.gr.f <- unique(fc.data.gr[-queryHits(findOverlaps(fc.data.gr,ctcf.sites.gr,maxgap = 500))])
fc.data.gr.filtered <- fc.data.gr.f[-queryHits(findOverlaps(fc.data.gr.f,rad21.sites.gr,maxgap = 500))]

fc.data.ctcfrad21.filtered <- as.data.frame(fc.data.gr.filtered)

f3.rad21ctcf_filtered.scatter <- ggplot() + geom_point(data = fc.data.ctcfrad21.filtered[fc.data.ctcfrad21.filtered$status== 'unchanged',], 
  aes(x=log2FoldChange.x,y=log2FoldChange.y),size=0.25,alpha=0.5) + 
  geom_point(data=fc.data.ctcfrad21.filtered[fc.data.ctcfrad21.filtered$status != 'unchanged',],
  aes(x=log2FoldChange.x,y=log2FoldChange.y,color=status),size=0.5,alpha=0.7) + 
  scale_colour_manual(values = c('red3','blue3',"dodgerblue")) + 
  geom_hline(yintercept = 0,linetype='dashed',size=0.25) + 
  geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) + 
  xlab('Log2 Fold Change CD69-DP vs CD69+SP') + ylab('Log2 Fold Change CD69-DP vs CD69-DP ΔCTCF/ΔRAD21') +
  theme(legend.position = 'bottom',legend.spacing.x = unit(0,'cm'))

#Δ
#see that it isn't dependent 
#and now the same again, but color on whether point has ctcf/rad21 or not?
fc.data.gr$mark <- 'None'
fc.data.gr[queryHits(findOverlaps(fc.data.gr,ctcf.sites.gr,maxgap = 500))]$mark <- 'CTCF'
fc.data.gr[queryHits(findOverlaps(fc.data.gr,rad21.sites.gr,maxgap = 500))]$mark <- 'RAD21'

both.gr <- ctcf.sites.gr[queryHits(findOverlaps(ctcf.sites.gr,rad21.sites.gr,maxgap = 250))]
fc.data.gr[queryHits(findOverlaps(fc.data.gr,both.gr))]$mark <- 'Both'

fc.data.m <- as.data.frame(fc.data.gr)

"""ggplot(fc.data.m,aes(x=log2FoldChange.x,y=log2FoldChange.y,color=mark)) + geom_point(size=0.25,alpha=0.5) +
  theme_classic() + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + 
  geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) + 
  xlab('Log2 Fold Change CD69-DP vs CD69+SP') + ylab('Log2 Fold Change CD69-DP vs CD69-DP ΔCTCF/ΔRAD21')


ggplot() + geom_point(data=fc.data.m[fc.data.m$mark == 'None',],aes(x=log2FoldChange.x,y=log2FoldChange.y),size=0.25,alpha=0.5) +
  geom_point(data=fc.data.m[fc.data.m$mark != 'None',],
             aes(x=log2FoldChange.x,y=log2FoldChange.y,color=mark),size=0.25,alpha=0.5) + 
   scale_colour_manual(values = c('red3','blue3',"dodgerblue")) + 
  theme_classic() + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + 
  geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) + 
  xlab('Log2 Fold Change CD69-DP vs CD69+SP') + ylab('Log2 Fold Change CD69-DP vs CD69-DP ΔCTCF/ΔRAD21')
"""

f4.dev_vs_dko.scatter.split <- ggplot(fc.data.m,aes(x=log2FoldChange.x,y=log2FoldChange.y,color=mark)) + 
  geom_point(size=0.25,alpha=0.5) + facet_wrap(.~mark) + scale_color_jco() + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + 
  geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-5,5),ylim = c(-5,5)) +  
  geom_density_2d(color='black') +
  xlab('Log2 Fold Change CD69-DP vs CD69+SP') + ylab('Log2 Fold Change CD69-DP vs CD69-DP ΔCTCF/ΔRAD21')
#i think this works not so well as a scatter plot(in this specific instance)
#put 1d density plots over like seaborn scatters #could highlight difreg genes?
#ggplot(fc.data.m,aes(x=log2FoldChange.y,fill=mark)) + stat_density(alpha=0.5,color='black')  + theme_classic() + scale_fill_jco() +   geom_hline(yintercept = 0,linetype='dashed',size=0.25) + 
#  geom_vline(xintercept = 0,linetype='dashed',size=0.25)  + xlab('Log2 Fold Change CD69-DP vs CD69+SP') + ylab('Log2 Fold Change CD69-DP vs CD69-DP ΔCTCF/ΔRAD21') + coord_cartesian(xlim=c(-5,5))

f4.dp <- ggplot(fc.data.m,aes(x=log2FoldChange.x,color=mark)) + geom_density(alpha=0.5)  + theme_classic() + scale_color_jco() +   geom_hline(yintercept = 0,linetype='dashed',size=0.25) + 
  geom_vline(xintercept = 0,linetype='dashed',size=0.25)  + xlab('Log2 Fold Change CD69-DP vs CD69+SP') + coord_cartesian(xlim=c(-5,5)) ######get gene locations for regioneR for perm tests

#not a fundemental shift in fold change based on whether there is rad21 or ctcf binding at gene or not, although
#genes with rad21 at promoters actually slightly more centerally distributed?
#in all cases see a minor correlation between dko and dev fold change (from ~0.12-0.26)
#takeaway: even though rad21/ctcf binding influences gene deregulation, certainly not the 'cause' of similarity to dev genes
# and not a major explanation of dev or dko gene fluctuations
# ie a subset of genes are rad21/ ctcf sensitive, but not true for wider set of genes that has rad21/ctcf binding
#in dko as well as dev

mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl',
               host="may2012.archive.ensembl.org")

mm9_genes<-getBM(attributes=c("ensembl_gene_id",
                              "external_gene_id",
                              "chromosome_name",
                              "start_position",
                              "end_position",
                              "gene_biotype"),
                 mart=mart)

mm9_pc_genes<-mm9_genes[mm9_genes$gene_biotype=="protein_coding",]

mm9_pc_genes.gr<-GRanges(seqnames=mm9_pc_genes$chromosome_name,
                         ranges=IRanges(start=mm9_pc_genes$start_position,
                                        end=mm9_pc_genes$end_position),
                         ens_id=mm9_pc_genes$ensembl_gene_id,
                         name=mm9_pc_genes$external_gene_id)

seqlevelsStyle(mm9_pc_genes.gr)<-"UCSC"

##### perm tests for 
baseMean_filtered <- wtwt[wtwt$baseMean > 20,]
all_exp_genes.gr <- mm9_pc_genes.gr[mm9_pc_genes.gr$ens_id %in% baseMean_filtered$ensembl_gene_id] 
#^ might want to set a baseMean cutoff here <- yep improves things. Also assumes same genes in both sets, should be a collected univers
#as the distributions of expression values for developmental shifting are quite dif
# interesting, developmental changing genes tend to have a much higher baseMean than random
#expressed genes, whilst cohesin KO changed genes don't really
#so I feel less bad not accounting for it, but best to have a slightly higher basemean
#things like gene lengths might still be having an effect- bit longer in both KO and dev changed genes
#but some quick shuffling of my own on results not expected to be significant has allayed those worries
#also, maaaybe of interest-- the WTdev genes that really really don't change are sig underrepresented in KO (as in, WTdev pval > 0.8)
#also genes upregulated through development and *downregulated* with ctcf/rad21 ko are actually underrepresented!!
#and for wtvsrad21 vs wtvsctcf, the signal is *much* stronger for both downregulated than both up
#ctcf up doesn't match well super well any of the rad21 combs hmmm
#results stand (stronger!) when comparing high logFold changes not pvals, which is nice

set1_diffexp <- wtwt[wtwt$padj < 0.05,] 
set1_diffexp.gr <- mm9_pc_genes.gr[mm9_pc_genes.gr$ens_id %in% set1_diffexp$ensembl_gene_id]

set2_diffexp <- wtko[wtko$padj < 0.05,] 
set2_diffexp.gr <- mm9_pc_genes.gr[mm9_pc_genes.gr$ens_id %in% set2_diffexp$ensembl_gene_id]

pt <- permTest(A=set1_diffexp.gr, ntimes=100, randomize.function=resampleRegions, universe=all_exp_genes.gr,
                        evaluate.function=numOverlaps, B=set2_diffexp.gr, verbose=T)



#####
##### log(count)-log(count) plots
fold_change <- 0
ctcf_difexp <- wtctcf[wtctcf$padj < 0.05 & abs(wtctcf$log2FoldChange) > fold_change,]$ensembl_gene_id
ctcf_counts <- read.table('~/Documents/Projects/Thymocyte_HiC/RNA_seq_counts/RAWCOUNTS_CTCFKO.tsv',header=T)
ctcf_counts$status <- 'unchanged'
ctcf_counts[ctcf_counts$id %in% ctcf_difexp ,]$status <- 'Diff'
ctcf_counts$wt <- (ctcf_counts$RNACD69negDPWTR4 + ctcf_counts$RNACD69negDPWTR5) /2
ctcf_counts$ko <- (ctcf_counts$RNACD69negDPCTCFKOR1 + ctcf_counts$RNACD69negDPCTCFKOR2) /2

rad21_difexp <- wtrad21[wtrad21$padj < 0.05 & abs(wtrad21$log2FoldChange) > fold_change,]$ensembl_gene_id
rad21_counts <- read.table('~/Documents/Projects/Thymocyte_HiC/RNA_seq_counts/RAWCOUNTS_Rad21KO.tsv',header=T)
rad21_counts$status <- 'unchanged'
rad21_counts[rad21_counts$id %in% rad21_difexp ,]$status <- 'Diff'
rad21_counts$wt <- (rad21_counts$WT1 + rad21_counts$WT2) /2
rad21_counts$ko <- (rad21_counts$KO1 + rad21_counts$KO2) /2

dko_difexp <- wtko[wtko$padj < 0.05 & abs(wtko$log2FoldChange) > fold_change,]$ensembl_gene_id
dko_counts <- read.table('~/Documents/Projects/Thymocyte_HiC/RNA_seq_counts/RAWCOUNTS_WT_and_DKO.tsv',header=T)
dko_counts$status <- 'unchanged'
dko_counts[dko_counts$id %in% dko_difexp ,]$status <- 'Diff'
dko_counts$wt <- (dko_counts$WT_CD69nDP_R1 + dko_counts$WT_CD69nDP_R2 + dko_counts$WT_CD69nDP_R3) /3
dko_counts$ko <- (dko_counts$DKO_CD69nDP_R1 + dko_counts$DKO_CD69nDP_R2 + dko_counts$DKO_CD69nDP_R3) /3
#
dev_difexp <- wtwt[wtwt$padj < 0.05 & abs(wtwt$log2FoldChange) > fold_change,]$ensembl_gene_id
dev_counts <- read.table('~/Documents/Projects/Thymocyte_HiC/RNA_seq_counts/RAWCOUNTS_WT_and_DKO.tsv',header = T)
dev_counts$status <- 'unchanged'
dev_counts[dev_counts$id %in% dev_difexp ,]$status <- 'Diff'
dev_counts$wt1 <- (dev_counts$WT_CD69nDP_R1 + dev_counts$WT_CD69nDP_R1 + dev_counts$WT_CD69nDP_R3) /3
dev_counts$wt4 <- (dev_counts$WT_CD69pCD4SP_R1 + dev_counts$WT_CD69pCD4SP_R2 ) /2


ctcf_cc_plot <- ggplot(ctcf_counts,aes(wt,ko,color=status)) + geom_point(size=0.1,show.legend = F) + scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + scale_color_manual(values = c('red','black')) + 
  ylab('ΔCTCF log(counts)') + xlab('CD69-DP log(counts)')

rad21_cc_plot <- ggplot(rad21_counts,aes(wt,ko,color=status)) + geom_point(size=0.1,show.legend = F) + scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + scale_color_manual(values = c('red','black')) + 
  ylab('ΔRAD21 log(counts)') + xlab('CD69-DP log(counts)')

dko_cc_plot <- ggplot(dko_counts,aes(wt,ko,color=status)) + geom_point(size=0.1,show.legend = F) + scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10')  + scale_color_manual(values = c('red','black'))  + 
  ylab('ΔRAD21/ΔCTCF log(counts)') + xlab('CD69-DP log(counts)')
  #ylab(expression('RAD21'^'-/-'*'/CTCF'^'-/- log(counts)')) + xlab('CD69negDP log(counts)')

dev_cc_plot <- ggplot(dev_counts,aes(wt1,wt4,color=status)) + geom_point(size=0.1,show.legend = F) + scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + scale_color_manual(values = c('red','black'))  + 
  ylab('CD69+SP log(counts)') + xlab('CD69-DP log(counts)')

f2.loglogcounts <- plot_grid(dev_cc_plot,dko_cc_plot,ctcf_cc_plot,rad21_cc_plot,nrow=1)

ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/loglog.counts.pdf',
        f2.loglogcounts,width=7,height=2,units = 'in',device = cairo_pdf)
#should have a legend for red=dereg, black=stable?

#another place where could put genes with rad21 or ctcf etc
####
#ecdf of expression levels

dev_difUp <- wtwt[wtwt$padj < 0.05 & abs(wtwt$log2FoldChange) > fold_change & wtwt$log2FoldChange > fold_change,]$ensembl_gene_id
dev_difDown <- wtwt[wtwt$padj < 0.05 & abs(wtwt$log2FoldChange) > fold_change & wtwt$log2FoldChange < -fold_change,]$ensembl_gene_id

dko_difUp <- wtko[wtko$padj < 0.05 & abs(wtko$log2FoldChange) > fold_change & wtko$log2FoldChange > fold_change,]$ensembl_gene_id
dko_difDown <- wtko[wtko$padj < 0.05 & abs(wtko$log2FoldChange) > fold_change & wtko$log2FoldChange < -fold_change,]$ensembl_gene_id

ctcf_difUp <- wtctcf[wtctcf$padj < 0.05 & abs(wtctcf$log2FoldChange) > fold_change & wtctcf$log2FoldChange > fold_change,]$ensembl_gene_id
ctcf_difDown <- wtctcf[wtctcf$padj < 0.05 & abs(wtctcf$log2FoldChange) > fold_change & wtctcf$log2FoldChange < -fold_change,]$ensembl_gene_id

rad21_difUp <- wtrad21[wtrad21$padj < 0.05 & abs(wtrad21$log2FoldChange) > fold_change & wtrad21$log2FoldChange > fold_change,]$ensembl_gene_id
rad21_difDown <- wtrad21[wtrad21$padj < 0.05 & abs(wtrad21$log2FoldChange) > fold_change & wtrad21$log2FoldChange < -fold_change,]$ensembl_gene_id

#or  
dev_counts[dev_counts$id %in% dev_difDown ,]$status <- 'Down'
dev_counts[dev_counts$id %in% dev_difUp ,]$status <- 'Up'
ggplot(dev_counts[dev_counts$wt1 >10,],aes(x=wt1,color=status)) + stat_ecdf() + scale_x_continuous(trans='log10') + theme_classic() + ylab('Proportion') + xlab('CD69negDP WT counts per gene (FPKM)')
ggplot(dev_counts[dev_counts$wt1 >10,],aes(x=wt4,color=status)) + stat_ecdf() + scale_x_continuous(trans='log10') + theme_classic() + ylab('Proportion') + xlab('CD69posSP WT counts per gene (FPKM)')

#for direct comparison
#dev_counts$status <- 'unchanged'
#dev_counts[dev_counts$id %in% dev_difDown & -(dev_counts$id %in% dko_difDown),]$status <- 'Down'
#dev_counts[dev_counts$id %in% dev_difUp & -(dev_counts$id %in% dko_difUp),]$status <- 'Up'
#
dko_counts[dko_counts$id %in% dko_difDown ,]$status <- 'Down'
dko_counts[dko_counts$id %in% dko_difUp ,]$status <- 'Up'
ggplot(dko_counts[dko_counts$wt >10,],aes(x=wt,color=status)) + stat_ecdf() + scale_x_continuous(trans='log10') + theme_classic() + ylab('Proportion') + xlab('CD69negDP WT counts per gene (FPKM)')
ggplot(dko_counts[dko_counts$wt >10,],aes(x=ko,color=status)) + stat_ecdf() + scale_x_continuous(trans='log10') + theme_classic() + ylab('Proportion') + xlab('CD69negDP CTCF-/RAD21- counts per gene (FPKM)')

#
ctcf_counts[ctcf_counts$id %in% ctcf_difDown ,]$status <- 'Down'
ctcf_counts[ctcf_counts$id %in% ctcf_difUp ,]$status <- 'Up'
ggplot(ctcf_counts[ctcf_counts$wt >10,],aes(x=wt,color=status)) + stat_ecdf() + scale_x_continuous(trans='log10') + theme_classic() + ylab('Proportion')
ggplot(ctcf_counts[ctcf_counts$wt >10,],aes(x=ko,color=status)) + stat_ecdf() + scale_x_continuous(trans='log10') + theme_classic() + ylab('Proportion')
#
rad21_counts[rad21_counts$id %in% rad21_difDown ,]$status <- 'Down'
rad21_counts[rad21_counts$id %in% rad21_difUp ,]$status <- 'Up'
ggplot(rad21_counts[rad21_counts$wt >10,],aes(x=wt,color=status)) + stat_ecdf() + scale_x_continuous(trans='log10') + theme_classic() + ylab('Proportion')
ggplot(rad21_counts[rad21_counts$wt >10,],aes(x=ko,color=status)) + stat_ecdf() + scale_x_continuous(trans='log10') + theme_classic() + ylab('Proportion')



# getting overlaps (for venn diagraming)

#get gene_ids common across all studies, as means the universe for the fishers or permutation test is consistent
gene_ids_all <- unique(c(as.character(wtwt$ensembl_gene_id),as.character(wtko$ensembl_gene_id),as.character(wtctcf$ensembl_gene_id),as.character(wtrad21$ensembl_gene_id)))
gene_ids_all <- gene_ids_all[gene_ids_all %in% wtwt$ensembl_gene_id & gene_ids_all %in% wtko$ensembl_gene_id & gene_ids_all %in% wtctcf$ensembl_gene_id & gene_ids_all %in% wtrad21$ensembl_gene_id]

get_contingency_table <- function(set1upID,set1downID,set2upID,set2downID, gene_ids_all = gene_ids_all){
  no_change_set1 <- gene_ids_all[!(gene_ids_all %in% set1upID | gene_ids_all %in% set1downID)]
  change_set1 <- gene_ids_all[gene_ids_all %in% set1upID | gene_ids_all %in% set1downID]
  no_change_set2 <- gene_ids_all[!(gene_ids_all %in% set2upID | gene_ids_all %in% set2downID)]
  change_set2 <- gene_ids_all[gene_ids_all %in% set2upID | gene_ids_all %in% set2downID]
  
  neither_change <- sum(no_change_set1 %in% no_change_set2)
  set2_only_change <- sum(change_set2 %in% no_change_set1)
  set1_only_change <- sum(change_set1 %in% no_change_set2)
  both_change <- sum(change_set1 %in% change_set2)
  return(c(neither_change,set2_only_change,set1_only_change,both_change))
  
}

#so above isn't looking at direction of change. I assume is because liz didn't
get_contingency_table.updown.concordance <- function(set1upID,set1downID,set2upID,set2downID){

  up.up <- length(set1upID[set1upID %in% set2upID])
  down.up <- length(set1downID[set1downID %in% set2upID])
  up.down <- length(set1upID[set1upID %in% set2downID])
  down.down <- length(set1downID[set1downID %in% set2downID])
  
  return(c(up.up,down.up,up.down,down.down))
  
}

#dev_ctcf.con <- get_contingency_table(dev_difUp,dev_difDown,ctcf_difUp,ctcf_difDown)

##getting fisher exact test values for overlaps across sets
#TeaTesting <-  matrix(c(dev_ctcf.con[1],dev_ctcf.con[2],dev_ctcf.con[3] , dev_ctcf.con[4]),
#        nrow=2, dimnames = list(Dev = c("Unchanged", "Changed"),
#                                 KO = c("Unchanged", "Changed")))

#fisher.test(TeaTesting,alternative = 'greater')

ttest_genes <- function(xUp,xDown, yUp, yDown, gene_ids=gene_ids_all){
  con.table <-  get_contingency_table(xUp,xDown, yUp, yDown, gene_ids)
  TeaTesting <-  matrix(c(con.table[1],con.table[2],con.table[3] , con.table[4]),
                      nrow=2, dimnames = list(A = c("Unchanged", "Changed"),
                                              B = c("Unchanged", "Changed")))
  return(fisher.test(TeaTesting,alternative = 'greater'))
  }


if(filtered == T){
  # fc.data.gr.filtered is genes that do *NOT* have proximal RAD21 or CTCF. SO the X are those that do not overlap <- check
  gene_ids_allX <- gene_ids_all[gene_ids_all %in% fc.data.gr.filtered$ensembl_gene_id]
  dev_difUpX <- dev_difUp[dev_difUp %in% fc.data.gr.filtered$ensembl_gene_id] 
  dev_difDownX <- dev_difDown[dev_difDown %in% fc.data.gr.filtered$ensembl_gene_id] 
  ctcf_difUpX <- ctcf_difUp[ctcf_difUp %in% fc.data.gr.filtered$ensembl_gene_id] 
  ctcf_difDownX <- ctcf_difDown[ctcf_difDown %in% fc.data.gr.filtered$ensembl_gene_id] 
  rad21_difUpX <- rad21_difUp[rad21_difUp %in% fc.data.gr.filtered$ensembl_gene_id]
  rad21_difDownX <- rad21_difDown[rad21_difDown %in% fc.data.gr.filtered$ensembl_gene_id] 
  dko_difUpX <- dko_difUp[dko_difUp %in% fc.data.gr.filtered$ensembl_gene_id] 
  dko_difDownX <- dko_difDown[dko_difDown %in% fc.data.gr.filtered$ensembl_gene_id]

  dev_difUp <- dev_difUp[!dev_difUp %in% fc.data.gr.filtered$ensembl_gene_id] 
  dev_difDown <- dev_difDown[!dev_difDown %in% fc.data.gr.filtered$ensembl_gene_id] 
  ctcf_difUp <- ctcf_difUp[!ctcf_difUp %in% fc.data.gr.filtered$ensembl_gene_id] 
  ctcf_difDown <- ctcf_difDown[!ctcf_difDown %in% fc.data.gr.filtered$ensembl_gene_id] 
  rad21_difUp <- rad21_difUp[!rad21_difUp %in% fc.data.gr.filtered$ensembl_gene_id]
  rad21_difDown <- rad21_difDown[!rad21_difDown %in% fc.data.gr.filtered$ensembl_gene_id] 
  dko_difUp <- dko_difUp[!dko_difUp %in% fc.data.gr.filtered$ensembl_gene_id] 
  dko_difDown <- dko_difDown[!dko_difDown %in% fc.data.gr.filtered$ensembl_gene_id]
  gene_ids_all <- gene_ids_all[!gene_ids_all %in% fc.data.gr.filtered$ensembl_gene_id]

  
}

dev.vs.ctcf.odds <- ttest_genes(dev_difUp,dev_difDown,ctcf_difUp,ctcf_difDown)
dev.vs.rad21.odds <- ttest_genes(dev_difUp,dev_difDown,rad21_difUp,rad21_difDown)
dev.vs.dko.odds <- ttest_genes(dev_difUp,dev_difDown,dko_difUp,dko_difDown)

rad21.vs.ctcf.odds <- ttest_genes(rad21_difUp,rad21_difDown,ctcf_difUp,ctcf_difDown)
rad21.vs.dko.odds <- ttest_genes(rad21_difUp,rad21_difDown,dko_difUp,dko_difDown)
ctcf.vs.dko.odds <- ttest_genes(ctcf_difUp,ctcf_difDown,dko_difUp,dko_difDown)



#all knockouts genes are enriched for dev. Higher if use higher thresholds, and with each other
######

#########Fishers test for convergence

ttest_genes.concordant <- function(xUp,xDown, yUp, yDown){
  con.table <-  get_contingency_table.updown.concordance(xUp,xDown, yUp, yDown)
  TeaTesting <-  matrix(c(con.table[1],con.table[2],con.table[3] , con.table[4]),
                        nrow=2, dimnames = list(A = c("Unchanged", "Changed"),
                                                B = c("Unchanged", "Changed")))
  return(fisher.test(TeaTesting,alternative = 'greater'))
}

dev.vs.ctcf.odds.c <- ttest_genes.concordant(dev_difUp,dev_difDown,ctcf_difUp,ctcf_difDown)
dev.vs.rad21.odds.c <- ttest_genes.concordant(dev_difUp,dev_difDown,rad21_difUp,rad21_difDown)
dev.vs.dko.odds.c <- ttest_genes.concordant(dev_difUp,dev_difDown,dko_difUp,dko_difDown)

rad21.vs.ctcf.odds.c <- ttest_genes.concordant(rad21_difUp,rad21_difDown,ctcf_difUp,ctcf_difDown)
rad21.vs.dko.odds.c <- ttest_genes.concordant(rad21_difUp,rad21_difDown,dko_difUp,dko_difDown)
ctcf.vs.dko.odds.c <- ttest_genes.concordant(ctcf_difUp,ctcf_difDown,dko_difUp,dko_difDown)

if(filtered == T){
  dev.vs.ctcf.oddsX <- ttest_genes(dev_difUpX,dev_difDownX,ctcf_difUpX,ctcf_difDownX,gene_ids = gene_ids_allX)
  dev.vs.rad21.oddsX <- ttest_genes(dev_difUpX,dev_difDownX,rad21_difUpX,rad21_difDownX,gene_ids = gene_ids_allX)
  dev.vs.dko.oddsX <- ttest_genes(dev_difUpX,dev_difDownX,dko_difUpX,dko_difDownX,gene_ids = gene_ids_allX)
  
  rad21.vs.ctcf.oddsX <- ttest_genes(rad21_difUpX,rad21_difDownX,ctcf_difUpX,ctcf_difDownX,gene_ids = gene_ids_allX)
  rad21.vs.dko.oddsX <- ttest_genes(rad21_difUpX,rad21_difDownX,dko_difUpX,dko_difDownX,gene_ids = gene_ids_allX)
  ctcf.vs.dko.oddsX <- ttest_genes(ctcf_difUpX,ctcf_difDownX,dko_difUpX,dko_difDownX,gene_ids = gene_ids_allX)
  
  dev.vs.ctcf.odds.cX <- ttest_genes.concordant(dev_difUpX,dev_difDownX,ctcf_difUpX,ctcf_difDownX)
  dev.vs.rad21.odds.cX <- ttest_genes.concordant(dev_difUpX,dev_difDownX,rad21_difUpX,rad21_difDownX)
  dev.vs.dko.odds.cX <- ttest_genes.concordant(dev_difUpX,dev_difDownX,dko_difUpX,dko_difDownX)
  
  rad21.vs.ctcf.odds.cX <- ttest_genes.concordant(rad21_difUpX,rad21_difDownX,ctcf_difUpX,ctcf_difDownX)
  rad21.vs.dko.odds.cX <- ttest_genes.concordant(rad21_difUpX,rad21_difDownX,dko_difUpX,dko_difDownX)
  ctcf.vs.dko.odds.cX <- ttest_genes.concordant(ctcf_difUpX,ctcf_difDownX,dko_difUpX,dko_difDownX)

}
#ok so when Matthias and Ediem say cohesin and ctcf ko are uncorrelated, they're wrong
#but direction of change isn't correlated

###### making a small table of odds ratios (concordance and not for has/doesn't have ctcf and rad21)

odds.table <- DataFrame(Condition = c('CTCFko','RAD21ko','DKO','CTCFko','RAD21ko','DKO'),
          'Proximal Marks' = c('RAD21 & CTCF','RAD21 & CTCF','RAD21 & CTCF','None','None','None'),
          Odds_olap = c(dev.vs.ctcf.odds$estimate , dev.vs.rad21.odds$estimate , dev.vs.dko.odds$estimate,
                   dev.vs.ctcf.oddsX$estimate , dev.vs.rad21.oddsX$estimate , dev.vs.dko.oddsX$estimate),
          Odds_concordance = c(dev.vs.ctcf.odds.c$estimate , dev.vs.rad21.odds.c$estimate , dev.vs.dko.odds.c$estimate,
                               dev.vs.ctcf.odds.cX$estimate , dev.vs.rad21.odds.cX$estimate , dev.vs.dko.odds.cX$estimate)
)

write.table(odds.table,'~/Documents/Projects/Thymocyte_HiC/concordance.table.tsv',sep='\t',row.names = F,quote = F)


#######

#library(VennDiagram)

dev.vdp <- c(as.character(dev_difDown),as.character(dev_difUp))
dko.vdp <- c(as.character(dko_difDown),as.character(dko_difUp))

ctcf.vdp <- c(as.character(ctcf_difUp), as.character(ctcf_difDown))
rad21.vdp <- c(as.character(rad21_difUp), as.character(rad21_difDown))

#overlap <- calculate.overlap(x = list("Dev" = dev.vdp,"dko" = dko.vdp))
venn.diagram(list("Dev Reg genes" = dev.vdp,"DKO sensitive genes" = dko.vdp),'~/Documents/woofers.tiff')
#venn.diagram(list("Developmental" = dev.vdp,"CTCF/RAD21" = dko.vdp,'CTCF KO' = ctcf.vdp,'RAD21 KO'=rad21.vdp),'~/Documents/woofers.tiff')
#
devNot.vdp <- as.character(wtwt[wtwt$padj > 0.05,]$ensembl_gene_id)
dkoNot.vdp <- as.character(wtko[wtko$padj > 0.05,]$ensembl_gene_id)
venn.diagram(list("Dev maintained genes" = devNot.vdp,"DKO sensitive genes" = dko.vdp),'~/Documents/woofers2.tiff')
venn.diagram(list("Dev Reg genes" = dev.vdp,"DKO insensitive genes" = dkoNot.vdp),'~/Documents/woofers3.tiff')



# this is cluttered. Trying upset plot. SHould also do this for concordant genes
library(UpSetR)
library(Cairo)
CairoPDF('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/upsetplots.pdf')


o <- list("Developmental" = dev.vdp,"ΔCTCF/ΔRAD21" = dko.vdp,'ΔCTCF' = ctcf.vdp,'ΔRAD21'=rad21.vdp)
upset(fromList(o), order.by = "freq",text.scale = 0.9)

##
o <- list("ΔCTCF/ΔRAD21 Up" = dko_difUp, "ΔCTCF/ΔRAD21 Down" = dko_difDown,'ΔCTCF Up' = ctcf_difUp,
          'ΔCTCF Down' = ctcf_difDown,'ΔRAD21 Up'=rad21_difUp,'ΔRAD21 Down'=rad21_difDown)
upset(fromList(o), order.by = "freq",text.scale = 0.9)  #this doesn't actually give all the overlaps, so don't use

o <- list('ΔCTCF Up' = ctcf_difUp,'ΔCTCF Down' = ctcf_difDown,'ΔRAD21 Up'=rad21_difUp,'ΔRAD21 Down'=rad21_difDown)
upset(fromList(o), order.by = "freq",text.scale = 0.9)  #this doesn't actually give all the overlaps, so don't use


dev.off()
#########
########barplots showing proportion of devreg/nondev reg genes down/upregulated by ctcfko, radko21 etc

#start with disreg ctcfko/ disreg rad21ko
dev.vdp <- c(as.character(dev_difDown),as.character(dev_difUp))
nondev.vdp <- gene_ids_all[!gene_ids_all %in% dev.vdp]


ctcf.d.in <- 100 * length(ctcf_difDown[ctcf_difDown %in% dev.vdp]) / length(dev.vdp)
ctcf.u.in <- 100 * length(ctcf_difUp[ctcf_difUp %in% dev.vdp])  / length(dev.vdp)
 
rad21.d.in <- 100 * length(rad21_difDown[rad21_difDown %in% dev.vdp])   / length(dev.vdp)
rad21.u.in <- 100 * length(rad21_difUp[rad21_difUp %in% dev.vdp])  / length(dev.vdp)

dko.d.in <- 100 * length(dko_difDown[dko_difDown %in% dev.vdp])  / length(dev.vdp)
dko.u.in <- 100 * length(dko_difUp[dko_difUp %in% dev.vdp])  / length(dev.vdp)
 
#
ctcf.d.out <- 100 * length(ctcf_difDown[ctcf_difDown %in% nondev.vdp]) / length(nondev.vdp)
ctcf.u.out <- 100 * length(ctcf_difUp[ctcf_difUp %in% nondev.vdp])  / length(nondev.vdp)
 
rad21.d.out <- 100 * length(rad21_difDown[rad21_difDown %in% nondev.vdp])   / length(nondev.vdp)
rad21.u.out <- 100 *length(rad21_difUp[rad21_difUp %in% nondev.vdp])  / length(nondev.vdp)
 
dko.d.out <- 100 * length(dko_difDown[dko_difDown %in% nondev.vdp])  / length(nondev.vdp)
dko.u.out <- 100 * length(dko_difUp[dko_difUp %in% nondev.vdp])  / length(nondev.vdp)

dev.olap.df <- data.frame('dev' = c(rep('Developmentally Regulated',6),rep('Stably Regulated',6)),
           'KO_direction' = rep(c('Down','Up'), 6), 
           'Percentage' = c(ctcf.d.in,ctcf.u.in,rad21.d.in,rad21.u.in,dko.d.in,dko.u.in, 
                            ctcf.d.out,ctcf.u.out,rad21.d.out,rad21.u.out,dko.d.out,dko.u.out),
           'KO' = c(rep(c('ΔCTCF','ΔCTCF','ΔRAD21','ΔRAD21','ΔCTCF/ΔRAD21','ΔCTCF/ΔRAD21'),2))
           )
 
library(ggsci)  

f4.withmarks.bars <- ggplot(dev.olap.df,aes(x=dev,y=Percentage,fill=KO_direction)) + geom_bar(stat='identity') + 
  coord_flip() + theme_classic() + scale_fill_jco() + facet_grid(.~KO) + ylim(0,26) + ylab('% of Genes Deregulated with Proximal RAD21&CTCF')

if(filtered == F){
  f4.withmarks.bars <- ggplot(dev.olap.df,aes(x=dev,y=Percentage,fill=KO_direction)) + geom_bar(stat='identity') + 
  coord_flip() + theme_classic() + scale_fill_jco() + facet_grid(.~KO) + ylim(0,26) + ylab('% of Genes Deregulated')

  ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/barcharts_koolapdev.pdf',
          f4.withmarks.bars,width=7,height=1,units = 'in',device = cairo_pdf)
}

#export 1200x300

#when filtered==T, note how CTCFko upregulated genes lose differences between stable and Dev Reg genes 
# in fact, dev.vs.ctcf.odds has a pval of 0.125, whilst dev.vs.ctcf.oddsX is < 2.2e-16
###
####

dev.vdpX <- c(as.character(dev_difDownX),as.character(dev_difUpX))
nondev.vdpX <- gene_ids_allX[!gene_ids_allX %in% dev.vdpX]


ctcf.d.in <- 100 * length(ctcf_difDownX[ctcf_difDownX %in% dev.vdpX]) / length(dev.vdpX)
ctcf.u.in <- 100 * length(ctcf_difUpX[ctcf_difUpX %in% dev.vdpX])  / length(dev.vdpX)

rad21.d.in <- 100 * length(rad21_difDownX[rad21_difDownX %in% dev.vdpX])   / length(dev.vdpX)
rad21.u.in <- 100 * length(rad21_difUpX[rad21_difUpX %in% dev.vdpX])  / length(dev.vdpX)

dko.d.in <- 100 * length(dko_difDownX[dko_difDownX %in% dev.vdpX])  / length(dev.vdpX)
dko.u.in <- 100 * length(dko_difUpX[dko_difUpX %in% dev.vdpX])  / length(dev.vdpX)

#
ctcf.d.out <- 100 * length(ctcf_difDownX[ctcf_difDownX %in% nondev.vdpX]) / length(nondev.vdpX)
ctcf.u.out <- 100 * length(ctcf_difUpX[ctcf_difUpX %in% nondev.vdpX])  / length(nondev.vdpX)

rad21.d.out <- 100 * length(rad21_difDownX[rad21_difDownX %in% nondev.vdpX])   / length(nondev.vdpX)
rad21.u.out <- 100 *length(rad21_difUpX[rad21_difUpX %in% nondev.vdpX])  / length(nondev.vdpX)

dko.d.out <- 100 * length(dko_difDownX[dko_difDownX %in% nondev.vdpX])  / length(nondev.vdpX)
dko.u.out <- 100 * length(dko_difUpX[dko_difUpX %in% nondev.vdpX])  / length(nondev.vdpX)

dev.olap.dfX <- data.frame('dev' = c(rep('Developmentally Regulated',6),rep('Stably Regulated',6)),
                          'KO_direction' = rep(c('Down','Up'), 6), 
                          'Percentage' = c(ctcf.d.in,ctcf.u.in,rad21.d.in,rad21.u.in,dko.d.in,dko.u.in, 
                                           ctcf.d.out,ctcf.u.out,rad21.d.out,rad21.u.out,dko.d.out,dko.u.out),
                          'KO' = c(rep(c('ΔCTCF','ΔCTCF','ΔRAD21','ΔRAD21','ΔCTCF/ΔRAD21','ΔCTCF/ΔRAD21'),2))
)


f4.nomarks.bars <- ggplot(dev.olap.dfX,aes(x=dev,y=Percentage,fill=KO_direction)) + geom_bar(stat='identity') + 
  coord_flip() + theme_classic() + scale_fill_jco() + facet_grid(.~KO) + ylim(0,26) + ylab('% of Genes Deregulated no Proximal RAD21&CTCF')

kos_with_bars <- plot_grid(f4.withmarks.bars,f4.nomarks.bars,ncol = 1) #these are all different, despite how similar the 
#note that I'm splitting out genes that have *either* CTCF or RAD21. Looking separately is also interesting
#and gives some odd looking results
if(filtered == T){
  ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/barcharts_koolapdev.marks.pdf',
        kos_with_bars,width=7,height=2,units = 'in',device = cairo_pdf)
  }

