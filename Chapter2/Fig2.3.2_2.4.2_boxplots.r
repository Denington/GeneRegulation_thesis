#matched boxplots
#make sure to note that this is looking at concordance of genes deregulated in both conditions

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(cowplot)

theme_set(theme_classic(base_size = 8))


min_fold <- 1.2  #liz' value of 2  . 4000 diff genes through a couple of stages of development seems so high, and Rao and Liz use 2 fold cutoff on top of pval. I like having a cutoff of 2 for dev at least. 1/3 of genes is too many :/
pval <- 0.05
fpkm_thresh <- 1
ymax <- 2500

fpkms <- read.table('~/Documents/Projects/Thymocyte_HiC/Ya_RNA_seq/DPdev_WTvsDKO_average_FPKM_normalized.tsv',header = T,sep='\t')
fpkms <- fpkms[,c('Ensembl_gene_id','WT_CD69.DP','WT_CD69.CD4SP','DKO_CD69.DP','DKO_CD69.CD4SP')] #I think, names are less clear
fpkm_rad21 <- read.table('~/Documents/Projects/Thymocyte_HiC/FPKM_Rad21KO.tsv',sep='\t',header = T)
fpkm_ctcf <- read.table('~/Documents/Projects/Thymocyte_HiC/FPKM_CTCFKO.tsv',sep='\t',header = T)


setwd('~/Documents/Projects/Thymocyte_HiC/DEseq2_in_vivo_Development_YaRNAseq/')

wtko <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_DKO_CD69nDPWT1-3_vs_CD69nDPDKO1-3.csv',header=T)
wtko <- wtko[!is.na(wtko$padj),] #any replicates? weird how the numbers of up/down genes changes so much
wtwt <- read.csv('04_CD69nDPWT1-3_vs_CD69pCD4SPWT1-2.csv',header=T)
wtwt <- wtwt[!is.na(wtwt$padj),]  #genes that are na here are na as their basemean is super low

wtrad21 <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_Rad21KO_DPWT1-2_vs_DPRad21KO1-2.csv',header = T)
wtrad21 <- wtrad21[!is.na(wtrad21$padj),]  #genes that are na here are na as their basemean is super low
wtctcf <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_CTCFKO_CD69nDPWT4-5_vs_CD69nDPCTCFKO1-2.csv',header = T)
wtctcf <- wtctcf[!is.na(wtctcf$padj),]  #genes that are na here are na as their basemean is super low

### fig3
#### dev vs dko

fold_changes_dev <- merge(wtko,wtwt,by='ensembl_gene_id')
fold_changes_dev <- merge(fold_changes_dev,fpkms,by.x='ensembl_gene_id',by.y='Ensembl_gene_id')
xy <- fold_changes_dev[fold_changes_dev$padj.x < pval & fold_changes_dev$padj.y < pval & abs(fold_changes_dev$log2FoldChange.x) > min_fold & abs(fold_changes_dev$log2FoldChange.y) > min_fold,]
xy$direction <- 'ΔCTCF/ΔRAD21 Down'
xy[xy$log2FoldChange.x > 0,]$direction <- 'ΔCTCF/ΔRAD21 Up'
#ggpaired(xx,cond1= 'WT_CD69.DP',cond2 = 'DKO_CD69.DP',facet.by = 'direction',title = 'Dev vs DKO',point.size = 0,line.size = 0.25) + theme_bw() + scale_y_continuous(trans='log10') + ylab('Gene Expression (average FPKM)')   #limits = c(1,10000) in scale_y

#xx <- xx[ (xx$WT_CD69.DP + xx$DKO_CD69.DP )/2 > fpkm_thresh,] #fixed!
xy <- xy[xy$WT_CD69.DP > fpkm_thresh & xy$WT_CD69.CD4SP > fpkm_thresh,] #fixed! kind of

yx <- xy[,c('direction','WT_CD69.DP','WT_CD69.CD4SP')]
colnames(yx) <- c('direction','CD69DP','CD69SP')
yx$woofers <- seq(1:length(yx$direction))
yx <- gather(yx,key='Condition',value = 'FPKM',CD69DP,CD69SP) 
yx$Condition[yx$Condition == 'CD69SP'] <- 'SP'
yx$Condition[yx$Condition == 'CD69DP'] <- 'DP'

yx$Condition <- factor(yx$Condition,levels = c('DP','SP'))
f3.dko_vs_dev <- ggplot(yx,aes(x=Condition,y=FPKM,fill=Condition)) + geom_line(aes(group = woofers),alpha=0.4) + geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.7) +
  facet_wrap(~direction) + scale_y_continuous(trans = 'log10')  + coord_cartesian(ylim=c(1,ymax)) + 
  ggpubr::stat_compare_means(paired = T,label = 'p.format') +scale_fill_manual(values=c('#868686FF','#4A6990FF')) +
  theme(strip.background = element_rect(fill='#CD534CFF')) + theme(legend.position = 'None') # + theme(strip.text = element_text(colour='white'))

#0073C2FF

####


# rad21ko vs dev


fold_changes_rad21 <- merge(wtrad21,wtwt,by='ensembl_gene_id')
#so foldchange.x is ko, foldchange.y is developement
fold_changes_rad21 <- merge(fold_changes_rad21,fpkms,by.x='ensembl_gene_id',by.y='Ensembl_gene_id')
xx <- fold_changes_rad21[fold_changes_rad21$padj.x < pval & fold_changes_rad21$padj.y < pval & abs(fold_changes_rad21$log2FoldChange.x) > min_fold & abs(fold_changes_rad21$log2FoldChange.y) > min_fold,]
xx$direction <- 'woof'
xx[xx$log2FoldChange.x > 0,]$direction <- 'ΔRAD21 Up'
xx[xx$log2FoldChange.x < 0,]$direction <- 'ΔRAD21 Down'


xx <- xx[xx$WT_CD69.DP > fpkm_thresh & xx$WT_CD69.CD4SP > fpkm_thresh,] #fixed! kind of

yy <- xx[,c('direction','WT_CD69.DP','WT_CD69.CD4SP')]
colnames(yy) <- c('direction','WT','WT4')
yy$woofers <- seq(1:length(yy$direction))
yy <- gather(yy,key='Condition',value = 'FPKM',WT,WT4) 
yy$Condition[yy$Condition == 'WT4'] <- 'SP'
yy$Condition[yy$Condition == 'WT'] <- 'DP'


yy$Condition <- factor(yy$Condition,levels = c('DP','SP'))

f3.rad21ko_vs_dev <-  ggplot(yy,aes(x=Condition,y=FPKM,fill=Condition)) + geom_line(aes(group = woofers),alpha=0.4) + geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.7) + 
  facet_wrap(~direction) + scale_y_continuous(trans = 'log10') +coord_cartesian(ylim=c(1,ymax)) +
  ggpubr::stat_compare_means(paired = T,label = 'p.format') + scale_fill_manual(values=c('#868686FF','#4A6990FF')) +
  theme(strip.background = element_rect(fill='#EFC000E5')) + theme(legend.position = 'None') #+ theme(strip.text = element_text(colour='white'))

#Δ
####

fold_changes_ctcf <- merge(wtctcf,wtwt,by='ensembl_gene_id')
#so foldchange.x is ko, foldchange.y is developement
fold_changes_ctcf <- merge(fold_changes_ctcf,fpkms,by.x='ensembl_gene_id',by.y='Ensembl_gene_id')
#quite a few lost here, is weird, about half of all genes????
#fixed! Mad reason- Ya had mgi symbol as the first column, whenever there was a none r was ignoring the row!
xx <- fold_changes_ctcf[fold_changes_ctcf$padj.x < pval & fold_changes_ctcf$padj.y < pval & abs(fold_changes_ctcf$log2FoldChange.x) > min_fold & abs(fold_changes_ctcf$log2FoldChange.y) > min_fold,]
#probably do want the cut off here
xx$direction <- 'woof'
xx[xx$log2FoldChange.x > 0,]$direction <- 'ΔCTCF Up'
xx[xx$log2FoldChange.x < 0,]$direction <- 'ΔCTCF Down'


xx <- xx[xx$WT_CD69.DP > fpkm_thresh & xx$WT_CD69.CD4SP > fpkm_thresh,] #fixed! kind of

yy <- xx[,c('direction','WT_CD69.DP','WT_CD69.CD4SP')]
colnames(yy) <- c('direction','WT','WT4')
yy$woofers <- seq(1:length(yy$direction))
yy <- gather(yy,key='Condition',value = 'FPKM',WT,WT4) 
yy$Condition[yy$Condition == 'WT4'] <- 'SP'
yy$Condition[yy$Condition == 'WT'] <- 'DP'


yy$Condition <- factor(yy$Condition,levels = c('DP','SP'))


f3.ctcfko_vs_dev <- ggplot(yy,aes(x=Condition,y=FPKM,fill=Condition)) + geom_line(aes(group = woofers),alpha=0.4) + 
  geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.7) + facet_wrap(~direction) + scale_y_continuous(trans = 'log10') + 
  coord_cartesian(ylim=c(1,ymax)) +  ggpubr::stat_compare_means(paired = T,label = 'p.format') + 
  scale_fill_manual(values=c('#868686FF','#4A6990FF')) +
  theme(strip.background = element_rect(fill='#0073C2CC')) + theme(legend.position = 'None') #+ theme(strip.text = element_text(colour='white'))



######
######
#fig2
######
#######
min_fold <- 0.5

#Δ
###
#### dko vs rad21ko

fold_changes_dev <- merge(wtko,wtrad21,by='ensembl_gene_id')
#so foldchange.x is ko, foldchange.y is developement
fold_changes_dev <- merge(fold_changes_dev,fpkms,by.x='ensembl_gene_id',by.y='Ensembl_gene_id')
#quite a few lost here, is weird, about half of all genes????
xx <- fold_changes_dev[fold_changes_dev$padj.x < pval & fold_changes_dev$padj.y < pval & abs(fold_changes_dev$log2FoldChange.x) > min_fold & abs(fold_changes_dev$log2FoldChange.y) > min_fold,]
#probably do want the cut off here
xx$direction <- 'ΔRAD21 Down'
xx[xx$log2FoldChange.y > 0,]$direction <- 'ΔRAD21 Up'
#ggpaired(xx,cond1= 'WT_CD69.DP',cond2 = 'DKO_CD69.DP',facet.by = 'direction',title = 'Dev vs DKO',point.size = 0,line.size = 0.25) + theme_bw() + scale_y_continuous(trans='log10') + ylab('Gene Expression (average FPKM)')   #limits = c(1,10000) in scale_y

#xx <- xx[ (xx$WT_CD69.DP + xx$DKO_CD69.DP )/2 > fpkm_thresh,] #fixed!
xx <- xx[xx$WT_CD69.DP > fpkm_thresh & xx$DKO_CD69.DP > fpkm_thresh,] #fixed! kind of

yy <- xx[,c('direction','WT_CD69.DP','DKO_CD69.DP')]
colnames(yy) <- c('direction','CD69','ΔCTCFΔRAD21')
yy$woofers <- seq(1:length(yy$direction))
yy <- gather(yy,key='Condition',value = 'FPKM',CD69,ΔCTCFΔRAD21) 
yy$Condition[yy$Condition == 'ΔCTCFΔRAD21'] <- 'ΔCTCF/ΔRAD21'
yy$Condition[yy$Condition == 'CD69'] <- 'DP'

yy$Condition <- factor(yy$Condition,levels = c('DP','ΔCTCF/ΔRAD21'))
f2.dko_vs_Rad21ko <- ggplot(yy,aes(x=Condition,y=FPKM,fill=Condition)) + geom_line(aes(group = woofers),alpha=0.4) +
  geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.7) +facet_wrap(~direction) + scale_y_continuous(trans = 'log10')  + 
  coord_cartesian(ylim=c(1,ymax)) + ggpubr::stat_compare_means(paired = T,label = 'p.format') + scale_fill_manual(values=c('#868686FF','#CD534CFF')) +
  theme(strip.background = element_rect(fill='#EFC000E5')) + theme(legend.position = 'None') #+ theme(strip.text = element_text(colour='white'))




###
#### ctcfko vs  dko

fold_changes_dev <- merge(wtko,wtctcf,by='ensembl_gene_id')
#so foldchange.x is ko, foldchange.y is developement
fold_changes_dev <- merge(fold_changes_dev,fpkms,by.x='ensembl_gene_id',by.y='Ensembl_gene_id')
#quite a few lost here, is weird, about half of all genes????
xx <- fold_changes_dev[fold_changes_dev$padj.x < pval & fold_changes_dev$padj.y < pval & abs(fold_changes_dev$log2FoldChange.x) > min_fold & abs(fold_changes_dev$log2FoldChange.y) > min_fold,]
#probably do want the cut off here
xx$direction <- 'ΔCTCF Down'
xx[xx$log2FoldChange.y > 0,]$direction <- 'ΔCTCF Up'
#ggpaired(xx,cond1= 'WT_CD69.DP',cond2 = 'DKO_CD69.DP',facet.by = 'direction',title = 'Dev vs DKO',point.size = 0,line.size = 0.25) + theme_bw() + scale_y_continuous(trans='log10') + ylab('Gene Expression (average FPKM)')   #limits = c(1,10000) in scale_y

#xx <- xx[ (xx$WT_CD69.DP + xx$DKO_CD69.DP )/2 > fpkm_thresh,] #fixed!
xx <- xx[xx$WT_CD69.DP > fpkm_thresh & xx$DKO_CD69.DP > fpkm_thresh,] #fixed! kind of

yy <- xx[,c('direction','WT_CD69.DP','DKO_CD69.DP')]
colnames(yy) <- c('direction','CD69','ΔCTCFΔRAD21')
yy$woofers <- seq(1:length(yy$direction))
yy <- gather(yy,key='Condition',value = 'FPKM',CD69,ΔCTCFΔRAD21) 
yy$Condition[yy$Condition == 'ΔCTCFΔRAD21'] <- 'ΔCTCF/ΔRAD21'
yy$Condition[yy$Condition == 'CD69'] <- 'DP'

yy$Condition <- factor(yy$Condition,levels = c('DP','ΔCTCF/ΔRAD21'))
f2.dko_vs_ctcfko <- ggplot(yy,aes(x=Condition,y=FPKM,fill=Condition)) + geom_line(aes(group = woofers),alpha=0.4) + 
  geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.7) +facet_wrap(~direction) + scale_y_continuous(trans = 'log10')+
  coord_cartesian(ylim=c(1,ymax)) + ggpubr::stat_compare_means(paired = T,label = 'p.format')+ scale_fill_manual(values=c('#868686FF','#CD534CFF')) +
  theme(strip.background = element_rect(fill='#0073C2CC')) + theme(legend.position = 'None') #+ theme(strip.text = element_text(colour='white'))


#####

#### rad21ko vs ctcfko
#

fold_changes_rad21 <- merge(wtrad21,wtctcf,by='ensembl_gene_id')
#so foldchange.x is ko, foldchange.y is developement
fold_changes_rad21 <- merge(fold_changes_rad21,fpkm_rad21,by.x='ensembl_gene_id',by.y='Ensembl_gene_id')
#quite a few lost here, is weird, about half of all genes????
#fixed! Mad reason- Ya had mgi symbol as the first column, whenever there was a none r was ignoring the row!
xx <- fold_changes_rad21[fold_changes_rad21$padj.x < pval & fold_changes_rad21$padj.y < pval & abs(fold_changes_rad21$log2FoldChange.x) > min_fold & abs(fold_changes_rad21$log2FoldChange.y) > min_fold,]
#probably do want the cut off here
xx$direction <- 'woof'
xx[xx$log2FoldChange.y > 0,]$direction <- 'ΔCTCF Up'
xx[xx$log2FoldChange.y < 0,]$direction <- 'ΔCTCF Down'

xx$WT_fpkm <- (xx$WT1 + xx$WT2)/2
xx$RAD21_fpkm <- (xx$Rad21KO1 + xx$Rad21KO2)/2

#xx <- xx[(xx$WT_fpkm + xx$RAD21_fpkm)/2 > fpkm_thresh,] #fixed! kind of
#xx <- xx[xx$WT_fpkm > fpkm_thresh,] #fixed! kind of
xx <- xx[xx$WT_fpkm > fpkm_thresh & xx$RAD21_fpkm > fpkm_thresh,] #fixed! kind of

#ggpaired(xx,cond1= 'WT_fpkm',cond2 = 'RAD21_fpkm',facet.by = 'direction',title = 'Dev vs RAD21ko',point.size = 0,line.size = 0.25) + theme_bw() + scale_y_continuous(trans='log10') + ylab('Gene Expression (average FPKM)') 

#

yy <- xx[,c('direction','WT_fpkm','RAD21_fpkm')]
colnames(yy) <- c('direction','WT','RAD21ko')
yy$woofers <- seq(1:length(yy$direction))
yy <- gather(yy,key='Condition',value = 'FPKM',WT,RAD21ko) 
yy$Condition[yy$Condition == 'RAD21ko'] <- 'ΔRAD21'
yy$Condition[yy$Condition == 'WT'] <- 'DP'


yy$Condition <- factor(yy$Condition,levels = c('DP','ΔRAD21'))


f2.rad21ko_vs_ctcfko <-  ggplot(yy,aes(x=Condition,y=FPKM,fill=Condition)) + geom_line(aes(group = woofers),alpha=0.4) + 
  geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.7) +  facet_wrap(~direction) + scale_y_continuous(trans = 'log10') + 
  coord_cartesian(ylim=c(1,ymax)) + ggpubr::stat_compare_means(paired = T,label = 'p.format') + scale_fill_manual(values=c('#868686FF','#EFC000FF')) +
  theme(strip.background = element_rect(fill='#0073C2CC')) + theme(legend.position = 'None') #+ theme(strip.text = element_text(colour='white'))


kovsdev <- plot_grid(f3.ctcfko_vs_dev, f3.rad21ko_vs_dev,f3.dko_vs_dev, labels = c('A', 'B','C'), label_size = 12,nrow=1)


kovsko <- plot_grid(f2.rad21ko_vs_ctcfko, f2.dko_vs_ctcfko,f2.dko_vs_Rad21ko, labels = c('A', 'B','C'), label_size = 12,nrow=1)

ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/kovsdev.matchbplots.pdf',
        kovsdev,width=7,height=3,units = 'in',device = cairo_pdf)

ggsave2('~/Documents/Projects/Thymocyte_HiC/Thesis/presentation_figs/kovsko.matchbplots.pdf',
        kovsko,width=7,height=3,units = 'in',device = cairo_pdf)


