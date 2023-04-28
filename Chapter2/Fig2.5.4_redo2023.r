source('~/Documents/Projects/Thymocyte_HiC/Thesis/scripts/general_functions.R')

library(ggplot2)
library(RColorBrewer)
library(regioneR)
library(cowplot)
library(ggsci)
library(ggpubr)

setwd('~/Documents/Projects/Thymocyte_HiC/DEseq2_in_vivo_Development_YaRNAseq/')

theme_set(theme_classic(base_size = 6))

fold_change <- 0
max_dist <- 500
min_dist <- 1

wtko <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_DKO_CD69nDPWT1-3_vs_CD69nDPDKO1-3.csv',header=T)
wtko <- wtko[!is.na(wtko$padj),]

wtwt <- read.csv('04_CD69nDPWT1-3_vs_CD69pCD4SPWT1-2.csv',header=T)
wtwt <- wtwt[!is.na(wtwt$padj),]  #genes that are na here are na as their basemean is super low


wtctcf  <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_CTCFKO_CD69nDPWT4-5_vs_CD69nDPCTCFKO1-2.csv',header=T)
wtctcf <- wtctcf[!is.na(wtctcf$padj),]


wtrad21 <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_Rad21KO_DPWT1-2_vs_DPRad21KO1-2.csv',header=T)
wtrad21 <- wtrad21[!is.na(wtrad21$padj),]

#remove duplicates
wtko <- wtko[order(wtko$ensembl_gene_id, abs(wtko$pvalue) ), ]
wtko <-  wtko[ !duplicated(wtko$ensembl_gene_id), ] 

wtwt <- wtwt[order(wtwt$ensembl_gene_id, abs(wtwt$pvalue) ), ]
wtwt <-  wtwt[ !duplicated(wtwt$ensembl_gene_id), ] 

wtctcf <- wtctcf[order(wtctcf$ensembl_gene_id, abs(wtctcf$pvalue) ), ]
wtctcf <-  wtctcf[ !duplicated(wtctcf$ensembl_gene_id), ] 

wtrad21 <- wtrad21[order(wtrad21$ensembl_gene_id, abs(wtrad21$pvalue) ), ]
wtrad21 <-  wtrad21[ !duplicated(wtrad21$ensembl_gene_id), ] 

#########################################################################################################
#########################################################################################################


merge_diffgene_exps <- function(diffgene1,diffgene2,pvalcut1,pvalcut2,logfoldcut1,logfoldcut2,cond1,cond2){
  merged_genes <- merge(diffgene1,diffgene2,by='ensembl_gene_id')
  #fold_changes[is.na(fold_changes)] <- 1 #remove the padj which go to NA, upset stuff
  merged_genes$status <- 'unchanged'
  merged_genes[merged_genes$padj.x <= pvalcut1 & abs(merged_genes$log2FoldChange.x) >logfoldcut1,]$status <- cond1
  merged_genes[merged_genes$padj.y <= pvalcut2 & abs(merged_genes$log2FoldChange.y) >logfoldcut2,]$status <- cond2
  merged_genes[merged_genes$padj.y <= pvalcut2 & merged_genes$padj.x <= pvalcut1 & abs(merged_genes$log2FoldChange.x) >logfoldcut1 & abs(merged_genes$log2FoldChange.y) >logfoldcut2,]$status <- 'Both'
  return(merged_genes)
}

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
  
  changes.gr$log2FoldChange.x <- changes$log2FoldChange.x
  changes.gr$padj.x  <- changes$padj.x
  changes.gr$log2FoldChange.y <- changes$log2FoldChange.y
  changes.gr$padj.y  <- changes$padj.y
  return(changes.gr)
}

dev_vs_dko <- difexp_to_granges(merge_diffgene_exps(wtwt,wtko,0.05,0.05,1.5,1.5,'Development','ΔCTCF/ΔRAD21'))
dev_vs_ctcf <- difexp_to_granges(merge_diffgene_exps(wtwt,wtctcf,0.05,0.05,1.5,1.5,'Development','ΔCTCF'))
dev_vs_rad21 <- difexp_to_granges(merge_diffgene_exps(wtwt,wtrad21,0.05,0.05,1.5,1.1,'Development','ΔRAD21'))
ctcf_vs_dko <- difexp_to_granges(merge_diffgene_exps(wtctcf,wtko,0.05,0.05,1.5,1.5,'ΔCTCF','ΔCTCF/ΔRAD21'))
ctcf_vs_rad21 <- difexp_to_granges(merge_diffgene_exps(wtctcf,wtrad21,0.05,0.05,1.5,1.1,'ΔCTCF','ΔRAD21'))
rad21_vs_dko <- difexp_to_granges(merge_diffgene_exps(wtrad21,wtko,0.05,0.05,1.1,1.5,'ΔRAD21','ΔCTCF/ΔRAD21'))



#####################################################################################
#####################################################################################
#all the below until stated should be put int oa functin so can switch out which comp fc.data.gr refers to


#fc.data <- dev_vs_dko
get_props <- function(fc.data.gr,max_dist, KO_name,min_dist = 0){

  fc.data.gr <- promoters(fc.data.gr)
  ctcf.sites.gr <- rtracklayer::import.bed('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/1_CTCFCD69negDPWTR1/1_CTCFCD69negDPWTR1_peaks.subpeaks.bed')
  rad21.sites.gr <- rtracklayer::import.bed('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/Rad21_ChIPseq_in_CD69negDP/WT_CD69negDP/Rad21CD69negDPWTR1_mapped_sorted_RemoveDuplicates_peaks.subpeaks.bed')
  
  #########################################################################################################################
  #########################################################################################################################
  #KO_name='ΔCTCF/ΔRAD21'  ΔCTCF ΔRAD21
  
  #maybe make the distances non-overlapping
  
  
  fc.data.gr$mark <- 'None'
  fc.data.gr[queryHits(findOverlaps(fc.data.gr,ctcf.sites.gr,maxgap = max_dist))]$mark <- 'CTCF'
  #fc.data.gr[queryHits(findOverlaps(fc.data.gr,rad21.sites.gr,maxgap = max_dist))]$mark <- 'RAD21'
  hasCTCF.gr <- fc.data.gr[fc.data.gr$mark == 'CTCF'] 
  hasCTCF.gr[queryHits(findOverlaps(hasCTCF.gr,rad21.sites.gr,maxgap = max_dist))]$mark <- 'Both'
  noCTCF.gr <- fc.data.gr[fc.data.gr$mark != 'CTCF'] 
  noCTCF.gr[queryHits(findOverlaps(noCTCF.gr,rad21.sites.gr,maxgap = max_dist))]$mark <- 'RAD21'
  
  if(min_dist > 0){
    fc.min <- fc.data.gr
    fc.min$mark <- 'None'
    fc.min[queryHits(findOverlaps(fc.min,ctcf.sites.gr,maxgap = min_dist))]$mark <- 'CTCF'
    #fc.min[queryHits(findOverlaps(fc.min,rad21.sites.gr,maxgap = min_dist))]$mark <- 'RAD21'
    hasMinCTCF.gr <- fc.min[fc.min$mark == 'CTCF'] 
    hasMinCTCF.gr[queryHits(findOverlaps(hasMinCTCF.gr,rad21.sites.gr,maxgap = max_dist))]$mark <- 'Both'
    hasCTCF.gr <- hasCTCF.gr[-unique(queryHits(findOverlaps(hasCTCF.gr,hasMinCTCF.gr)))]
  }
  
  fc.data.gr <- c(noCTCF.gr,hasCTCF.gr)
  table(fc.data.gr$mark)
  
  genes_no_binding.gr <- fc.data.gr[fc.data.gr$mark == 'None']
  genes_with_binding.gr <- fc.data.gr[fc.data.gr$mark == 'Both'] # can set to both if prefer
  
  ###### Get developmentally regulated / stable genes with or without proximal binding
  
  #dev_down_nobind <- genes_no_binding.gr[genes_no_binding.gr$padj.x < 0.05 & genes_no_binding.gr$log2FoldChange.x < 0]
  #dev_up_nobind <- genes_no_binding.gr[genes_no_binding.gr$padj.x < 0.05 & genes_no_binding.gr$log2FoldChange.x > 0]
  dev_stab_nobind <- genes_no_binding.gr[genes_no_binding.gr$padj.x > 0.05 ]
  dev_reg_nobind <- genes_no_binding.gr[genes_no_binding.gr$padj.x < 0.05 ]
  
  #dev_down_Ybind <- genes_with_binding.gr[genes_with_binding.gr$padj.x < 0.05 & genes_with_binding.gr$log2FoldChange.x < 0]
  #dev_up_Ybind <- genes_with_binding.gr[genes_with_binding.gr$padj.x < 0.05 & genes_with_binding.gr$log2FoldChange.x > 0]
  dev_stab_Ybind <- genes_with_binding.gr[genes_with_binding.gr$padj.x > 0.05 ]
  dev_reg_Ybind <- genes_with_binding.gr[genes_with_binding.gr$padj.x < 0.05 ]
  
  
  length(dev_down_nobind) ; length(dev_up_nobind) ; length(dev_stab_nobind)
  length(dev_down_Ybind) ; length(dev_up_Ybind) ; length(dev_stab_Ybind)
  
  #######
  
  split_by_kofate <- function(genes){
    down <- genes[genes$padj.y < 0.05 & genes$log2FoldChange.y < 0]
    up <- genes[genes$padj.y < 0.05 & genes$log2FoldChange.y > 0]
    stab <- genes[genes$padj.y > 0.05 ]  
    return(list(KO.downprop=100*length(down)/length(genes),KO.upprop=100*length(up)/length(genes), nChange = length(up) + length(down), nNoChange = length(genes) - length(up) - length(down)    ))
    #return(list(KO.down = down,KO.up=up,KO.stab=stab,KO.downprop=100*length(down)/length(genes),KO.upprop=100*length(up)/length(genes)))
  }
  
  ###
  
  
  dev_stab_nobindProps <- split_by_kofate(dev_stab_nobind)
  dev_reg_nobindProps <- split_by_kofate(dev_reg_nobind)
  
  print(fisher.test(matrix(c(dev_stab_nobindProps$nChange,dev_stab_nobindProps$nNoChange,
                  dev_reg_nobindProps$nChange, dev_reg_nobindProps$nNoChange),ncol = 2)))
  
  dev_stab_YbindProps <- split_by_kofate(dev_stab_Ybind)
  dev_reg_YbindProps <- split_by_kofate(dev_reg_Ybind)
  
  print(fisher.test(matrix(c(dev_stab_YbindProps$nChange,dev_stab_YbindProps$nNoChange,
                  dev_reg_YbindProps$nChange, dev_reg_YbindProps$nNoChange),ncol = 2)))
  
  #
  
  #
  
  x<- data.frame(dev=c('Developmentally Regulated','Developmentally Regulated','Stably Regulated','Stably Regulated'),
             KO_direction=rep(c('Down','Up'),2),
             Percentage= c(dev_reg_nobindProps[[1]],dev_reg_nobindProps[[2]],
                           dev_stab_nobindProps[[1]],dev_stab_nobindProps[[2]]),
             KO=KO_name
             )
  
  
  y<- data.frame(dev=c('Developmentally Regulated','Developmentally Regulated','Stably Regulated','Stably Regulated'),
                 KO_direction=rep(c('Down','Up'),2),
                 Percentage= c(dev_reg_YbindProps[[1]],dev_reg_YbindProps[[2]],
                               dev_stab_YbindProps[[1]],dev_stab_YbindProps[[2]]),
                 KO=KO_name
  )
  return(list(x,y))
}

dev_ctcfko_dfs <- get_props(dev_vs_ctcf,500,"ΔCTCF")
dev_dko_dfs <- get_props(dev_vs_dko,500,"ΔCTCF/ΔRAD21")
dev_rad21ko_dfs <- get_props(dev_vs_rad21,500,"ΔRAD21")

noBinding <- rbind(dev_ctcfko_dfs[[1]],dev_dko_dfs[[1]],dev_rad21ko_dfs[[1]])
wBinding <- rbind(dev_ctcfko_dfs[[2]],dev_dko_dfs[[2]],dev_rad21ko_dfs[[2]])

noBinding$mark <- 'No proximal binding'
wBinding$mark <- 'CTCF & RAD21'
#
bp500 <- ggplot(rbind(noBinding,wBinding),aes(x=dev,y=Percentage,fill=KO_direction)) + geom_bar(stat='identity') + 
  coord_flip() + theme_classic() + scale_fill_jco() + facet_grid(mark~KO) + ylab('% of Genes Deregulated no Proximal RAD21&CTCF')

#ggplot(noBinding,aes(x=dev,y=Percentage,fill=KO_direction)) + geom_bar(stat='identity') + 
#  coord_flip() + theme_classic() + scale_fill_jco() + facet_grid(.~KO) + ylab('% of Genes Deregulated no Proximal RAD21&CTCF')

#ggplot(wBinding,aes(x=dev,y=Percentage,fill=KO_direction)) + geom_bar(stat='identity') + 
#  coord_flip() + theme_classic() + scale_fill_jco() + facet_grid(.~KO) + ylab('% of Genes Deregulated no Proximal RAD21&CTCF')

ggsave2('~/Desktop/Fig2.5.4.pdf',
        bp500,width=7,height=3.5,units = 'in',device = cairo_pdf)

