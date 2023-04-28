source('~/Documents/Projects/Thymocyte_HiC/Thesis/scripts/general_functions.R')
library("clusterProfiler")



wt.vs.rad21 <- as.data.frame(conds.gr.list$rad21_ko)[,c('entrez_id','log2FoldChange','padj')]
colnames(wt.vs.rad21) <- c('entrez_id','log2FC.rad21ko','padj.rad21ko')

wt.vs.ctcf <- as.data.frame(conds.gr.list$ctcf_ko)[,c('entrez_id','log2FoldChange','padj')]
colnames(wt.vs.ctcf) <- c('entrez_id','log2FC.ctcfko','padj.ctcfko')
#####

###
ego <- enrichGO(gene          = conds.gr.list$rad21_ko[conds.gr.list$rad21_ko$padj < 0.05 & conds.gr.list$rad21_ko$log2FoldChange > 0]$entrez_id,
                universe      = as.character(conds.gr.list$rad21_ko$entrez_id),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

dotplot(ego,showCategory=20)
###
test <- as.list(List(rad21ko.up = conds.gr.list$rad21_ko[conds.gr.list$rad21_ko$padj < 0.05 & conds.gr.list$rad21_ko$log2FoldChange > 0]$entrez_id,
     rad21ko.down = conds.gr.list$rad21_ko[conds.gr.list$rad21_ko$padj < 0.05 & conds.gr.list$rad21_ko$log2FoldChange < 0]$entrez_id,
     ctcfko.up = conds.gr.list$ctcf_ko[conds.gr.list$ctcf_ko$padj < 0.05 & conds.gr.list$ctcf_ko$log2FoldChange > 0]$entrez_id,
     ctcfko.down = conds.gr.list$ctcf_ko[conds.gr.list$ctcf_ko$padj < 0.05 & conds.gr.list$ctcf_ko$log2FoldChange < 0]$entrez_id)
)

uni=unique(c(as.character(conds.gr.list$rad21_ko$entrez_id),as.character(conds.gr.list$ctcf_ko$entrez_id,universe=uni)))

ck <- compareCluster(geneCluster = test, fun = "enrichGO",OrgDb=org.Mm.eg.db,ont= "BP",universe=uni)

dotplot(ck, showCategory=15)

######3
#######

uni=unique(c(as.character(conds.gr.list$rad21_ko$entrez_id),as.character(conds.gr.list$ctcf_ko$entrez_id,universe=uni),as.character(conds.gr.list$dev$entrez_id,universe=uni)))

test <- as.list(List(rad21ko.up = conds.gr.list$rad21_ko[conds.gr.list$rad21_ko$padj < 0.05 & conds.gr.list$rad21_ko$log2FoldChange > 0]$entrez_id,
                     rad21ko.down = conds.gr.list$rad21_ko[conds.gr.list$rad21_ko$padj < 0.05 & conds.gr.list$rad21_ko$log2FoldChange < 0]$entrez_id,
                     ctcfko.up = conds.gr.list$ctcf_ko[conds.gr.list$ctcf_ko$padj < 0.05 & conds.gr.list$ctcf_ko$log2FoldChange > 0]$entrez_id,
                     ctcfko.down = conds.gr.list$ctcf_ko[conds.gr.list$ctcf_ko$padj < 0.05 & conds.gr.list$ctcf_ko$log2FoldChange < 0]$entrez_id,
                     dev.up = conds.gr.list$dev[conds.gr.list$dev$padj < 0.05 & conds.gr.list$dev$log2FoldChange > 0]$entrez_id,
                     dev.down = conds.gr.list$dev[conds.gr.list$dev$padj < 0.05 & conds.gr.list$dev$log2FoldChange < 0]$entrez_id))

ck <- compareCluster(geneCluster = test, fun = "enrichGO",OrgDb=org.Mm.eg.db,ont= "BP",universe=uni)

dotplot(ck, showCategory=12)


