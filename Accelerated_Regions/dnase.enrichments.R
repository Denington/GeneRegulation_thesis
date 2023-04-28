#more specific dnase questions- what tissues are enriched... etc
#takes ~25 human/ 25 mouse, ask which are most active in accs compared to cnes
#do for gcbias only as well- interesting results

#ok- remaining- remove duplicates, and remove some more human cell lines
setwd('~/Documents/Projects/Accelerated_Regions/Acceleration_DNase_work/')

#use the collated matrices here, which are already nicely set as all cell lines/tissues overlapping all cnes
# -x (e.g. "cerebellar.cortex-adult" ) in each string says if are fetal/embryonic, adult, postnatal, or other ('unknown')
#so can subset easily
#may remove all postnatal
library(LOLA)
library(data.table)
library(rtracklayer)
library(ComplexHeatmap)

species_v <- c(hg19='Human',panTro4='Chimp',rheMac3='Macaque',ponAbe2='Orangutan',mm10='Mouse',rn5='Rat',cavPor3='Guinea Pig',oryCun2='Rabbit',susScr3='Pig',
               bosTau7='Cow',turTru2='Dolphin',orcOrc1='Killer Whale',canFam3='Dog',felCat5='Cat',lepWed1='Seal',odoRosDiv1='Walrus',equCab2='Horse',myoLuc2='Brown Bat',
               pteVam1='Vampire Bat',loxAfr3='Elephant',triMan1='Manatee')


mouse.dnase <- fread('mouse_collated.matrix.any.bed'); mouse.dnase <- makeGRangesFromDataFrame(mouse.dnase,keep.extra.columns = T)
human.dnase <- fread('human.collated.matrix.any.bed'); human.dnase <- makeGRangesFromDataFrame(human.dnase,keep.extra.columns = T)
#i can't remember what the differences is stringent vs not. probably a minimum score of 0.5?
#i should look on ax4 for how i strctured this, i think it's just average dnase accesibility score between 2 replicates

cnes.gr <- import.bed('~/Documents/Projects/Accelerated_Regions/CNEs.22.bed')
mouse.dnase <- mouse.dnase[unique(queryHits(findOverlaps(mouse.dnase,cnes.gr)))]
human.dnase <- human.dnase[unique(queryHits(findOverlaps(human.dnase,cnes.gr)))]


hg.acc <- import.bed('../generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/hg19.bed')
mm.acc <- import.bed('../generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/mm10.bed')

human.dnase.hgacc <- human.dnase[unique(queryHits(findOverlaps(human.dnase,hg.acc)))]
human.dnase.cne <- human.dnase[-unique(queryHits(findOverlaps(human.dnase,hg.acc)))]

mouse.dnase.mmacc <- mouse.dnase[unique(queryHits(findOverlaps(mouse.dnase,mm.acc)))]
mouse.dnase.cnes <- mouse.dnase[-unique(queryHits(findOverlaps(mouse.dnase,mm.acc)))]


get.enrichment <- function(bed,accs,mark,min=0.5){
  ac <- bed[unique(queryHits(findOverlaps(bed,accs)))]
  no <- bed[-unique(queryHits(findOverlaps(bed,accs)))]
  
  ac.df <- as.data.frame(ac)
  acc.l <- length(ac.df[,mark][ac.df[,mark] > min])
  acc.nl <- length(ac) - acc.l

  no.df <- as.data.frame(no)
  no.l <- length(no.df[,mark][no.df[,mark] > min])
  no.nl <- length(no) - no.l
  m <- matrix(c(acc.l,acc.nl,no.l,no.nl),ncol = 2)
  return(m)
}
#there's an argument to using cnes accelerated in none as the other in get.enrichment.. maybe?
#e.g get.enrichment(human.dnase,hg.acc,'midbrain.adult') - note it's '.' not '-', gets converted by R

all.hg.marks <- function(mark){
  out <- fisher.test(get.enrichment(human.dnase,hg.acc,mark))
  return(out$estimate)
}

hg.enrichments <- as.data.frame(sapply(colnames(as.data.frame(human.dnase))[-c(1,2,3,4,5)], all.hg.marks)) #this works

accel_beds <- Sys.glob('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/*.bed') #GC_bias/noGCbias/*bed')
for(i in accel_beds){
  sp_accel <- import.bed(i)
  sp <- strsplit(i,"\\/")
  sp <- sp[[1]][length(sp[[1]])]
  sp <- strsplit(sp,"\\.")[[1]][1]
  
  all.hg.marks.i <- function(mark){
    out <- fisher.test(get.enrichment(human.dnase,sp_accel,mark))
    return(out$estimate)
  }
  activemarks <- sapply(colnames(as.data.frame(human.dnase))[-c(1,2,3,4,5)], all.hg.marks.i) 
  hg.enrichments[,sp] <- activemarks
}

hg.enrichments.bu <- hg.enrichments #backup
hg.enrichments <- hg.enrichments[,-1]
#hg.enrichments <- hg.enrichments[,-c(6,12,18,22,25)] #if using the all species, this removes internal nodes
hg.enrichments <- hg.enrichments[ !row.names(hg.enrichments) %in% c("smooth.muscle.cell.of.the.brain.vasculature.unknown.odds ratio","PC.unknown.odds ratio" ,"NCI.unknown.odds ratio",
                  "dermis.blood.vessel.endothelial.cell.unknown.odds ratio","HCT116.unknown.odds ratio","fibroblast.of.dermis.unknown.odds ratio","EL.unknown.odds ratio",
                  "dermis.blood.vessel.endothelial.cell.postnatal.odds ratio","SK.adult.odds ratio",
                  "foreskin.fibroblast.postnatal.odds ratio",
                  "lung.adult.odds ratio",
                  "lung.unknown.odds ratio",
                  "spinal.cord.unknown.odds ratio"),]
rownames(hg.enrichments) <- stringr::str_remove(rownames(hg.enrichments),'\\.odds ratio')

hg.enrichments <- hg.enrichments[!duplicated(row.names(hg.enrichments)),]


#row.order <- species_v[c("hg19","panTro4","rheMac3","ponAbe2","mm10","rn5","cavPor3","oryCun2","susScr3","bosTau7",
#               "turTru2","orcOrc1","equCab2","lepWed1","odoRosDiv1","canFam3","felCat5","myoLuc2","pteVam1",
#               "triMan1","loxAfr3")]
row.order <- species_v[c("hg19","panTro4","rheMac3","mm10","rn5","cavPor3","susScr3","bosTau7",
                         "turTru2","equCab2","canFam3","felCat5","myoLuc2","loxAfr3")]

ctest <- rep('Cell line',length(rownames(hg.enrichments)))
ctest[grepl('fetal',rownames(hg.enrichments))] <- 'Embryonic'
ctest[grepl('postnatal',rownames(hg.enrichments))] <- 'Postnatal'
ctest[grepl('adult',rownames(hg.enrichments))] <- 'Adult'

rownames(hg.enrichments) <- stringr::str_replace(rownames(hg.enrichments), '\\.unknown','\\.')
rownames(hg.enrichments) <- stringr::str_replace(rownames(hg.enrichments), '\\.fetal','\\,')
rownames(hg.enrichments) <- stringr::str_replace(rownames(hg.enrichments), '\\.adult','\\`')
rownames(hg.enrichments) <- stringr::str_replace(rownames(hg.enrichments), '\\.postnatal',';')

rownames(hg.enrichments)[8] <- "cerebellum.astrocyte." ; rownames(hg.enrichments)[9] <- "hippocampus.astrocyte."
rownames(hg.enrichments)[10] <- "spinal.cord.astrocyte." ; rownames(hg.enrichments)[14] <- "brain.endothelial.cell."


colnames(hg.enrichments) <- species_v[colnames(hg.enrichments)]
ctest <- factor(ctest,levels = c('Embryonic','Postnatal','Adult','Cell line'))
pS <- Heatmap(t(scale(hg.enrichments)),rect_gp = gpar(col = "white", lwd = 1),row_order = row.order,
        column_split = ctest,column_gap = unit(0, "mm"), border = TRUE,cluster_column_slices = F) #maybe scale the data first
p <- Heatmap(t(hg.enrichments),rect_gp = gpar(col = "white", lwd = 1),row_order = row.order,
        column_split = ctest,column_gap = unit(0, "mm"), border = TRUE,cluster_column_slices = F) #maybe scale the data first

#column_names_rot = 45
##
all.mm.marks <- function(mark){
  out <- fisher.test(get.enrichment(mouse.dnase,mm.acc,mark))
  return(out$estimate)
}

mouse.enrichments <- as.data.frame(sapply(colnames(as.data.frame(mouse.dnase))[-c(1,2,3,4,5)], all.mm.marks)) #this works
accel_beds <- Sys.glob('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/*.bed') #GC_bias/noGCbias/*bed')
for(i in accel_beds){
  sp_accel <- import.bed(i)
  sp <- strsplit(i,"\\/")
  sp <- sp[[1]][length(sp[[1]])]
  sp <- strsplit(sp,"\\.")[[1]][1]
  
  all.mm.marks.i <- function(mark){
    out <- fisher.test(get.enrichment(mouse.dnase,sp_accel,mark))
    return(out$estimate)
  }
  activemarks <- sapply(colnames(as.data.frame(mouse.dnase))[-c(1,2,3,4,5)], all.mm.marks.i) 
  mouse.enrichments[,sp] <- activemarks
}

mouse.enrichments <- mouse.enrichments[,-1]
mouse.enrichments.bu <- mouse.enrichments #backup

rownames(mouse.enrichments) <- stringr::str_remove(rownames(mouse.enrichments),'\\.odds ratio')
#mouse.enrichments <- mouse.enrichments[,-c(6,12,18,22,25)]
mouse.enrichments <- mouse.enrichments[ !row.names(mouse.enrichments) %in% c("A20.unknown","G1E.ER4.postnatal" ,"G1E.postnatal",
                                                                             "cKit.positive.CD71.negative.TER119.negative.erythroid.progenitor.cells.embryonic",
                                                                             "cKit.positive.CD71.positive.TER119.negative.erythroid.progenitor.cells.embryonic",
                                                                             "cKit.positive.CD71.positive.TER119.positive.erythroid.progenitor.cells.embryonic"),]
mouse.enrichments <- mouse.enrichments[!duplicated(row.names(mouse.enrichments)),]

ctestM <- rep('Cell line',length(rownames(mouse.enrichments)))
ctestM[grepl('embryonic',rownames(mouse.enrichments))] <- 'Embryonic'
ctestM[grepl('postnatal',rownames(mouse.enrichments))] <- 'Postnatal'
ctestM[grepl('adult',rownames(mouse.enrichments))] <- 'Adult'
ctestM <- factor(ctestM,levels = c('Embryonic','Postnatal','Adult','Cell line'))
rownames(mouse.enrichments) <- stringr::str_replace(rownames(mouse.enrichments), '\\.unknown','\\.')
rownames(mouse.enrichments) <- stringr::str_replace(rownames(mouse.enrichments), '\\.embryonic','\\,')
rownames(mouse.enrichments) <- stringr::str_replace(rownames(mouse.enrichments), '\\.adult','\\`')
rownames(mouse.enrichments) <- stringr::str_replace(rownames(mouse.enrichments), '\\.postnatal','\\;')

#rownames(mouse.enrichments)[10] <- "CD4+CD25+.alpha.beta.T.cell`"
#rownames(mouse.enrichments)[11] <- "CD4+.helper.T.cell`" 
#rownames(mouse.enrichments)[14] <- "erythroid.progenitor.cells,"
#rownames(mouse.enrichments)[49] <-  "naive.CD4+.alpha.beta.T.cell`"

colnames(mouse.enrichments) <- species_v[colnames(mouse.enrichments)]


pSM <- Heatmap(t(scale(mouse.enrichments)),rect_gp = gpar(col = "white", lwd = 1),row_order = row.order,
        column_split = ctestM,column_gap = unit(0, "mm"), border = TRUE,cluster_column_slices = F) #maybe scale the data first
pM <- Heatmap(t(mouse.enrichments),rect_gp = gpar(col = "white", lwd = 1),row_order = row.order,
        column_split = ctestM,column_gap = unit(0, "mm"), border = TRUE,cluster_column_slices = F) #maybe scale the data first

#column_names_rot = 30
#works, but need to shear the names of the rows. Partly should use fewer types
#A20, G18

######general CNE enrichment over random (with regioneR... this is going to be slow)
"""library(regioneR)
#im tempted just to do this on accels so it's not so slow
general.enrichments <- function(bed,accs,mark,min=0.5,genome='hg19',nshuffles=250){
  
  bed.df <- as.data.frame(bed)
  mark.bed <- makeGRangesFromDataFrame(bed.df[bed.df[,mark] > 1,] )
  #

  #
  pt <- overlapPermTest(A= accs, B=mark.bed, genome=genome,ntimes = nshuffles,
                 force.parallel = T)
  oe<- pt$numOverlaps$observed / mean(pt$numOverlaps$permuted)
  return(oe)
}

all.hg.marks.general <- function(mark){
  out <- general.enrichments(human.dnase,human.dnase,mark)
  return(out)
}

hg.general.enrichments <- as.data.frame(sapply(colnames(as.data.frame(human.dnase))[-c(1,2,3,4,5)], all.hg.marks.general)) 
#is going to be v slow!
#don;t need this for every species....?
#for(i in accel_beds){
#  sp_accel <- import.bed(i)
#  sp <- strsplit(i,'\\/')
#  sp <- sp[[1]][length(sp[[1]])]
#  sp <- strsplit(sp,'\\.')[[1]][1]
all.mm.marks.general <- function(mark){
  out <- general.enrichments(mouse.dnase,mouse.dnase,mark)
  return(out)
}

mm.general.enrichments <- as.data.frame(sapply(colnames(as.data.frame(mouse.dnase))[-c(1,2,3,4,5)], all.mm.marks.general)) 

  
#  all.mm.marks.i <- function(mark){
#    out <- all.hg.marks.general(human.dnase,sp_accel,mark)
#    return(out$estimate)
#  }
#  activemarks <- sapply(colnames(as.data.frame(mouse.dnase))[-c(1,2,3,4,5)], all.mm.marks.i) 
#  hg.general.enrichments[,sp] <- activemarks
#}
""" 
#the above needs to be redone from scratch, over all beds rather
#mouse/rat immune related stuff may reflect observed e.g. unsual immune cell distributions
#should look at bgc vs non-bgc

###################################
################################### getting average values each cell line
#it would be nice to get the below in, but it doesn't feel like a make or break
library(DescTools)
quan <- function(x){return(quantile(x,0.9))} # or do a winsorized mean
w.mean <- function(x){
  z = Winsorize(x,probs = c(0,0.99))
  return(mean(z))
}
non0.mean <- function(x){
  z = Winsorize(x,probs = c(0,0.99))
  z = z[z != 0]
  return(mean(z))
}

human.cell.order <- unlist(column_order(p))
mouse.cell.order <- unlist(column_order(pM))
#
hg.hAd <- human.dnase[unique(queryHits(findOverlaps(human.dnase,hg.acc)))]
hg.mAd <- human.dnase[unique(queryHits(findOverlaps(human.dnase,mm.acc)))]
hg.noAd <- human.dnase[-unique(queryHits(findOverlaps(human.dnase,c(hg.hAd,hg.mAd))))]
#want to try the 'no' as not accel in any

hg.hAd <- as.data.frame(hg.hAd)
hg.mAd <- as.data.frame(hg.mAd)
hg.noAd <- as.data.frame(hg.noAd)

hg.hgav <- apply(hg.hAd[,-c(1:5)],2,non0.mean)
hg.mmav <- apply(hg.mAd[,-c(1:5)],2,non0.mean)
hg.noav <- apply(hg.noAd[,-c(1:5)],2,non0.mean)

hgav.df <- data.frame(cell.line=names(hg.hgav),av.ac=hg.hgav,element='HAR')
mmav.df <- data.frame(cell.line=names(hg.mmav),av.ac=hg.mmav,element='MAR')
noav.df <- data.frame(cell.line=names(hg.noav),av.ac=hg.noav,element='CNE')
av.act <- rbind(hgav.df,mmav.df,noav.df)
av.act$cell.line <- stringr::str_replace(av.act$cell.line, '\\.unknown','\\.'); 
av.act$cell.line <-  stringr::str_replace(av.act$cell.line, '\\.fetal','\\,'); 
av.act$cell.line <-  stringr::str_replace(av.act$cell.line,'\\.adult','\\`') ; 
av.act$cell.line <-  stringr::str_replace(av.act$cell.line, '\\.postnatal','\\;')

av.act$cell.line <- factor(av.act$cell.line,levels = (rownames(hg.enrichments)[human.cell.order]))

ggplot(av.act[!is.na(av.act$cell.line),],aes(cell.line,av.ac,color=element)) + geom_point() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('Averaged DNase over elements')

####same for mouse

mm.hAd <- mouse.dnase[unique(queryHits(findOverlaps(mouse.dnase,hg.acc)))]
mm.mAd <- mouse.dnase[unique(queryHits(findOverlaps(mouse.dnase,mm.acc)))]
mm.noAd <- mouse.dnase[-unique(queryHits(findOverlaps(mouse.dnase,c(mm.hAd,mm.mAd))))]
#want to try the 'no' as not accel in any

mm.hAd <- as.data.frame(mm.hAd)
mm.mAd <- as.data.frame(mm.mAd)
mm.noAd <- as.data.frame(mm.noAd)

mm.hgav <- apply(mm.hAd[,-c(1:5)],2,non0.mean)
mm.mmav <- apply(mm.mAd[,-c(1:5)],2,non0.mean)
mm.noav <- apply(mm.noAd[,-c(1:5)],2,non0.mean)

hgavM.df <- data.frame(cell.line=names(mm.hgav),av.ac=mm.hgav,element='HAR')
mmavM.df <- data.frame(cell.line=names(mm.mmav),av.ac=mm.mmav,element='MAR')
noavM.df <- data.frame(cell.line=names(mm.noav),av.ac=mm.noav,element='CNE')
av.act <- rbind(hgavM.df,mmavM.df,noavM.df)
av.act$cell.line <- stringr::str_replace(av.act$cell.line, '\\.unknown','\\.'); 
av.act$cell.line <-  stringr::str_replace(av.act$cell.line, '\\.embryonic','\\,'); 
av.act$cell.line <-  stringr::str_replace(av.act$cell.line,'\\.adult','\\`') ; 
av.act$cell.line <-  stringr::str_replace(av.act$cell.line, '\\.postnatal','\\;')


av.act$cell.line <- factor(av.act$cell.line,levels = (rownames(mouse.enrichments)[mouse.cell.order]))

ggplot(av.act[!is.na(av.act$cell.line),],aes(cell.line,av.ac,color=element)) + geom_point() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('Averaged DNase over elements')
