#just ideogram stuff
library(rtracklayer)
# need to 
x <- import.bed('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/hg19.AcceleratedRegions.bed')
y <- import.bed('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/mm10.AcceleratedRegions.bed')

library(karyoploteR)
kp <- plotKaryotype(genome="hg19",chromosomes =  "autosomal")
kp <- kpPlotDensity(kp,y,col=transparent('grey'),window.size = 1000000)
kp <- kpPlotDensity(kp,x,col=transparent('red'),window.size = 1000000)


#kp <- kpLines(kp, data=kp$latest.plot$computed.values$windows, y=kp$latest.plot$computed.values$density, col="black", r0=0.5, r1=1, data.panel=2, ymax=20)


##

library(RIdeogram)
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")

library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)

get.tiles <- function(bed,wid){
  gr.windows <- tileGenome(seqinfo(Hsapiens), tilewidth=wid,cut.last.tile.in.chrom=TRUE)
  gr.windows$val <- countOverlaps(gr.windows,import.bed(bed)) 
  a <- as.data.frame(gr.windows)
  a$Chr <- stringr::str_remove(a$seqnames,'chr')
  a <- a[,c(7,2,3,6)]
  colnames(a) <- c('Chr','Start','End','Value')
  return(a)
}

cne_tile <- get.tiles('~/Documents/Projects/Accelerated_Regions/CNEs.July17.bed',1000000)
mm_acc <- get.tiles('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/mm10.AcceleratedRegions.bed',1000000)
mm_acc$Color1 <- 'fc8d62'

hg_acc <- get.tiles('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/hg19-panTro4.AcceleratedRegions.bed',1000000)
hg_acc$Color <- '8da0cb'
mm_acc$Value2 <- hg_acc$Value
mm_acc$Color2 <- hg_acc$Color ; colnames(mm_acc) <- c('Chr','Start','End','Value_1','Color_1','Value_2','Color_2')

ideogram(karyotype = human_karyotype, overlaid = cne_tile,label=mm_acc,label_type='line')
convertSVG("chromosome.svg", device = "png")
