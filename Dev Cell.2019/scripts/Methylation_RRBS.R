library(methylKit)
list.files()
file.list <- list("cpg.Mm_HepKOCon_plus1.mincov10.txt",    
                  "cpg.Mm_HepKOCon_plus4.mincov10.txt",   
                  "cpg.Mm_LivUhrf1Fl_flCon_1.mincov10.txt",
                  "cpg.Mm_LivUhrf1Fl_flCon_2.mincov10.txt" )
myobj <- read(file.list, sample.id=list('KO','KO','Con','Con'),
              assembly='mm10', treatment=c(1,1,0,0), context='CpG')
### discriptive statistics
getMethylationStats(myobj[[2]],plot=T,both.strands=F)
getCoverageStats(myobj[[2]],plot=T,both.strands=F)
meth <- unite(myobj, destrand=FALSE)
getCorrelation(meth,plot=T)
### Finding DM bases or regions
myDiff <- calculateDiffMeth(meth)
# get all differentially methylated bases (MD larger than 25%)
myDiff25p=get.methylDiff(myDiff,difference=25,qvalue=0.05)
# get hyper methylated bases 
myDiff25p.hyper=get.methylDiff(myDiff,difference=25,qvalue=0.05,type="hyper")
# get hypo methylated bases
myDiff25p.hypo=get.methylDiff(myDiff,difference=25,qvalue=0.05,type="hypo")
### read gene annotation file
library(genomation)
gene.obj <- readTranscriptFeatures("mm10_refseq.bed")
cpg.obj <- readFeatureFlank("CpG annotation.txt",
                            feature.flank.name = c("CpGi", "shores"))

diffAnn = annotateWithGenicParts(myDiff25p, gene.obj)
diffCpGann = annotateWithFeatureFlank(myDiff25p,cpg.obj$CpGi, 
                                      cpg.obj$shores, feature.name = "CpGi", flank.name = "shores")
# plot percentages of differentially methylated bases over all chromosomes
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.05,meth.cutoff=25)
plotTargetAnnotation()
plotTargetAnnotation(diffCpGann, col=c('coral3','cadetblue3','cornsilk1'))
# get bedgraph file for differentially methylated bases
# ready to upload to UCSC
bedgraph(myDiff25p,file.name="KO_Con_diff.bedgraph", col.name="meth.diff")
save.image('KO_Con.RData')

########## Retained Methylation ##########
### plot those CpGs are retained
no_change <- myDiff[abs(myDiff$meth.diff) <= 0.05, ]
noAnn <- annotateWithGeneParts(as(no_change,"GRanges"), gene.obj)
## colors order: pro, exon, intron, intergenic
plotTargetAnnotation(noAnn,precedence=TRUE, bty="n",
                     col=c('brown3', 'darkolivegreen3','goldenrod3','azure3'),
                     main="Retained methylation annotation WT_KO")

noCpGann <- annotateWithFeatureFlank(as(no_change,"GRanges"),
                                     cpg.obj$CpGi,cpg.obj$shores,
                                     feature.name="CpGi",flank.name="shores")
## colors order: CpGi, shores, other
plotTargetAnnotation(noCpGann,col=c('brown3','cyan3',"azure3"),
                     main="retained methylation Ann WT_KO")

########### Residual Methylation ##########
### copy meth varialb
meth_bi <- meth
### caculate meth in con and ko, and meth diff
meth_bi$Con <- (meth_bi$numCs3 + meth_bi$numCs4)/(meth_bi$coverage3 + meth_bi$coverage4)
meth_bi$Ko <- (meth_bi$numCs1 + meth_bi$numCs2)/(meth_bi$coverage1 + meth_bi$coverage2)
meth_bi$diff <- meth_bi$Con - meth_bi$Ko
### ratio
length(which(meth_bi$Con < 0.2))/nrow(meth_bi)
### get residual methylated CpG sites, which are >80% methlated in KO
meth_RMD <- meth_bi[which(meth_bi$Ko > 0.8), ]
myDiff_RMD <- myDiff[which(meth_bi$Ko > 0.8), ]
### Annotation
library(genomation)
gene.obj <- readTranscriptFeatures("mm10_refseq.bed")
diffAnn_RMD <- annotateWithGeneParts(as(myDiff_RMD, 'GRanges'), gene.obj)
plotTargetAnnotation(diffAnn_RMD)
### write out
write.csv(meth_RMD, 'Residual_Methyl_CpG.csv', quote = F)
### overlap RMD with TE
TE <- readBed('mm10_RepeatMasker.bed')
findOverlaps(as(meth_RMD, 'GRanges'), TE)
TE_RMD_idx <- findOverlaps(TE, as(meth_RMD, 'GRanges'))
TE_RMD_idx <- unique(queryHits(TE_RMD_idx))
TE_RMD <- TE[TE_RMD_idx, ]$name
write.table(TE_RMD, 'Residual_TE.txt', quote = F, sep='\t')



