library(limma)
#BiocManager::install("minfi")
library(minfi)
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#BiocManager::install("IlluminaHumanMethylationEPICmanifest")
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
#BiocManager::install("missMethyl")
library(missMethyl)
library(matrixStats)
#BiocManager::install("minfiData")
library(minfiData)
library(Gviz)
#BiocManager::install("DMRcate")
library(DMRcate)
library(stringr)

targets=read.table("methy_sample_type.txt",header=T,sep="\t",row.names=1)
targets$type=factor(c(rep("IPF",24),rep("control",13)))
#需要先解压出idat文件
rgset <- read.metharray.exp(file.path("GSE173356_RAW"))
colnames(rgset)=targets[gsub("_.*$","",colnames(rgset)),"name"]

# get the 850k annotation data
ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

detP <- detectionP(rgset)
head(detP)

# examine mean detection p-values across all samples to identify any failed samples
pdf(file="mean_detection_pvalue.pdf",width=7)
pal <- brewer.pal(8,"Dark2")
par(
    #mfrow=c(1,2),
    mar=c(9,4,4,2)
    )
# barplot(colMeans(detP), col=pal[factor(targets$type)], las=2, 
#         cex.names=0.8, ylab="Mean detection p-values",yaxt = "n")
# abline(h=0.05,col="red")
# axis(2)
# legend("topleft", legend=levels(factor(targets$type)), fill=pal,
#        bg="white")

barplot(colMeans(detP), col=pal[factor(targets$type)], las=2, 
        cex.names=0.8, ylim=c(0,0.0016), ylab="Mean detection p-values",yaxt = "n")
abline(h=0.05,col="red")
axis(2)
legend("topright", legend=levels(factor(targets$type)), fill=pal, 
       bg="white",
       bty="n"
       )
dev.off()

##############################
#
#############################


## Not run: 

qcReport(rgset,sampNames=targets$name, sampGroups=targets$type, 
         pdf="qcReport.pdf")


###################################################
#去掉不合格的样本
###################################################
# remove poor quality samples
keep <- colMeans(detP) < 0.05
RGset <- rgset[,keep]

targets <- targets[keep,]
#targets[,1:2]

detP <- detP[,keep]
dim(detP)

#################################################################
#归一化
#################################################################
# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(RGset)

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgset)


# visualise what the data looks like before and after normalisation
pdf(file="before_after_norm.pdf",width=15)
par(mfrow=c(1,2))
densityPlot(RGset, sampGroups=targets$type,main="Raw", legend=FALSE
)
legend("top", legend = levels(factor(targets$type)),
       text.col=brewer.pal(8,"Dark2"))

densityPlot(getBeta(mSetSq), sampGroups=targets$type,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$type)),
       text.col=brewer.pal(8,"Dark2"))
dev.off()
#######################################################################
#去掉不合格的探针
#######################################################################
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

keep <- !(featureNames(mSetSqFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

#install.packages("devtools")
library(devtools)
#devtools::install_github("markgene/maxprobes")
#移除匹配在多个基因组上的探针（如果是850k）
library(maxprobes) 
xloci <- maxprobes::xreactive_probes(array_type = "EPIC")
length(xloci)
mSetSqFlt <- maxprobes::dropXreactiveLoci(mSetSqFlt) 
#########################################################
#PCA
#########################################################

#par(mfrow=c(1,2))
pdf(file="PCA_P1_P2.pdf",width=10)
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$type)], cex=0.8)
legend("top", legend=levels(factor(targets$type)), text.col=pal,
       cex=0.65, bg="white")
dev.off()

pdf(file="PCA_others.pdf",width=15)
par(mfrow=c(1,3))
# Examine higher dimensions to look at other sources of variation
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$type)], dim=c(1,3))
legend("top", legend=levels(factor(targets$type)), text.col=pal,
       cex=0.7, bg="white")
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$type)], dim=c(2,3))
legend("topright", legend=levels(factor(targets$type)), text.col=pal,
       cex=0.7, bg="white")
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$type)], dim=c(3,4))
legend("right", legend=levels(factor(targets$type)), text.col=pal,
       cex=0.7, bg="white")
dev.off()
###############################################
#计算M值和beta值
###############################################
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])

bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])

pdf(file="mvalue_betavalue.pdf",width=15)
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$type, main="Beta values",
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$type)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals,sampGroups=targets$type, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$type)),
       text.col=brewer.pal(8,"Dark2"))
dev.off()


###########################################################
#DMP
###########################################################
# this is the factor of interest
#cellType <- factor(targets$type)
# this is the individual effect that we need to account for
#individual <- factor(targets$Sample_Source)

# use the above to create a design matrix
design <- model.matrix(~0+type,data=targets)
colnames(design) <- levels(targets$type)

# fit the linear model
fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(IPF-control,levels=design)

contMatrix
# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

# get the table of results for the first contrast (naive - rTreg)
ann850kSub <- ann850k[match(rownames(mVals),ann850k$Name),
                      c(1:4,12:19,24:ncol(ann850k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann850kSub)
head(DMPs)
dim(DMPs)
write.table(DMPs[DMPs$adj.P.Val<0.01,], file="case_DMPs.csv", sep=",", row.names=FALSE)
#####################################
#volcano plot
#####################################
library(ggplot2)
VolcanoPlot <- function(deg.all, fc=2, pval=0.01) {
    geneList <- deg.all
    geneList$threshold <- c()
    geneList$threshold[geneList$logFC>log(fc,2) & geneList$FDR<pval] <- 1
    geneList$threshold[geneList$logFC>=-log(fc,2) & geneList$logFC<=log(fc,2) 
                       | geneList$FDR>=pval] <- 2
    geneList$threshold[geneList$logFC < -log(fc,2) & geneList$FDR<pval] <- 3
    
    geneList$threshold <- as.factor(geneList$threshold)
    
    lim <- max(max(geneList$logFC), abs(min(geneList$logFC)))+0.5
    
    volcano <- ggplot(data=geneList, aes(x=geneList$logFC, 
                                         y = -log10(geneList$FDR)))
    volcano+geom_point(aes(color=geneList$threshold), alpha=1, size=0.8) + 
        xlab("log2(Fold Change)") + ylab("-log10(FDR)") +
        scale_colour_manual(breaks = geneList$threshold, 
                            values = c('red','black','green3')) + xlim(c(-lim,lim)) +
        geom_vline(xintercept = c(-log(fc,2),log(fc,2)), 
                   color='darkgreen', linetype=3) + 
        geom_hline(yintercept = -log(pval,10), color='darkgreen',linetype=3)+
        theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour='black'),
              panel.background = element_blank()) +
        theme(legend.position="none") +
        theme(axis.text=element_text(size=14), 
              axis.title=element_text(size=16))
}

names(DMPs)[40]="FDR"
pdf(file="IPF_DMP_volcano.pdf")
VolcanoPlot(DMPs,fc=1.5,pval=0.05)
dev.off()
#####################################
#heatmap
######################################
heatmap(mVals[rownames(DMPs[DMPs$FDR<0.01,]),])
##############################################
#################################################
###热图 根据差异表达RNA对样本进行聚类
##################################################
library(gplots)
pdf(file="IPF_heatmap_0.01.pdf",width=12,height=7)
DMPsig=mVals[rownames(DMPs[DMPs$FDR<0.01,]),]
degName = rownames(DMPsig)
deg.id = degName 
#rna.expr = DMPsig
degDa <- DMPsig
sampleCol <- ifelse(targets$type == "IPF", 
                    "red", "blue")
########################
#labRow控制行名
#lmat控制四块的摆放顺序
#lwid控制宽度
#lhei控制高度
#margins控制列名和行名的留白大小
########################
lmat = rbind(c(4, 3), c(2, 1))
lwid = c(1, 5)
lhei = c(1, 5)
heatmap.2(as.matrix(degDa),col = bluered(75), trace = "none", 
          cexCol = 0.32, cexRow = 1, dendrogram = "both", srtCol = 90, 
          adjCol = c(0.8, 0.15), density.info = "none", labRow = NA, 
          key.title = NA, na.color = NA, lwid = lwid, lhei = lhei, 
          margins = c(1, 10), labCol = NA, key.xlab = "Normalized intensity", 
          scale = "row", ColSideColors = sampleCol)
#####
#可以得到鼠标在图上的位置，然后根据这个位置摆放legend的位置
#coords <- locator(1)
#
legend(y=0.46, x=0.95, xpd=TRUE,     
       legend = c("IPF","control"),
       col = c("red","blue"), 
       #lty= 1,             
       #lwd = 5,
       pch=15 ,          
       cex=1
)
dev.off()
#####################################
#
#####################################
pdf(file="top4_DMC.pdf",width=10)
# plot the top 4 most significantly differentially methylated CpGs
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
    plotCpg(bVals, cpg=cpg, pheno=targets$type, ylab = "Beta values")
})
dev.off()
#################################################
#ggplot top DMPs
#################################################
#bVals[rownames(DMPs)[1:4],]
#DMPs[rownames(DMPs)[1:4],]
library(ggplot2)
num=20
nrow=ceiling(num/5)
data=data.frame(t(bVals[rownames(DMPs)[1:num],]),group=targets$type)
expr_gather <- tidyr::gather(data,key = probe, value = beta,-c(group))
pdf(file=paste("top",num,"_DMP.pdf",sep=""),width=10)
ggplot(
    expr_gather, 
    aes(x = group, y = beta,color=group)
) + 
    geom_boxplot()+scale_color_manual(values=c("blue", "red"))+
    facet_wrap(~probe,nrow = nrow) 
dev.off()
######################################
#DMR
######################################
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M",
                             analysis.type = "differential", design = design,
                             contrasts = TRUE, cont.matrix = contMatrix,
                             coef = "IPF - control", arraytype = "EPIC")
#Your contrast returned 190503 individually significant probes. 
#We recommend the default setting of pcutoff in dmrcate().


str(myAnnotation)
##########################################
#all DMRs
##########################################
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2,min.cpgs=2)
results=extractRanges(DMRs,genome="hg19")
#dim(DMRs$results)
write.table(file="DMR_all.txt",results,quote=F,sep="\t")

#############################################
#betacutoff>0.15
#############################################

DMRs.1=results[abs(results$meandiff)>0.1,]
write.table(file="DMR_0.1.txt",DMRs.1,quote=F,sep="\t",row.names=F)
#87 DMRs

# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$type))]
names(groups) <- levels(factor(targets$type))
cols <- groups[as.character(factor(targets$type))]
#names(cols)=rownames(targets)
samps <- 1:nrow(targets)

# draw the plot for the top DMR
#beta=getBeta(mSetSq)
# pdf(file="DMR1.pdf",width=10,height=10)
# DMR.plot(ranges=DMRs.1, dmr=1, CpGs=beta, phen.col=as.character(cols), what = "Beta",
#          arraytype = "EPIC",
#          genome="hg19")
# dev.off()


colnames(bVals)=rownames(targets)
for(i in 7:10){
    i=71
pdf(file=paste0("DMR_",i,".pdf"),width=10)
DMR.plot(ranges=DMRs.1, dmr=i, CpGs=bVals, phen.col=cols, what = "Beta",
         arraytype = "EPIC", 
         #toscale=TRUE, 
         #plotmedians=TRUE,
         genome="hg19", samps=samps)
dev.off()

}
#############################
#Gene Ont testing
#############################
# Get the significant CpG sites at less than 5% FDR
cutoff=0.05
sigCpGs <- DMPs$Name[DMPs$FDR<cutoff]
# First 10 significant CpGs
sigCpGs[1:10]

# Total number of significant CpGs at 5% FDR
length(sigCpGs)

# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name
# Total number of CpG sites tested
length(all)

par(mfrow=c(1,1))
go <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=T)
## Warning in alias2SymbolTable(flat$symbol): Multiple symbols ignored for one
## or more aliases

# Top 10 GO categories
go=go[order(go$FDR,decreasing =F),]
go=go[go$FDR<0.05,]
write.table(file=paste("GO_",cutoff,"_.txt",sep=""), go,quote=F,sep="\t")
#go=read.table("GO_1e-05_.txt",header=T,sep="\t",quote="")


kegg <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=F,collection="KEGG" )
kegg=kegg[order(kegg$FDR,decreasing =F),]
kegg=kegg[kegg$FDR<0.05,]
write.table(file=paste("KEGG_",cutoff,"_.txt",sep=""), kegg,quote=F,sep="\t")
#kegg=read.table("KEGG_1e-05_.txt",header=T,sep="\t",quote="")
##############################################
#GO KEGG plot
##############################################
library(clusterProfiler)
library(ggplot2)

enrichBubblePlotFun <- function(kegg, pval=0.01) {
    keggBasic = ggplot(data=kegg, mapping=aes(x=Terms, 
                                              y=fe,color=FDR,size=Count))
    keggBasic+geom_point()+ coord_flip() + 
        scale_x_discrete(limits=rev(kegg$Terms)) + 
        scale_colour_gradientn(limits=c(0,pval),
                               colors= c("red","yellow","green")) + 
        xlab('')+ylab('Fold enrichment') +
        guides(shape = guide_legend(order=1),
               colour = guide_colourbar(order=2)) + 
        theme_bw()+theme(axis.line = element_line(colour = "black"),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour='black'),
                         panel.background = element_blank()) +
        ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size=20)) +
        theme(axis.text.y=element_text(size=14), 
              axis.title=element_text(size=15),
              axis.text.x=element_text(size=14)) + 
        theme(legend.text = element_text(size = 14),
              legend.title = element_text(size = 14)) +
        theme(strip.text = element_text(size = 14), 
              legend.key.size = unit(0.8,'cm'))
}


enrichBarPlotFun <- function(kegg, type='single', bar.color='black') {
    
    if (type=='single') {
        keggBasic = ggplot(data=kegg, mapping=aes(x=Term, y=-log(FDR,10)))
        keggBasic + geom_bar(stat='identity',fill=bar.color) + 
            scale_x_discrete(limits=rev(kegg$Term)) + 
            ylim(0, max(-log(kegg$FDR,10))) + 
            theme(legend.title=element_blank())+ylab('-log10(FDR)')+
            xlab('') + coord_flip() + 
            theme_bw()+theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_rect(colour='white'),
                             panel.background = element_blank()) +
            theme(axis.text.y=element_text(size=14), 
                  axis.title=element_text(size=15),
                  axis.text.x=element_text(size=14)) + 
            theme(legend.position = 'none')
        
    } else if (type=='multiple') {
        keggBasic = ggplot(data=kegg, mapping=aes(x=Term, y=-log(FDR,10), 
                                                  fill=Ont))
        keggBasic + geom_bar(stat='identity') + 
            scale_x_discrete(limits=rev(kegg$Term)) + 
            ylim(0, max(-log(kegg$FDR,10))) + 
            theme(legend.title=element_blank())+
            ylab('-log10(FDR)')+xlab('') + coord_flip() + 
            scale_fill_hue(name='',breaks=kegg$Ont, 
                           labels=kegg$Ont) +
            theme_bw()+theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_rect(colour='white'),
                             panel.background = element_blank()) +
            theme(axis.text=element_text(size=14), 
                  axis.title=element_text(size=15)) + 
            theme(legend.text = element_text(size=16))
    }
}

goPlot<-function (enrichment, type = "bar", category = "GO", num.terms = 10, 
                  bar.color = "black") 
{
    goBP <- enrichment[enrichment$Ont == "BP", ]
    goCC <- enrichment[enrichment$Ont == "CC", ]
    goMF <- enrichment[enrichment$Ont == "MF", ]
    kegg <- enrichment
    #DO <- enrichment[enrichment$Category == "DO", ]
    if (nrow(goBP) > num.terms) {
        goBP <- goBP[seq_len(num.terms), ]
    }
    if (nrow(goCC) > num.terms) {
        goCC <- goCC[seq_len(num.terms), ]
    }
    if (nrow(goMF) > num.terms) {
        goMF <- goMF[seq_len(num.terms), ]
    }
    if (nrow(kegg) > num.terms) {
        kegg <- kegg[seq_len(num.terms),]
    }
    
    if (category == "GO") {
        go <- rbind(goBP, goCC, goMF)
        if (type == "bubble") {
            go$Terms <- paste(go$Term, "[", go$Ont, "]", 
                              sep = "")
            enrichBubblePlotFun(go)
        }
        else if (type == "bar") {
            enrichBarPlotFun(go, type = "multiple")
        }
    }
    else if (category %in% c("BP", "CC", "MF","KEGG")) {
        goList <- list(goBP, goCC, goMF,kegg)
        names(goList) <- c("BP", "CC", "MF","KEGG")
        go <- goList[[category]]
        if (type == "bubble") {
            go$Terms <- paste(go$Term, "[", go$Ont, "]", 
                              sep = "")
            enrichBubblePlotFun(go)
        }
        else if (type == "bar") {
            enrichBarPlotFun(go, bar.color = bar.color)
        }
    }
}
pdf(file="go_barplot.pdf",width=15)
names(go)[1]="Ont"
names(go)[2]="Term"
goPlot(go, type='bar', category='GO', num.terms = 10)
dev.off()
pdf(file="kegg_barplot.pdf",width=15)
names(kegg)[1]="Term"
goPlot(kegg,type='bar', category='KEGG', num.terms = 10,bar.color = "blue")
dev.off()
#################################################################################
#KEGG plot
##################################################################################
library(DOSE)
font.size=12
num.terms = 10

barplotkegg<-function(data,font.size=12,num.terms = 10){
    if(nrow(data)>=num.terms){
        data=data[1:num.terms,]
    }
    #p <- ggplot(data, aes(x = Term, y = -log(FDR,10)))
    p <- ggplot(data, aes(x = Term, y = N))
    p <- p + geom_bar(stat = "identity") + coord_flip() + theme_dose(font.size)
    p <- p + aes(fill=FDR) +
        scale_fill_continuous(low="red", high="blue")
    p <- p + ggtitle("KEGG Enrichment") + xlab("") + ylab("Counts")
    return(p)
}
pdf(file="kegg_barplot2.pdf",width=15)
barplotkegg(kegg)
dev.off()
