library(ELMER)
library(parallel)
library(MultiAssayExperiment)
distal.probes <- get.feature.probe(genome = "hg19", met.platform = "EPIC",promoter = T)
?createMAE
#save(IPF.met,dataFilt,inter,IPF.met.dmc,DEG_expr_rmdp,dataDEGs,file="TCGAbiolinks.RData")
load("TCGAbiolinks.RData")
#IPF.met需要用beta value
#DEG_expr_rmdp需要用log(count)
mae <- createMAE(exp = logexpr[,inter], 
                 met = IPF.met[,inter],
                 colData=data.frame(primary=inter,group=c(rep("IPF",23),rep("control",13))),
                 save = TRUE,
                 linearize.exp = F,  #不要在对RNA expr取log了
                 filter.probes = distal.probes,
                 save.filename = "mae_IPF.rda",
                 met.platform = "EPIC",
                 genome = "hg19",
                 TCGA = F)

group.col <- "group"
group1 <-  "IPF"
group2 <- "control"
#can be "hypo", "hyper" or "both", 
direction <- "both"
dir.out <- file.path("ELMER",direction)
dir.create(dir.out, recursive = TRUE)

sig.diff <- get.diff.meth(data = mae, 
                          group.col = group.col,
                          group1 =  group1,
                          group2 = group2,
                          minSubgroupFrac = 0.2,
                          sig.dif = 0.1,  #0.25
                          diff.dir = direction, # Search for hypomethylated probes in group1
                          cores = 2, 
                          dir.out = dir.out, 
                          pvalue = 0.05   #0.05
)

#############################################################
#获取显著的甲基化探针和基因对
#############################################################
# Collect nearby 20 genes for Sig.probes

nearGenes <- GetNearGenes(data = mae, 
                          probes = sig.diff$probe, 
                          numFlankingGenes = 20, #20， 10 upstream and 10 dowstream genes
                          #cores = 1
)

#############################################################
#利用散点图来展示甲基化水平和基因表达水平
##############################################################
scatter.plot(data = mae,
             byProbe = list(probe = sig.diff$probe[2], numFlankingGenes = 20), 
             category = "group", 
             dir.out = "plots",
             lm = TRUE, # Draw linear regression curve
             save = TRUE) 
sig.diff$probe[1]
# #"cg16269199"
# nearGenes <- GetNearGenes(data = mae, 
#                           probes = sig.diff$probe[1], 
#                           numFlankingGenes = 20, #20， 10 upstream and 10 dowstream genes
#                           #cores = 1
# )
# symbol <- getSymbol(mae, geneID = "ENSG00000187583")
# #PLEKHN1
# probe=sig.diff$probe[1]
# data=mae
# gene="ENSG00000187583"
# plot(assay(getMet(data)[probe, ]), 
#      assay(getExp(data)[gene, ]))


# scatter(meth = assay(getMet(data)[probe, ]), 
#              exp = assay(getExp(data)[gene, ]), 
#         category = category, 
#              ylim = ylim, dots.size = dots.size, legend.title = legend.title, 
#              correlation = correlation, xlab = sprintf("DNA methylation at %s", 
#                                                        probe), ylab = sprintf("%s gene expression", 
#                                                                               symbol), title = sprintf("%s_%s", probe, 
#                                                                                                       symbol), ...)
# 1 cg16269199 ENSG00000187583 PLEKHN1    891026 L1   
# 2 cg16269199 ENSG00000188290 HES4       923491 R1   
# 3 cg16269199 ENSG00000187608 ISG15      937952 R2   
# 4 cg16269199 ENSG00000162572 SCNN1D    1204965 R3   
# 5 cg16269199 ENSG00000127054 CPSF3L    1236114 R4   
# 6 cg16269199 ENSG00000169962 TAS1R3    1255843 R5   
# 7 cg16269199 ENSG00000162576 MXRA8     1277218 R6   
# 8 cg16269199 ENSG00000221978 CCNL2     1310240 R7   
# 9 cg16269199 ENSG00000242485 MRPL20    1326437 R8   
# 10 cg16269199 ENSG00000235098 ANKRD65   1342949 R9   
# 11 cg16269199 ENSG00000248333 CDK11B    1559752 R10  
# 12 cg16269199 ENSG00000008128 CDK11A    1623318 R11  
# 13 cg16269199 ENSG00000215790 SLC35E2   1645426 R12  
# 14 cg16269199 ENSG00000162585 C1orf86   2105052 R13  
# 15 cg16269199 ENSG00000116151 MORN1     2241841 R14  
# 16 cg16269199 ENSG00000157916 RER1      2312416 R15  
# 17 cg16269199 ENSG00000149527 PLCH2     2346568 R16  
# 18 cg16269199 ENSG00000215912 TTC34     2556564 R17  
# 19 cg16269199 ENSG00000130762 ARHGEF16  3360139 R18  
# 20 cg16269199 ENSG00000078900 TP73      3558233 R19
#################################################################
######################################################
# scatter.plot(data = mae,
#              byProbe = list(probe = "cg00100948", numFlankingGenes = 20), 
#              category = "definition", 
#              dir.out = "plots",
#              lm = TRUE, # Draw linear regression curve
#              save = TRUE) 
#######################################################
#批量画甲基化探针跟附近20个基因表达水平的关系
#######################################################
starburst=read.csv("starburst_results.csv",stringsAsFactors = F)
starburst_sig=starburst[starburst$starburst.status=="Down regulated & Hyper methylated" |starburst$starburst.status=="Up regulated & Hypo methylated", ]
for (cpg in starburst_sig$probeID){
  scatter.plot(data = mae,
               byProbe = list(probe = cpg, numFlankingGenes = 20), 
               category = "group", 
               dir.out = "probe_nearby20_plots",
               lm = TRUE, # Draw linear regression curve
               save = TRUE) 
}

for (i in 1:nrow(starburst_sig)){
  scatter.plot(data = mae,
               byPair  = list(probe = starburst_sig[i,1], gene=starburst_sig[i,3]), 
               category = "group", 
               dir.out = "pairs_plots",
               lm = TRUE, # Draw linear regression curve
               save = TRUE) 
}
########################################################
#批量画出受甲基化影响的基因的表达及甲基化水平
########################################################
dir.create("cpg_expression")
for(i in 1:nrow(starburst_sig)){
    #i=1
    gene=as.character(starburst_sig[i,3])
    symbol=getSymbol(mae,gene)
  pdf(file=paste("cpg_expression/",as.character(starburst_sig[i,1]),"_",symbol,".pdf",sep=""),width=9)
  par(mfrow=c(1,2))
  meth=assay(mae[[1]])[as.character(starburst_sig[i,1]),]
  type=factor(colData(mae)$group)
  ddf=data.frame(meth=meth,type=type)
  bp=boxplot(meth ~ type, data = ddf, lwd = 2, ylab = 'Beta value',main=as.character(starburst_sig[i,1]),xlab=NULL
             #lines(c(1,2),c(3,3),lwd = 2
  )
  stripchart(meth ~ type, vertical = TRUE, data = ddf, 
             method = "jitter", add = TRUE, pch = 20, col = c('blue','brown'))
  #text(1.5,3.4,"***",col="red",cex = 2)
  
  
  expr=assay(mae[[2]])[as.character(starburst_sig[i,3]),]
  type=factor(colData(mae)$group)
  ddf=data.frame(expr=expr,type=type)
  bp=boxplot(expr ~ type, data = ddf, lwd = 2, ylab = 'Expression',main=symbol,xlab=NULL
             #lines(c(1,2),c(3,3),lwd = 2
  )
  stripchart(expr ~ type, vertical = TRUE, data = ddf, 
             method = "jitter", add = TRUE, pch = 20, col = c('blue','brown'))
  #text(1.5,3.4,"***",col="red",cex = 2)
  dev.off()
}
