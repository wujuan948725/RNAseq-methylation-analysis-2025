load("PolyA_Counts2.RData")
library(edgeR)
ployA=read.table("GSE173355_PolyA_Counts.txt",header=T,sep=" ",check.names = F)
dim(ployA)
ployA[1:3,1:3]

ids=strsplit(ployA$target_id,"\\|")
id_df=data.frame(do.call(rbind,ids))
sum(id_df$X8=="protein_coding")
#[1] 80930
pc_index=id_df$X8=="protein_coding"
pc_expr=ployA[pc_index,-1]
rownames(pc_expr)=id_df[pc_index,1]
boxplot(log2(pc_expr+1))
length(unique(id_df[pc_index,1]))

RNAexpr = DGEList(counts = pc_expr)
keep <- rowSums(cpm(RNAexpr) > 1) >= 0.5*ncol(pc_expr)
RNAexpr <- RNAexpr[keep,,keep.lib.sizes = TRUE]

RNAexpr = calcNormFactors(RNAexpr)
RNAexpr <- voom(RNAexpr, design=NULL, plot = FALSE)$E
type=factor(c(rep("IPF",23),rep("control",14)))
pdf(file="normalized_mRNA_expression.pdf",width=12)
par(mar=c(10,4,4,2))
boxplot(RNAexpr,las=2,ylim=c(-10,20),
        col=c("blue","red")[type]
        )
legend("top",legend=(levels(type)),col=c("blue","red"),horiz = T,bty="n",pch=15)
dev.off()

#####################################################
#
#####################################################
#定义做差异表达分析的函数
#@counts是合并之后的表达矩阵
#@group是样本的类型，PrimaryTumor还是SolidTissueNormal
#@comparison是比较的方式，eg. PrimaryTumor-SolidTissueNormal,就是肿瘤比上癌旁正常对照
#@method，差异表达分析方法，默认为limma，也可以设定为edgeR或DESeq2
#@filter，逻辑值，是否根据基因的cpm值对基因进行过滤
DEA <- function(counts, group, comparison, method='limma', filter=TRUE) {
  #加载注释文件，这个文件中包含了基因的注释信息，有ensembl gene ID，基因名字以及基因的编码类型
  #是mRNA还是长非编码RNA
  load("annotation.rda")
  #基于表达谱矩阵创建一个DGEList的对象
  dge = DGEList(counts = counts)
  #计算每一个gene的cpm（count per million）值，要求在一半的样本里面cpm>1
  #这个基因才会被保留，否则这个基因会被删掉，不进行后续分析
  #这里的cpm>1，以及0.5的样本数，可以根据自己的需求修改
  keep <- rowSums(cpm(dge) > 1) >= 0.5*length(group)
  
  #判断使用哪一种方法做差异表达分析
  #如果使用DESeq2
  if (method == 'DESeq2') {
    #根据样本类型创建一个数据框
    coldata <- data.frame(group)
    #这个是创建DESeq2特有的数据格式，包括表达谱矩阵，样本类型，以及分组
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = coldata, design = ~ group)
    
    #dds对象中增加一个group分组信息，是因子
    #因子的顺序根据comparison的设置来确定
    dds$group <- factor(dds$group, levels = 
                          rev(strsplit(comparison, '-', fixed=TRUE)[[1]]))
    
    #如果设置了过滤，那么根据cpm值进行过滤
    #keep中保存的是通过了过滤标准的基因
    if (filter==TRUE) {
      dds <- dds[keep, ]
    }
    
    #使用DESeq函数做差异表达分析
    dds <- DESeq(dds)
    #提取差异表达分析的结果
    res <- results(dds)
    #转换成数据框
    DEGAll <- data.frame(res)
    #设置数据框的列明
    colnames(DEGAll) <- c('baseMean', 'logFC', 'lfcSE', 'stat', 
                          'PValue', 'FDR')
    #如果设置的差异表达分析方法为edgeR或者limma
  } else if (method %in% c('edgeR', 'limma')) {
    #将分组信息转换成因子
    group <- factor(group)
    #创建差异表达分析的设计，也就是如何对样本进行分组
    design <- model.matrix(~0+group)
    #修改这个设计的列名
    colnames(design) <- levels(group)
    #创建比较方法，eg.PrimaryTumor-SolidTissueNormal，肿瘤比上癌旁正常对照
    contrast.matrix <- makeContrasts(contrasts=comparison, 
                                     levels=design)
    #如果设置了过滤，那么根据cpm值进行过滤
    #keep中保存的是通过了过滤标准的基因
    if (filter==TRUE) {
      dge <- dge[keep,,keep.lib.sizes = TRUE]
    }
    #有基因被过滤之后，每个样本的文库大小会发生变化
    #这里需要重新计算归一化因子
    dge <- calcNormFactors(dge)
    
    #进一步区分，如果使用的是edgeR方法
    if (method == 'edgeR') {
      #先估计基因表达的离散程度
      dge <- estimateDisp(dge, design)
      #对每一个基因进行负二项式广义对数线性模型拟合
      fit <- glmFit(dge, design)
      lrt <- glmLRT(fit, contrast=contrast.matrix)
      #提取差异表达分析的结果
      DEGAll <- lrt$table
      #如果使用limma方法
    } else if (method == 'limma') {
      #首先使用voom函数对表达谱数据进行归一化
      v <- voom(dge, design=design, plot = FALSE)
      #对每一个基因进行线性模型拟合
      fit <- lmFit(v, design)
      #计算给定的一组对比的系数和标准误差
      fit2 <- contrasts.fit(fit, contrast.matrix)
      #利用贝叶斯方法计算统计t值，p值
      fit2 <- eBayes(fit2)
      #导出差异表达分析的结果
      DEGAll <- topTable(fit2, coef=1, n = Inf)
      #修改列名
      colnames(DEGAll) <- c('logFC', 'AveExpr', 't', 
                            'PValue', 'FDR', 'B')
    }
  }
  
  #对差异表达分析得到的原始p值进行FDR校正
  DEGAll$FDR <- p.adjust(DEGAll$PValue, method = 'fdr')
  #对FDR进行有小到大排序
  o <- order(DEGAll$FDR)
  DEGAll <- DEGAll[o,]
  
  #将Ensembl基因ID转换成基因名字
  #添加基因的类型，mRNA还是长的非编码RNA，或其他RNA类型
  if (startsWith(rownames(counts)[1], 'ENSG')) {
    degList <- biotype[match(rownames(DEGAll), biotype$ensemblID),]
    degOutput <- data.frame(symbol=degList$geneSymbol, 
                            group=degList$group, DEGAll)
    #如果找不到对应的基因名字，删除这个基因
    keep <- which(! is.na(degOutput$symbol))
    degOutput <- degOutput[keep,]
    return(degOutput)
  } else {
    #如果不是ENSG开头的基因，例如miRNA，则不作任何处理，直接输出
    return (DEGAll)
  }
}

#######
dim(pc_expr)
#[1] 80930    37
DEGAll=DEA(pc_expr,type,comparison="IPF-control")
dim(DEGAll)

library(dplyr)
DEGAll %>% filter(FDR<0.05 & abs(logFC)>1) %>% dim
#[1] 5328    6
deALL=DEGAll %>% filter(FDR<0.05 & abs(logFC)>1)
dim(deALL)
DEG_symbol=id_df[match(rownames(deALL),id_df$X1),6]
sum(DEG_symbol=="-")
deALL$symbol=DEG_symbol
write.csv(file="DEG.csv",deALL,row.names=T)
dim(deALL)
#[1] 5328    7

length(unique(symbol))
#[1] 4083
#######################################
##  miRNA火山图
#######################################

library(ggplot2)
gdcVolcanoPlot2<-function (deg.all, fc = 2, pval = 0.01,size=0.8) {
  geneList <- deg.all
  geneList$threshold <- c()
  geneList$threshold[geneList$logFC > log(fc, 2) & geneList$FDR < 
                       pval] <- 1
  geneList$threshold[geneList$logFC >= -log(fc, 2) & geneList$logFC <= 
                       log(fc, 2) | geneList$FDR >= pval] <- 2
  geneList$threshold[geneList$logFC < -log(fc, 2) & geneList$FDR < 
                       pval] <- 3
  geneList$threshold <- as.factor(geneList$threshold)
  lim <- max(max(geneList$logFC), abs(min(geneList$logFC))) + 
    0.5
  volcano <- ggplot(data = geneList, aes(x = logFC, 
                                         y = -log10(FDR)))
  volcano + geom_point(aes(color = threshold), alpha = 1,
                       size =size) + xlab("log2(Fold Change)") + ylab("-log10(FDR)") + 
    scale_colour_manual( values = c("red", 
                                    "black", "green3")) + xlim(c(-lim, lim)) + geom_vline(xintercept = c(-log(fc, 
                                                                                                              2), log(fc, 2)), color = "darkgreen", linetype = 3) + 
    geom_hline(yintercept = -log(pval, 10), color = "darkgreen", 
               linetype = 3) + theme_bw() + theme(axis.line = element_line(colour = "black"), 
                                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.border = element_rect(colour = "black"), panel.background = element_blank()) + 
    theme(legend.position = "none") + theme(axis.text = element_text(size = 14), 
                                            axis.title = element_text(size = 16))
}

##################
#这里的fc，pval和size可以根据自己的需要进行修改
#size控制火山图中点的大小，默认值是0.8
##################
pdf("mRNA_volcano_default.pdf")
gdcVolcanoPlot2(DEGAll,size=1)
dev.off()

##############################################
#heatmap
##############################################
library(gplots)
deALL=DEGAll %>% filter(FDR<0.05 & abs(logFC)>1)
degName = rownames(deALL)
#rna.expr = DMPsig
degDa <- RNAexpr[degName,]
sampleCol <- ifelse(type == "IPF", "red", "blue")
pdf(file="mRNA_heatmap.pdf",width=12)
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

legend(y=0.46, x=0.92, xpd=TRUE,     
       legend = c("IPF","Normal"),
       col = c("red","blue"), 
       #lty= 1,             
       #lwd = 5,
       pch=15 ,          
       cex=1
       )
dev.off()
###################################################
#DEG PCA
##################################################
library(limma)
pdf(file="mRNA_PCA.pdf",width=7)
sampleCol <- ifelse(type == "IPF", "red", "blue")
pch=ifelse(type=="IPF",15,20)
plotMDS(degDa,labels=NULL,pch=pch,col=sampleCol)
legend("bottomleft",legend=c("IPF","control"),col=c("red","blue"),pch=c(15,20))
dev.off()

############################################
#GO KEGG enrichment
############################################
head(deALL)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
ego <- enrichGO(gene = deALL$symbol,   #要进行GO富集分析的基因的名字
                OrgDb=org.Hs.eg.db,     #注释数据库
                ont = "ALL",            #一次性对BP，MF和CC做富集分析,如果只对其中一种做富集，写对应的BP，MF或CC
                pAdjustMethod = "BH",   #FDR校正p值
                minGSSize = 10,         #富集的GO条目至少要包含10个基因，可以修改
                pvalueCutoff = 0.05,    #p值阈值<0.05，可以修改
                qvalueCutoff = 0.2,    #q值阈值<0.2，可以修改
                keyType='SYMBOL'      #跟第一个参数对应，可以为ENSEMBL，SYMBOL或ENTREZID
) 
write.csv(as.data.frame(ego),"GO_all_enrich.csv",row.names =F)
#######################################################
##柱形图分三个框显示BP,MF和CC的富集分析结果
########################################################
pdf(file="GO_all_bar.pdf",width=10)
barplot(ego, split="ONTOLOGY",showCategory = 10,label_format = 100) + facet_grid(ONTOLOGY~., scale="free")
dev.off()
#######################################################
##气泡图分三个框显示BP,MF和CC的富集分析结果
########################################################
#showCategory = 10可以控制显示几条，默认BP,MF和CC各显示10条
pdf(file="GO_all_dot.pdf",width=10,height=7)
dotplot(ego, split="ONTOLOGY",showCategory = 10,label_format = 100) + facet_grid(ONTOLOGY~., scale="free")
dev.off()

###############################################
#KEGG
###############################################
eg = bitr(deALL$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#进行KEGG富集分析
kk <- enrichKEGG(gene = eg[[2]],    #显著差异表达基因的ENTREZID
                 organism     = 'hsa',   #富集分析的物种，物种缩写可以参考https://www.genome.jp/kegg/catalog/org_list.html
                 pAdjustMethod = "BH",   #FDR校正p值
                 pvalueCutoff = 0.05,    #p值阈值<0.05
                 qvalueCutoff = 0.2,    #q值阈值<0.2
                 minGSSize = 10          #富集的GO条目至少包含10个基因
)
#将富集分析结果中的ENTREZID再转换成基因名字
kk=setReadable(kk,OrgDb=org.Hs.eg.db,keyType  ="ENTREZID")
#############################
##保存KEGG富集结果
#############################
write.csv(as.data.frame(kk),"KEGG-enrich.csv",row.names =F)
###############################################
##气泡图和柱状图显示KEGG富集结果
###############################################
#showCategory = 15可以控制显示几条
pdf(file="KEGG_bar.pdf",width=10)
barplot(kk, showCategory=15,title="KEGG Enrichment")
dev.off()
pdf(file="KEGG_dot.pdf",width=10)
dotplot(kk,showCategory=15,title="KEGG Enrichment")
dev.off()
###############################################
#cibersortx
###############################################
head(deALL)
DEG_expr=data.frame(symbol=deALL$symbol,RNAexpr[rownames(deALL),],stringsAsFactors = F)
head(DEG_expr)
write.table(file="DEG_expr.txt",DEG_expr,quote=F,row.names=F,sep="\t")

DEG_expr_rmdp=aggregate(.~symbol,data=DEG_expr,mean)
write.table(file="DEG_expr_rmdp.txt",DEG_expr_rmdp,quote=F,row.names=F,sep="\t")

#######################

cibersort=read.csv("CIBERSORTx_IPF_Results.csv",header=T,row.names=1,check.names=F)
library(dplyr)
library(tibble)
data=cibersort %>% rownames_to_column("name") %>% 
  #filter(`P-value`<0.05) %>% 
  select(-`P-value`,-Correlation,-RMSE) %>%
  column_to_rownames("name")

#rownames(data)=gsub("\\-\\w{3}$","",rownames(data))
#data=na.omit(data[ESCC_sample,])


###############################
#immune cell composition
##############################

library(tidyr)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
pdf(file="mRNA_cibersort_composition.pdf",width=12)
ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),
        axis.ticks.x = element_blank(),
        legend.position = "right") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))
dev.off()
###########################################
#boxplot group cell composition difference
###########################################
library(ggpubr)
names(type)=rownames(cibersort)

dat$group = type[match(dat$Sample,names(type))]
pdf("cibersort_boxplot_group.pdf",width = 12)
ggplot(dat,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = c("black")) + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1)
        )+
  scale_fill_manual(values = c("blue","red"))+
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "t.test")
dev.off()


########################################
#小提琴图CIBERSORTx vioplot RVA
########################################
#BiocManager::install("vioplot")
library(vioplot) 

#index1=clinical$type %in% c("Control","RVA")
RVA_group=factor(type,levels=c("IPF","control"))
vioplotdata = data
#all(rownames(vioplotdata)==names(RVA_group))

#ciborsort_index=order(RVA_group)
RVA=as.numeric(table(RVA_group)[1])
Control=as.numeric(table(RVA_group)[2])

pdf(paste0("./","CIBERSORTx_allgene_vioplot.pdf"),height=8,width=15)#保存图片的文件名称
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(vioplotdata))
y=c(1:ncol(vioplotdata))
plot(x,y,
     xlim=c(0,(ncol(vioplotdata)-1)*3),#xlim=c(0,63)这个与细胞个数(列数)有关
     #xlim=c(0,57)可以设置成列数-1，再×3
     ylim=c(min(vioplotdata),max(vioplotdata)+0.02),
     
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
#对每个免疫细胞循环，绘制vioplot，正常用green表示，肿瘤用blue表示
for(i in 1:ncol(vioplotdata)){
  tumorData=vioplotdata[1:RVA,i]
  normalData=vioplotdata[(RVA+1):(Control+RVA),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",p)),cex = 0.8)
  text(seq(1,(ncol(vioplotdata)-1)*3+1,3),-0.05,xpd = NA,labels=colnames(vioplotdata),cex = 1,srt = 45,pos=2)
}
legend("topright",legend=c("IPF","control"),col=c("red","blue"),pch=15,bty="n")
dev.off()

######################################
#all 22 cell type heatmap
#####################################
anno = data.frame(group = type, 
                  row.names = rownames(data))
library(pheatmap)
pdf(file="22celltype_heatmap.pdf",width=12)
pheatmap(t(data), #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = T,#列聚类，可以看出样本之间的区分度
         annotation_col =anno, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100)#调色
         #filename = "heatmap_F.pdf",#是否保存
)
dev.off()
##################################################
#sig_cell_type  heatmap  didn't work
#################################################
anno = data.frame(group = type, 
                  row.names = rownames(data))
pval=compare_means(Proportion ~ group,dat,group.by = "Cell_type",method="t.test")
sig_cell_type=pval[pval$p<0.05,"Cell_type"]
library(pheatmap)
pdf(file="22celltype_heatmap.pdf",width=12)
pheatmap(t(data[,unlist(sig_cell_type)]), #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = T,#列聚类，可以看出样本之间的区分度
         annotation_col =anno, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100)#调色
         #filename = "heatmap_F.pdf",#是否保存
)
dev.off()

######################################################
#Estimate
######################################################
#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)

filterCommonGenes(input.f="DEG_expr_rmdp.txt", 
                  output.f="commonGenes.txt", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.txt",
              output.ds="estimateScore.txt")
scores=read.table("estimateScore.txt",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
#rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="ESTIMATE_scores.txt",sep="\t",quote=F,col.names=F)

################################################
#TME estimate score
################################################
library(ggpubr)
TME=read.table("ESTIMATE_scores.txt",header=T,sep="\t",row.names=1)
TME=TME[,-4]
TME$type=type
library(reshape2)
TME_df=melt(TME)
#TME_df$AAGcluster=factor(TME_df$AAGcluster)
pdf(file="TME.pdf",width=10)
ggplot(TME_df, aes(x=variable, y=value, fill=type)) +
  geom_violin(width=0.6) +
  geom_boxplot(width=0.05, color="white", alpha=0.2,position=position_dodge(0.6)) +
  theme_classic()+
  theme(
    legend.position="top",
    text = element_text(size=20),
    plot.title = element_text(size=11)
  ) +
  scale_fill_manual(values=c("blue","red"))+
  #ggtitle("A Violin wrapping a boxplot") +
  xlab(NULL)+
  ylab("TME score")+
  stat_compare_means(aes(group = type),
                     label = "p.signif",
                     method="t.test"
  )
dev.off()
#####################################################
#GSVA 29 immune cell
#####################################################
library(GSVA)
library(clusterProfiler)

immu29<-read.gmt("immune29.gmt")
immu29=unstack(immu29,gene~term)
#load("../1_TCGA_Expr_DEA_survival_code/CHOL1.RData")

rownames(DEG_expr_rmdp)=DEG_expr_rmdp$symbol
DEG_expr_rmdp=DEG_expr_rmdp[,-1]
ssgsea.res <-gsva(
  as.matrix(DEG_expr_rmdp),
  immu29,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = T
)

colnames(ssgsea.res)
################################
#29immune  group boxplot
################################
library(dplyr)
library(tidyr)
library(tibble)
dat <- ssgsea.res  %>% t() %>% as.data.frame() %>% 
  mutate(group = type) %>%
  rownames_to_column("sample") %>% 
  gather(key = Cell_type,value = Proportion,-c(sample, group))
View(dat)

library(ggplot2)
library(ggpubr)
pdf(file="immune29_group_boxplot.pdf", width = 15)
ggplot(dat,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = c("black")) + 
  theme_bw() + 
  labs(x = "Cell Type", y = "GSVA enrichment score") +
  theme(legend.position = "top") + 
  theme(text = element_text(size=12),
    axis.text.x = element_text(angle=45,vjust = 1, hjust = 1))+
  scale_fill_manual(values = c("blue", "red"))+
  labs(x="")+
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "t.test")
dev.off()

################################
#all 26 cell heatmap not good
################################
#identical(colnames(ssgsea.res), rownames(risk))
anno = data.frame(group = type, 
                  row.names = colnames(DEG_expr_rmdp))
library(pheatmap)
pheatmap(ssgsea.res, #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         annotation_col =anno, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100)#调色
         #filename = "heatmap_F.pdf",#是否保存
)
#############################
#19sig cell type heatmap
#############################
pval=compare_means(Proportion ~ group,dat,group.by = "Cell_type",method="t.test")
sig_cell_type=pval[pval$p<0.05,"Cell_type"]
anno = data.frame(group = type, 
                  row.names = colnames(DEG_expr_rmdp))
library(pheatmap)
pheatmap(ssgsea.res[unlist(sig_cell_type),], #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         annotation_col =anno, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100)#调色
         #filename = "heatmap_F.pdf",#是否保存
)
##################################
#
##################################
#library(devtools)
#install_github('ebecht/MCPcounter',ref='master', subdir='Source',force = TRUE)
library(MCPcounter)
genes=read.table("Signatures/genes.txt",sep='\t',stringsAsFactors=FALSE,header=TRUE,colClasses='character',check.names=FALSE)
results<- MCPcounter.estimate(DEG_expr_rmdp,
                              featuresType= "HUGO_symbols",
                              genes=genes
)
results2 <- results[1:9,]#去除最后一行fbro细胞的行
#colnames(results2) = gsub("\\.", "-", colnames(results2))
colnames(results2)
#identical(colnames(results2), rownames(risk))
anno = data.frame(group = type, 
                  row.names = colnames(DEG_expr_rmdp))
pdf(file="MCPcounter_heatmap.pdf",height=5)
library(pheatmap)
pheatmap(results2, #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = T,#列聚类，可以看出样本之间的区分度
         annotation_col =anno, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100)#调色
         #filename = "heatmap_F.pdf",#是否保存
)
dev.off()


library(dplyr)
library(tidyr)
library(tibble)
dat <- results2  %>% t() %>% as.data.frame() %>% 
  mutate(group = type) %>%
  rownames_to_column("sample") %>% 
  gather(key = Cell_type,value = Proportion,-c(sample, group))
View(dat)

library(ggplot2)
library(ggpubr)
pdf("mcpcounter_boxplot.pdf", width = 12)
ggplot(dat,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = c("black")) + 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,vjust = 1, hjust = 1))+
  scale_fill_manual(values = c("blue", "red"), 
                    )+
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test")
dev.off()

###############################################
#xCell failed
###############################################
#remotes::install_github("dviraran/xCell")
library(xCell)
xcell.result = xCellAnalysis(DEG_expr_rmdp,rnaseq = T)
# [1] "Num. of genes: 2570"
# [1] "ERROR: not enough genes"
# Error in return - 1 : non-numeric argument to binary operator


##############################
# epic
###############################
#library(devtools)
#'testit', 'limSolve', 'ComICS'
#remotes::install_github("icbi-lab/immunedeconv")
#安装之后有一些编码的warning，
library(immunedeconv)
result=immunedeconv::deconvolute(DEG_expr_rmdp, "epic")

results2 <- data.frame(result,check.names = F)#去除最后一行fbro细胞的行
rownames(results2)=result$cell_type
results2=results2[,-1]
#colnames(results2) = gsub("\\.", "-", colnames(results2))
colnames(results2)
#identical(colnames(results2), rownames(risk))
anno = data.frame(group = type, 
                  row.names = names(type))
pdf(file="EPIC_heatmap.pdf",height=4)
library(pheatmap)
pheatmap(results2, #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = T,#列聚类，可以看出样本之间的区分度
         annotation_col =anno, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100)#调色
         #filename = "heatmap_F.pdf",#是否保存
)
dev.off()

library(dplyr)
library(tidyr)
library(tibble)
dat <- results2  %>% t() %>% as.data.frame() %>% 
  mutate(group = type) %>%
  rownames_to_column("sample") %>% 
  gather(key = Cell_type,value = Proportion,-c(sample, group))
View(dat)

library(ggplot2)
library(ggpubr)
pdf(file="EPIC_boxplot.pdf", width = 12)
ggplot(dat,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = c("black")) + 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top",
        #axis.text = element_text(size = 14),
        text = element_text(size=12)
        ) + 
  theme(axis.text.x = element_text(angle=45,vjust = 1, hjust = 1))+
  scale_fill_manual(values = c("red", "blue"), 
                    )+
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "t.test")
dev.off()


#######################################
#TIMER
#######################################
results2=deconvolute_timer(DEG_expr_rmdp, indications = rep('luad',ncol(DEG_expr_rmdp)))

#colnames(results2) = gsub("\\.", "-", colnames(results2))
colnames(results2)
#identical(colnames(results2), rownames(risk))
anno = data.frame(group = type, 
                  row.names = names(type))
pdf(file="timer_heatmap.pdf",height=4)
library(pheatmap)
pheatmap(results2, #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = T,#列聚类，可以看出样本之间的区分度
         annotation_col =anno, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100)#调色
         #filename = "heatmap_F.pdf",#是否保存
)
dev.off()


library(dplyr)
library(tidyr)
library(tibble)
dat <- results2  %>% t() %>% as.data.frame() %>% 
  mutate(group = risk$risk) %>%
  rownames_to_column("sample") %>% 
  gather(key = Cell_type,value = Proportion,-c(sample, group))
View(dat)

library(ggplot2)
library(ggpubr)
pdf(file="TIMER_boxplot.pdf", width = 10)
ggplot(dat,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = c("black")) + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top",
        text = element_text(size=12)
        ) + 
  theme(axis.text.x = element_text(angle=45,vjust = 1, hjust = 1))+
  scale_fill_manual(values = c("blue", "red") 
                    )+
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "t.test")
dev.off()
##################################################
#GSVA KEGG
##################################################
##########################################
#GSVA 4083 DEG   (used)
##########################################

library(GSVA)
library(fgsea)
library(clusterProfiler)
#选取非零基因


#转换成矩阵
expr <- as.matrix(DEG_expr_rmdp)


#从GSEA官网下载GSEA分析需要的基因集
#http://www.gsea-msigdb.org/gsea/index.jsp

#下载免疫相关的基因集，h: hallmark gene sets
gmtfile ='genesets/c2.cp.kegg.v7.5.1.symbols.gmt'

#读取gmt文件中的pathway信息
pathway<-read.gmt(gmtfile)[,c(2,1)]
#去堆叠,转换成list
genesets=unstack(pathway)

#进行GSVA分析
gsva.res <- gsva(expr, genesets, method="ssgsea") 
#保存GSVA分析结果
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res2.csv", row.names = F)

###################################
#DEG   GSVA pathways
####################################
library(limma)
#sample_type=factor(cluster[index],levels=c(1,2),labels=c("A","B"))
#构建分组模型
design <- model.matrix(~ type + 0)
#修改列名
colnames(design) <- levels(type)
#进行线性拟合
fit <- lmFit(gsva.res, design)
#创建比较模型
cont.matrix <- makeContrasts(IPF-control, levels=design)
#进行比较，寻找差异表达基因
fit2 <- contrasts.fit(fit, cont.matrix)
#使用经验贝叶斯方法来调整估计的倍数变化的标准误差
fit2 <- eBayes(fit2)
#输出差异表达分析结果
DEG <- topTable(fit2,coef=1, number=Inf)
#通过探针名字获取相应的基因名字


#筛选显著差异表达的基因
index1=(DEG$adj.P.Val<0.05 & abs(DEG$logFC)>0.1) 
sum(index1)
DEG_sig=DEG[index1,]
#保存所有结果
View(DEG_sig)



#绘制热图
library(pheatmap)
pdf(file="KEGG_GSVA_heatmap.pdf",width=14,height=8)
annotation_col = data.frame(type=factor(type))
pheatmap(gsva.res[rownames(DEG_sig),], 
         show_colnames = F, scale = "row",
         annotation_col = annotation_col,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols = T,
         #treeheight_row = 0,
         cluster_rows = T
         
)
dev.off()
