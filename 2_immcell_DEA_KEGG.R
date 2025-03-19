setwd("C:/other/project/Dr.Kong/IPF_meth")
load("PolyA_Counts2.RData")
cibersort
dim(cibersort)
index=grepl("LTRC",rownames(cibersort))

IPF_cibersort=cibersort[index,]
cells=c("Plasma cells","Macrophages M0","Macrophages M1","Macrophages M2","Mast cells resting","Mast cells activated")

RNAexpr[1:3,1:3]


library(DOSE)
library(clusterProfiler)
library(limma)
IPF_expr=RNAexpr[,index]
dim(IPF_expr)
dir.create("immcell_DEA_KEGG")
all_sig=list()
for(cell in cells){
  #cell=cells[3]
  percent=IPF_cibersort[,cell]
  group=factor(ifelse(percent>median(percent),"high","low"))


  #构建分组模型
  design <- model.matrix(~ group + 0)
  #修改列名
  colnames(design) <- levels(group)
  #进行线性拟合
  fit <- lmFit(IPF_expr, design)
  #创建比较模型
  cont.matrix <- makeContrasts(high-low, levels=design)
  #进行比较，寻找差异表达基因
  fit2 <- contrasts.fit(fit, cont.matrix)
  #使用经验贝叶斯方法来调整估计的倍数变化的标准误差
  fit2 <- eBayes(fit2)
  #输出差异表达分析结果
  DEG <- topTable(fit2,coef=1, number=Inf)
  
  DEG_symbol=id_df[match(rownames(DEG),id_df$X1),6]
  DEG$symbol=DEG_symbol
  #筛选显著差异表达的基因
  index=(DEG$P.Value<0.05 & abs(DEG$logFC)>1) 
  sum(index)
  #保存所有结果
  #result=data.frame(DEG,symbol)
  write.table(file=paste0("immcell_DEA_KEGG/",cell,"_all_DEG.txt"),DEG,quote=F,sep="\t")
  
  #保存显著差异表达基因的结果
  sig_result=DEG[index,]
  #DEG_sig_symbol=id_df[match(rownames(sig_result),id_df$X1),6]
  #sig_result$symbol=DEG_sig_symbol
  dim(sig_result)
  write.table(file=paste0("immcell_DEA_KEGG/",cell,"_sig_DEG.txt"),sig_result,quote=F,sep="\t")
  all_sig[[cell]]=rownames(sig_result)
  
  eg = bitr(sig_result$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  gene_list=eg[[2]]
  kk <- enrichKEGG(gene = gene_list,    #这里需要输入一串gene，ID只能为ENTREZID
                   organism     = 'hsa',   #富集分析的物种，物种缩写可以参考https://www.genome.jp/kegg/catalog/org_list.html
                   pAdjustMethod = "BH",   #FDR校正p值
                   pvalueCutoff = 0.05,    #p值阈值<0.05
                   qvalueCutoff = 0.2,    #q值阈值<0.2
                   minGSSize = 10          #富集的GO条目至少包含10个基因
  )
  #将富集分析结果中的ENTREZID再转换成基因名字
  kk=setReadable(kk,OrgDb=org.Hs.eg.db,keyType  ="ENTREZID")
  if(nrow(kk)>0){
  write.csv(as.data.frame(kk),file=paste0("immcell_DEA_KEGG/",cell,"_KEGG.csv"),row.names =F)
  
  ###############################################
  ##气泡图和柱状图显示KEGG富集结果
  ###############################################
  
  pdf(file=paste0("immcell_DEA_KEGG/",cell,"_KEGG_dot.pdf"),width=10)
  print(dotplot(kk,showCategory=15,title="KEGG Enrichment"))
  dev.off()
  }
}

Reduce(intersect,all_sig)
sapply(all_sig,length)


