#TCGAbiolinks_2.14.0
library(TCGAbiolinks)
library(SummarizedExperiment)
library(minfi)
#BiocManager::install("ELMER.data")
library(ELMER.data)

library(edgeR)
ployA=read.table("GSE173355_PolyA_Counts.txt",header=T,sep=" ",check.names = F,stringsAsFactors = F)
dim(ployA)
ployA[1:3,1:3]

ids=strsplit(ployA$target_id,"\\|")
id_df=data.frame(do.call(rbind,ids))
sum(id_df$X8=="protein_coding")
#[1] 80930
pc_index=id_df$X8=="protein_coding"
pc_expr=ployA[pc_index,]

index=grepl("ENSG\\d+\\.\\d+$",id_df$X2[pc_index])
eng_id=gsub("\\.\\d+","",id_df$X2[pc_index][index])

pc_expr1=pc_expr[index,]
pc_expr1$target_id=eng_id

pc_expr2=aggregate(.~target_id,pc_expr1,mean)
rownames(pc_expr2)=pc_expr2$target_id
pc_expr2=pc_expr2[,-1]
logexpr=log2(pc_expr2+1)
boxplot(logexpr)
#logexpr_norm=sweep(logexpr,2, apply(logexpr,2,median,na.rm=T))
#boxplot(logexpr_norm)


load("TCGAbiolinks.RData")
#load("sysdata.rda")
#load("PolyA_Counts2.RData")
colData <- DataFrame(type=targets$type,
                     row.names=targets$name)
colnames(bVals)=targets$name
IPF.met <- SummarizedExperiment(assays=SimpleList(counts=bVals[,targets$name]),
                                 rowRanges = myAnnotation@ranges,
                                 colData=colData)
colData(IPF.met)
#这里不分case和control，case.vs.control和control.vs.case的p值都会算出来
#写法需要注意diffmean.control.IPF,是IPF-control的差值
IPF.met.dmc <- TCGAanalyze_DMR(IPF.met, groupCol = "type",
                                p.cut = 0.05, #10^-5
                                diffmean.cut = 0.1,
                                #legend = "State",
                                plot.filename = "tumorvsnormal_metvolcano.pdf")


#names(DEG_expr_rmdp)=gsub("\\.","-",names(DEG_expr_rmdp))
#ens=gsub("\\.\\d+$","",id_df[match(rownames(DEG_expr_rmdp),id_df$X6),2])
#rownames(DEG_expr_rmdp)=ens



#colnames(rgset)
inter=intersect(names(logexpr),colnames(bVals))
length(inter)
#36
dataFilt=logexpr[,inter]
idx=grepl("^LTRC",names(dataFilt))
idx2=!idx
#limma cond2-cond1
dataDEGs <- TCGAanalyze_DEA(mat2 = as.matrix(dataFilt[,idx]),
                            mat1 = as.matrix(dataFilt[,idx2]),
                            metadata = F,
                            pipeline = "limma",
                            Cond1type = "control",
                            Cond2type = "IPF",
                            #method = "exactTest",
                            #fdr.cut = 0.05,
                            #logFC.cut = 1
                            )


TCGAVisualize_volcano(x = dataDEGs$logFC,
                      y = dataDEGs$adj.P.Val,
                      filename = "IPF_volcanoexp.pdf",
                      x.cut = 1,
                      y.cut = 0.05,
                      #names = rownames(dataDEGs),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot (IPF vs control)",
                      width = 10)

#test=load("EPIC.hg19.manifest.rda")
#saveRDS(EPIC.hg19.manifest,file="EPIC.hg19.manifest.rds")

names(dataDEGs)[5]="FDR"
dataDEGs$ensembl_gene_id=rownames(dataDEGs)
starburst <- TCGAvisualize_starburst(met = IPF.met.dmc, 
                                     exp = dataDEGs,
                                     genome = "hg19",
                                     group1 = "control",
                                     group2 = "IPF",
                                     filename = "starburst.pdf",
                                     met.platform = "EPIC",
                                     met.p.cut = 0.05,
                                     exp.p.cut = 0.05,
                                     diffmean.cut = 0.1,
                                     logFC.cut = 1,
                                     names = F, 
                                     height=10,
                                     width=20,
                                     dpi=300)

save(IPF.met,inter,targets,dataFilt,IPF.met.dmc,logexpr,dataDEGs,file="TCGAbiolinks.RData")
