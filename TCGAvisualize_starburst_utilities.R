getMetPlatInfo <- function(genome, platform) {
  base <- "http://zwdzwd.io/InfiniumAnnotation/current"
  
  platform <-  switch(platform,
                      "450K" = "hm450",
                      "450k" = "hm450",
                      "27k"  = "hm27",
                      "27K"  = "hm27",
                      "EPIC" = "EPIC")
  path <- file.path(base,platform,paste(platform,"hg19.manifest.rds", sep ="."))
  if (grepl("hg38", genome)) path <- gsub("hg19","hg38",path)
  message(path)
  if(Sys.info()["sysname"] == "Windows") mode <- "wb" else  mode <- "w"
  if(!file.exists(basename(path))) downloader::download(path,basename(path), mode = mode)
  gr <- readRDS(basename(path))
  return(gr)
}

getTSS <- function(genome = "hg38", TSS = list(upstream = NULL, downstream = NULL)){
  
  host <- ifelse(genome == "hg19",  "grch37.ensembl.org","www.ensembl.org")
  ensembl <- tryCatch({
    useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host =  host)
  },  error = function(e) {
    useEnsembl("ensembl",
               dataset = "hsapiens_gene_ensembl",
               mirror = "uswest",
               host =  host)
  })
  attributes <- c("chromosome_name",
                  "start_position",
                  "end_position",
                  "strand",
                  "transcript_start",
                  "transcription_start_site",
                  "transcript_end",
                  "ensembl_transcript_id",
                  "ensembl_gene_id",
                  "entrezgene",
                  "external_gene_name")
  
  chrom <- c(1:22, "X", "Y","M","*")
  db.datasets <- listDatasets(ensembl)
  description <- db.datasets[db.datasets$dataset=="hsapiens_gene_ensembl",]$description
  message(paste0("Downloading transcripts information. Using: ", description))
  
  filename <-  paste0(gsub("[[:punct:]]| ", "_",description),"_tss.rda")
  if(!file.exists(filename)) {
    tss <- getBM(attributes = attributes, filters = c("chromosome_name"), values = list(chrom), mart = ensembl)
    tss <- tss[!duplicated(tss$ensembl_transcript_id),]
    save(tss, file = filename)
  } else {
    tss <- get(load(filename))
  }
  if(genome == "hg19") tss$external_gene_name <- tss$external_gene_id
  tss$chromosome_name <-  paste0("chr", tss$chromosome_name)
  tss$strand[tss$strand == 1] <- "+"
  tss$strand[tss$strand == -1] <- "-"
  tss <- makeGRangesFromDataFrame(tss,start.field = "transcript_start", end.field = "transcript_end", keep.extra.columns = TRUE)
  
  if(!is.null(TSS$upstream) & !is.null(TSS$downstream))
    tss <- promoters(tss, upstream = TSS$upstream, downstream = TSS$downstream)
  
  return(tss)
}

map.ensg <- function(genome = "hg38", genes) {
  gene.location <- get.GRCh.bioMart(genome)
  gene.location <- gene.location[match(genes,gene.location$ensembl_gene_id),]
  return(gene.location)
}


bquote(atop(Qua[0.99] == phantom(), .(round(3.1415927,4))))
plot(1:3,main=bquote(atop(Qua[0.99] == phantom(), .(round(3.1415927,4)))))
plot(1:3,main=bquote(atop("Starburst Plot", scriptstyle((list(Delta, 
                                                              bar(beta) >= .(diffmean.cut), group("|", 
                                                                                                  logFC, "|") >= .(logFC.cut), FDR[expression] <= 
                                                                .(exp.p.cut), FDR[DNAmethylation] <= .(met.p.cut)))))))
