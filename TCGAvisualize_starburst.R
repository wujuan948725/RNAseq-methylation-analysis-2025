
function (met, exp, group1 = NULL, group2 = NULL, exp.p.cut = 0.01, 
          met.p.cut = 0.01, diffmean.cut = 0, logFC.cut = 0, met.platform, 
          genome, names = FALSE, names.fill = TRUE, filename = "starburst.pdf", 
          return.plot = FALSE, ylab = expression(atop("Gene Expression", 
                                                      paste(Log[10], " (FDR corrected P values)"))), 
          xlab = expression(atop("DNA Methylation", paste(Log[10], 
                                                          " (FDR corrected P values)"))), title = "Starburst Plot", 
          legend = "DNA Methylation/Expression Relation", color = NULL, 
          label = c("Not Significant", "Up regulated & Hypo methylated", 
                    "Down regulated & Hypo methylated", "hypo methylated", 
                    "hyper methylated", "Up regulated", "Down regulated", 
                    "Up regulated & Hyper methylated", "Down regulated & Hyper methylated"), 
          xlim = NULL, ylim = NULL, height = 10, width = 20, dpi = 600) 
{
  
  
  library(biomaRt)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  color = NULL
  names.fill = TRUE
  title = "Starburst Plot"
  label = c("Not Significant", "Up regulated & Hypo methylated", 
            "Down regulated & Hypo methylated", "hypo methylated", 
            "hyper methylated", "Up regulated", "Down regulated", 
            "Up regulated & Hyper methylated", "Down regulated & Hyper methylated")
  met = IPF.met.dmc
  exp = dataDEGs
  genome = "hg19"
  group1 = "control"
  #group2会显示在图上,label
  group2 = "IPF"
  filename = "starburst.pdf"
  met.platform = "EPIC"
  met.p.cut = 0.05
  exp.p.cut = 0.05
  diffmean.cut = 0.1
  logFC.cut = 1
  names = F
  height=10
  width=20
  dpi=300
  title = "Starburst Plot"
  xlim = NULL
  ylim = NULL
  return.plot = FALSE
  ylab = expression(atop("Gene Expression", 
                         paste(Log[10], " (FDR corrected P values)")))
  xlab = expression(atop("DNA Methylation", paste(Log[10], 
                                                  " (FDR corrected P values)")))
  legend = "DNA Methylation/Expression Relation"
  
  
  
  
  if (missing(genome)) 
    stop("Please set genome (hg19 or hg38)")
  if (missing(met.platform)) 
    stop("Please set met.platform (EPIC, 450K or 27K)")
  .e <- environment()
  group1.col <- gsub("[[:punct:]]| ", ".", group1)
  group2.col <- gsub("[[:punct:]]| ", ".", group2)
  if (title == "Starburst Plot") {
    if (diffmean.cut != 0 & logFC.cut == 0) {
      title <- bquote(atop("Starburst Plot", scriptstyle((list(Delta ~ 
                                                                 bar(beta) >= .(diffmean.cut), FDR[expression] <= 
                                                                 .(exp.p.cut), FDR[DNAmethylation] <= .(met.p.cut))))))
    }
    else if (logFC.cut != 0 & diffmean.cut == 0) {
      title <- bquote(atop("Starburst Plot", scriptstyle((list(group("|", 
                                                                     logFC, "|") >= .(logFC.cut), FDR[expression] <= 
                                                                 .(exp.p.cut), FDR[DNAmethylation] <= .(met.p.cut))))))
    }
    else if (logFC.cut != 0 & diffmean.cut != 0) {
      title <- bquote(atop("Starburst Plot", scriptstyle((list(Delta, 
                                                               bar(beta) >= .(diffmean.cut), group("|", 
                                                                                                   logFC, "|") >= .(logFC.cut), FDR[expression] <= 
                                                                 .(exp.p.cut), FDR[DNAmethylation] <= .(met.p.cut))))))
    }
  }
  if (is.null(color)) 
    color <- c("#000000", "#E69F00", "#56B4E9", 
               "#009E73", "red", "#0072B2", "#D55E00", 
               "#CC79A7", "purple")
  names(color) <- as.character(1:9)
  names(label) <- as.character(1:9)
  names.color <- color
  names(names.color) <- label
  if (is.null(group1) || is.null(group2)) {
    message("Please, set the group1 and group2 parameters")
    return(NULL)
  }
  if (class(met) == class(as(SummarizedExperiment(), "RangedSummarizedExperiment"))) {
    met <- values(met)
  }
  pcol <- gsub("[[:punct:]]| ", ".", paste("p.value.adj", 
                                           group1, group2, sep = "."))
  if (!(pcol %in% colnames(met))) {
    pcol <- gsub("[[:punct:]]| ", ".", paste("p.value.adj", 
                                             group2, group1, sep = "."))
  }
  if (!(pcol %in% colnames(met))) {
    stop("Error! p-values adjusted not found. Please, run TCGAanalyze_DMR")
  }
  fctr.cols <- sapply(exp, is.factor)
  exp[, fctr.cols] <- sapply(exp[, fctr.cols], as.character)
  idx <- grep("ENSG", exp[1, ])
  if (length(idx) > 0) {
    colnames(exp)[idx] <- "ensembl_gene_id"
  }else{
    if (grepl("ENSG", rownames(exp)[1])) {
      exp$ensembl_gene_id <- rownames(exp)
    }
    else {
      gene.info <- get.GRCh.bioMart(genome)
      if (any(sapply(rownames(exp)[1:10], function(y) any(grepl(y, 
                                                                gene.info[, grep("external_gene_", colnames(gene.info))]))))) {
        idx <- match(rownames(exp), gene.info[, grep("external_gene_", 
                                                     colnames(gene.info))])
        exp <- exp[!is.na(idx), ]
        idx <- match(rownames(exp), gene.info[, grep("external_gene_", 
                                                     colnames(gene.info))])
        exp$ensembl_gene_id <- gene.info[idx, "ensembl_gene_id"]
      }
    }
  }
  if (!"probeID" %in% colnames(met)) 
    met$probeID <- rownames(IPF.met.dmc)
  if (!"ensembl_gene_id" %in% colnames(exp)) {
    stop("Column ensembl_gene_id was not found")
  }
  message("o Fetching auxiliary information")
  message("oo Fetching probes genomic information")
  met.info <- getMetPlatInfo(genome = genome, platform = met.platform)
  values(met.info) <- NULL
  message("oo Fetching TSS information")
  tss <- getTSS(genome = genome)
  tss <- promoters(tss, upstream = 0, downstream = 0)
  message("o Mapping probes to nearest TSS")
  dist <- distanceToNearest(tss, met.info)
  g <- suppressWarnings(as.data.frame(tss[queryHits(dist)]))
  g$start_position <- NULL
  g$end_position <- NULL
  colnames(g)[1:5] <- paste0("gene_", colnames(g)[1:5])
  m <- suppressWarnings(as.data.frame(met.info[subjectHits(dist)], 
                                      row.names = NULL))
  colnames(m) <- paste0("probe_", colnames(m))
  m$probeID <- names(met.info[subjectHits(dist)])
  nearest <- cbind(m, g)
  nearest$distance_TSS <- values(dist)$distance
  nearest$id <- paste0(nearest$probeID, nearest$ensembl_gene_id)
  nearest <- nearest[order(nearest$id, -abs(nearest$distance_TSS)), 
  ]
  nearest <- nearest[!duplicated(nearest$id), ]
  nearest$id <- NULL
  nearest <- nearest[, c("distance_TSS", "probeID", 
                         "ensembl_gene_id")]
  message("o Mapping results information")
  volcano <- plyr::join(nearest, exp, by = "ensembl_gene_id")
  volcano <- merge(volcano, met, by = "probeID")
  #99454    
  volcano$ID <- paste(volcano$ensembl_gene_id, volcano$probeID, 
                      sep = ".")
  volcano <- volcano[!is.na(volcano$FDR), ]
  volcano$geFDR <- log10(volcano$FDR)
  volcano$geFDR2 <- volcano$geFDR
  volcano[volcano$logFC > 0, "geFDR2"] <- -1 * volcano[volcano$logFC > 
                                                         0, "geFDR"]
  #这里写法需要注意diffmean.control.IPF,是IPF-control的差值
  #group2-group1的差值
  diffcol <- gsub("[[:punct:]]| ", ".", paste("diffmean", 
                                              group1, group2, sep = "."))
  volcano$meFDR <- log10(volcano[, pcol])
  volcano$meFDR2 <- volcano$meFDR
  idx <- volcano[, diffcol] > 0
  idx[is.na(idx)] <- FALSE
  volcano[idx, "meFDR2"] <- -1 * volcano[idx, "meFDR"]
  label[2:9] <- paste(label[2:9], "in", group2)
  met.lowerthr <- log10(met.p.cut)
  met.upperthr <- (-met.lowerthr)
  exp.lowerthr <- log10(exp.p.cut)
  exp.upperthr <- (-exp.lowerthr)
  volcano <- suppressWarnings(tibble::as.tibble(volcano))
  # Group 1: not sifnificant
  # Group 2:up regulated and hypomethylated
  a <- dplyr::filter(volcano, volcano$geFDR2 > exp.upperthr & 
                       volcano$meFDR2 < met.lowerthr & volcano[, diffcol] < 
                       (-diffmean.cut) & volcano$logFC > logFC.cut)
  # Group 3: down regulated and hypomethylated
  b <- dplyr::filter(volcano, volcano$geFDR2 < exp.lowerthr & 
                       volcano$meFDR2 < met.lowerthr & volcano[, diffcol] < 
                       (-diffmean.cut) & volcano$logFC < (-logFC.cut))
  # Group 4: hypomethylated
  c <- dplyr::filter(volcano, volcano$geFDR2 > exp.lowerthr & 
                       volcano$geFDR2 < exp.upperthr & volcano$meFDR2 < met.lowerthr & 
                       volcano[, diffcol] < (-diffmean.cut))
  # Group 5: hypermethylated
  d <- dplyr::filter(volcano, volcano$geFDR2 > exp.lowerthr & 
                       volcano$geFDR2 < exp.upperthr & volcano$meFDR2 > met.upperthr & 
                       volcano[, diffcol] > diffmean.cut)
  # Group 6: upregulated
  e <- dplyr::filter(volcano, volcano$geFDR2 > exp.upperthr & 
                       volcano$meFDR2 < met.upperthr & volcano$meFDR2 > met.lowerthr & 
                       volcano$logFC > logFC.cut)
  # Group 7: downregulated
  f <- dplyr::filter(volcano, volcano$geFDR2 < exp.lowerthr & 
                       volcano$meFDR2 < met.upperthr & volcano$meFDR2 > met.lowerthr & 
                       volcano$logFC < (-logFC.cut))
  # Group 8: upregulated and hypermethylated
  g <- dplyr::filter(volcano, volcano$geFDR2 > exp.upperthr & 
                       volcano$meFDR2 > met.upperthr & volcano[, diffcol] > 
                       diffmean.cut & volcano$logFC > logFC.cut)
  # Group 9: downregulated and hypermethylated
  h <- dplyr::filter(volcano, volcano$geFDR2 < exp.lowerthr & 
                       volcano$meFDR2 > met.upperthr & volcano[, diffcol] > 
                       diffmean.cut & volcano$logFC < (-logFC.cut))
  groups <- as.character(seq(2, 9))
  volcano$starburst.status <- "Not Significant"
  volcano$shape <- "1"
  volcano$threshold.starburst <- "1"
  volcano$threshold.size <- "1"
  state <- c("Up regulated & Hypo methylated", "Down regulated & Hypo methylated", 
             "hypo methylated", "hyper methylated", "Up regulated", 
             "Down regulated", "Up regulated & Hyper methylated", 
             "Down regulated & Hyper methylated")
  s <- list(a, b, c, d, e, f, g, h)
  for (i in seq_along(s)) {
    idx <- s[[i]]$ID
    if (length(idx) > 0) {
      volcano[volcano$ID %in% idx, "threshold.starburst"] <- groups[i]
      volcano[volcano$ID %in% idx, "starburst.status"] <- state[i]
    }
  }
  s <- list(a, b, g, h)
  significant <- NULL
  for (i in seq_along(s)) {
    idx <- s[[i]]$ID
    if (length(idx) > 0) {
      significant <- rbind(significant, volcano[volcano$ID %in% 
                                                  idx, ])
    }
  }
  message("o Plotting figure")
  volcano.aux <- volcano
  p <- ggplot(data = volcano.aux, environment = .e, aes(x = volcano.aux$meFDR2, 
                                                        y = volcano.aux$geFDR2, colour = volcano.aux$threshold.starburst)) + 
    geom_point()
  if (names == TRUE & !is.null(significant)) {
    message("oo Adding names to genes")
    if (names.fill) {
      p <- p + geom_label_repel(data = significant, aes(x = significant$meFDR2, 
                                                        y = significant$geFDR2, label = map.ensg(genome, 
                                                                                                 significant$ensembl_gene_id)$external_gene_name, 
                                                        fill = as.factor(significant$starburst.status)), 
                                size = 4, show.legend = FALSE, fontface = "bold", 
                                color = "black", box.padding = unit(0.35, 
                                                                    "lines"), point.padding = unit(0.3, "lines")) + 
        scale_fill_manual(values = names.color)
    }else {
      p <- p + geom_text_repel(data = significant, aes(x = significant$meFDR2, 
                                                       y = significant$geFDR2, label = map.ensg(genome, 
                                                                                                significant$ensembl_gene_id)$external_gene_name, 
                                                       fill = significant$starburst.status), size = 4, 
                               show.legend = FALSE, fontface = "bold", 
                               color = "black", point.padding = unit(0.3, 
                                                                     "lines"), box.padding = unit(0.5, "lines"))
    }
  }
  if (!is.null(xlim)) {
    p <- p + xlim(xlim)
  }
  if (!is.null(ylim)) {
    p <- p + ylim(ylim)
  }
  title=bquote(atop("Starburst Plot", scriptstyle((list(Delta, 
                                                        bar(beta) >= .(diffmean.cut), group("|", 
                                                                                            logFC, "|") >= .(logFC.cut), FDR[expression] <= 
                                                          .(exp.p.cut), FDR[DNAmethylation] <= .(met.p.cut))))))
  p <- p + ggtitle(title) + ylab(ylab) + xlab(xlab) + guides(size = FALSE)
  #p <- p + ggtitle("title") + ylab(ylab) + xlab(xlab) + guides(size = FALSE)
  
  p <- p + scale_color_manual(values = color, labels = label, 
                              name = legend) + guides(col = guide_legend(nrow = 3))
  p <- p + geom_hline(aes(yintercept = exp.lowerthr), colour = "black", 
                      linetype = "dashed") + geom_hline(aes(yintercept = exp.upperthr), 
                                                        colour = "black", linetype = "dashed") + 
    geom_vline(aes(xintercept = met.lowerthr), colour = "black", 
               linetype = "dashed") + geom_vline(aes(xintercept = met.upperthr), 
                                                 colour = "black", linetype = "dashed") + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), 
                       axis.line.y = element_line(colour = "black"), legend.position = "top", 
                       legend.key = element_rect(colour = "white"), plot.title = element_text(face = "bold", 
                                                                                              size = 16, hjust = 0.5), legend.text = element_text(size = 14), 
                       legend.title = element_text(size = 14), axis.text = element_text(size = 14), 
                       axis.title.x = element_text(face = "bold", size = 14), 
                       axis.text.x = element_text(vjust = 0.5, size = 14), axis.title.y = element_text(face = "bold", 
                                                                                                       size = 14), axis.text.y = element_text(size = 14))
  if (!return.plot) 
    ggsave(filename = filename, width = width, height = height, 
           dpi = dpi)
  volcano$shape <- NULL
  volcano$threshold.starburst <- NULL
  volcano$threshold.size <- NULL
  volcano <- dplyr::filter(volcano, volcano$geFDR <= exp.lowerthr & 
                             volcano$meFDR <= met.lowerthr)
  if (diffmean.cut != 0) {
    volcano <- dplyr::filter(volcano, abs(volcano[, diffcol]) > 
                               diffmean.cut)
  }
  if (logFC.cut != 0) {
    volcano <- dplyr::filter(volcano, abs(volcano$logFC) >= 
                               logFC.cut)
  }
  message("o Saving results")
  message("oo Saving significant results as: starburst_results.csv")
  message("oo It contains pair with changes both in the expression level of the nearest gene and  in the DNA methylation level")
  suppressWarnings(write_csv(x = volcano, path = "starburst_results.csv"))
  if (return.plot) {
    return(list(plot = p, starburst = volcano))
  }
  return(volcano)
}