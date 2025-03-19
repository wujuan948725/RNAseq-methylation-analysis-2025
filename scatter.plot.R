function (data, byPair = list(probe = c(), gene = c()), byProbe = list(probe = c(), 
                                                                       numFlankingGenes = 20), byTF = list(TF = c(), probe = c()), 
          category = NULL, ylim = NULL, dots.size = 0.9, correlation = FALSE, 
          width = 7, height = 6, dir.out = "./", save = TRUE, 
          ...) 
{
  dir.create(dir.out, recursive = TRUE, showWarnings = FALSE)
  simpleCap <- function(x) {
    if (is.na(x)) 
      return("NA")
    s <- x
    paste(toupper(substring(s, first = 1, last = 1)), tolower(substring(s, 
                                                                        2)), sep = "", collapse = " ")
  }
  if (missing(data)) 
    stop("A data object should be included.")
  if (!is.null(category) && length(category) == 1) {
    if (!category %in% colnames(colData(data))) {
      stop("category not found in the  phenotypic data (colData(data)) ")
    }
    if (is.null(category)) 
      stop("Please, set category argument")
    legend.title <- simpleCap(category)
    samples <- sampleMap(data)[sampleMap(data)$assay == "DNA methylation", 
                               "primary"]
    category <- colData(data)[samples, category]
    if (!"color.value" %in% names(list(...))) 
      category <- sapply(category, simpleCap)
  }
  if (length(byPair$probe) != 0) {
    if (length(byPair$probe) != length(byPair$gene)) 
      stop("In pairs, the length of probes should be the same with the length of genes.")
    pb <- txtProgressBar(min = 0, max = length(byPair$gene), 
                         title = "creating images", style = 3, initial = 0, 
                         char = "=")
    for (i in 1:length(byPair$probe)) {
      setTxtProgressBar(pb, i)
      probe <- byPair$probe[i]
      gene <- byPair$gene[i]
      symbol <- getSymbol(data, geneID = gene)
      P <- scatter(meth = assay(getMet(data)[probe, ]), 
                   exp = assay(getExp(data)[gene, ]), category = category, 
                   ylim = ylim, dots.size = dots.size, legend.title = legend.title, 
                   correlation = correlation, xlab = sprintf("DNA methylation at %s", 
                                                             probe), ylab = sprintf("%s gene expression", 
                                                                                    symbol), title = sprintf("%s_%s", probe, 
                                                                                                             symbol), ...)
      if (save) {
        filename <- sprintf("%s/%s_%s_bypair.pdf", 
                            dir.out, probe, symbol)
        ggsave(filename = filename, plot = P, useDingbats = FALSE, 
               width = width, height = height)
      }
    }
    close(pb)
  }
  if (length(byProbe$probe) != 0) {
    nearGenes <- GetNearGenes(data = data, probes = byProbe$probe, 
                              numFlankingGenes = byProbe$numFlankingGenes)
    for (i in byProbe$probe) {
      probe <- i
      gene <- nearGenes %>% filter(nearGenes$ID == i) %>% 
        pull("GeneID")
      symbol <- getSymbol(data, geneID = gene)
      exp <- assay(getExp(data)[gene, ])
      meth <- assay(getMet(data)[byProbe$probe, ])
      rownames(exp) <- symbol
      P <- scatter(meth = meth, exp = exp, ylim = ylim, 
                   category = category, dots.size = dots.size, legend.title = legend.title, 
                   xlab = sprintf("DNA methylation at %s", 
                                  probe), ylab = sprintf("Gene expression"), 
                   title = sprintf("%s nearby %s genes", probe, 
                                   byProbe$numFlankingGenes), ...)
      if (save) 
        ggsave(filename = sprintf("%s/%s_byprobe.pdf", 
                                  dir.out, probe), plot = P, useDingbats = FALSE, 
               width = width, height = height)
    }
  }
  if (length(byTF$TF) != 0) {
    probes <- byTF$probe[byTF$probe %in% rownames(assay(getMet(data)))]
    meth <- colMeans(assay(getMet(data)[probes, ]), na.rm = TRUE)
    gene <- getGeneID(data, symbol = byTF$TF)
    found <- NULL
    if (any(is.na(gene))) {
      found <- !is.na(gene)
      message("Gene not found: ", byTF$TF[!found])
      gene <- na.omit(gene)
    }
    exp <- assay(getExp(data)[gene, ])
    if (nrow(exp) > 0) {
      if (!is.null(found)) {
        rownames(exp) <- byTF$TF[found]
      }
      else {
        rownames(exp) <- byTF$TF
      }
    }
    P <- scatter(meth = meth, exp = exp, ylim = ylim, category = category, 
                 dots.size = dots.size, correlation = correlation, 
                 legend.title = legend.title, xlab = "Avg DNA methylation", 
                 ylab = sprintf("TF expression"), title = "TF vs avg DNA methylation", 
                 ...)
    if (save) 
      ggsave(filename = sprintf("%s/%s_byTF.pdf", 
                                dir.out, paste(byTF$TF, collapse = "_")), 
             plot = P, useDingbats = FALSE, width = max(6, 
                                                        3 * (length(byTF$TF)%%5)), height = max(4, 
                                                                                                3 * ceiling(length(byTF$TF)/5)))
  }
  return(P)
}