#' @export
#' @import knitr
#' @import RColorBrewer
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics par
#' @importFrom stats dist hclust
#' @importFrom edgeR filterByExpr
#' @importFrom limma plotMDS
#' @importFrom wordcloud wordcloud
revealRNA = function (counts, biotypes, report = FALSE, verbose = FALSE)
{
  counts = counts[grepl("^Sample", names(counts))]
  table = table(biotypes$data.gene_biotype)
  target = c("miRNA", "protein_coding", "lincRNA", "pseudogene",
             "snoRNA", "snRNA", "ribozyme")
  table = data.frame(table)
  index <- table$Var1 %in% target
  barplot = table[index, ]
  if (report) {
    File <- tempfile(fileext = ".pdf")
    warning("Temporary report at", File)
    dir.create(dirname(File))
    pdf(File, width = 15, height = 15)
    oldpar <- par(no.readonly=TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(2, 2))
    plotMDS(counts, main = "Multidimensional Scaling")
    sortbar = barplot[order(barplot$Freq, decreasing = TRUE),
    ]
    barplot(sortbar$Freq, names.arg = sortbar$Var1, col = 1:6,
            main = "Absolute Quantity")
    x = t(counts)
    x = hclust(dist(x))
    plot(x, main = "Cluster Dendrogram")
    wordcloud(table$Var1, table$Freq, colors = brewer.pal(5,"Dark2"), min.freq = 10)
    dev.off()
  }
  if (verbose) {
    cat("RAW TABLE")
    print(knitr::kable(table))
  }
  pos <- 1
  envir = as.environment(pos)
  BioInsight <- assign("BioInsight", BioInsight, envir = envir)
}

