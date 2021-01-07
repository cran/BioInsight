#' @export
#' @import knitr
#' @import RColorBrewer
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics par
#' @importFrom stats dist hclust
#' @importFrom edgeR filterByExpr
#' @importFrom limma plotMDS
#' @importFrom wordcloud wordcloud
describeRNA = function(counts, biotypes, groups, report=FALSE, verbose=FALSE, filter=1)
{
  table = table(biotypes$gene_biotype, exclude=NA)
  target = c("miRNA", "protein_coding", "lincRNA","pseudogene","snoRNA","snRNA", "ribozyme")
  table = data.frame(table)
  index <- table$Var1 %in% target 
  barplot = table[index, ] 
  if (report) {
    File <- tempfile(fileext=".pdf") 
    warning("\n Temporary report at  ", File, call.=FALSE, immediate. =TRUE) 
    dir.create(dirname(File))
    pdf(File,width=15, height=15)
    oldpar <- par(no.readonly=TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(2,2))
    plotMDS(counts, main="Multidimensional Scaling")
    sortbar = barplot[order(barplot$Freq, decreasing=TRUE), ]
    barplot(sortbar$Freq, names.arg=sortbar$Var1, col=1:6, main="Absolute Quantity")
    x = t(counts)
    x = hclust(dist(x))
    plot(x, main="Cluster Dendrogram")
    wordcloud(table$Var1, table$Freq, colors=brewer.pal(5, "Dark2"), min.freq = 10)
    dev.off()
  }
  if (filter == 1) {
    data_filtered = filterByExpr(counts, group=groups)
    data_filtered = counts[data_filtered, ]
  }
  if (filter == 2){
    data_filtered = filterByExpr(counts, group=groups, min.count=15, min.total.count=25)
    data_filtered = counts[data_filtered, ]
  }
  if (filter == 3){
    data_filtered = filterByExpr(counts, group=groups, min.count=25, min.total.count=40)
    data_filtered = counts[data_filtered, ]
  }
  if (verbose){
    print(knitr::kable(table))
    cat("\nGENES", sep='\n')
    cat("Total number of genes:", nrow(counts))
    cat("\nGenes remaining:", dim(data_filtered)[1])
  }
  pos <- 1
  envir = as.environment(pos)
  BioInsight <- assign("BioInsight", data_filtered, envir = envir)
}
