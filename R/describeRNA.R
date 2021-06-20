##' @title describeRNA
##' @description This function provide filterByExpr with two customized options. See 'Arguments'. Different strategies can be useful either if you are trying to compare different approaches of normalization or if you want to analyze a particular biotype where a different variation of expression is expected under certain conditions. A BioInsight data frame will return with your new count matrix where you can proceed with your Differential Expression Analysis.
##' @param counts data.frame where you have gene counts
##' @param biotypes data.frame where you have a gene_biotype column
##' @param groups factor with groups and samples (see Examples)
##' @param report if TRUE will generate a .pdf file in /tmp/ folder ('wordcloud' and 'RColorBrewer' are necessary)
##' @param verbose if TRUE will print your result in your console ('knitr' is necessary)
##' @param filter Your threshold (see edgeR::filterByExpr)
##' @import edgeR wordcloud knitr RColorBrewer grDevices graphics stats limma
##' @details For filter you may want to use "1" if you want the default option from filterByExpr. If 2 "Slightly" above the default will be applied. If "3" A more restrictive option will be applied. 
##' @note You can use trace(describeRNA, edit=T) to set different values as threshold for "filter" option.
##' @references Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26(1):139-140. doi:10.1093/bioinformatics/btp616
##' @examples 
##' counts = system.file("extdata", "count_matrix.tsv", package="BioInsight")
##' counts = read.table(counts, row.names=1, header=TRUE)
##' biotypes = system.file("extdata", "Rattus_Norvegicus_biomart.tsv", package="BioInsight")
##' biotypes = read.table(biotypes, row.names=1, header=TRUE)
##' 
##' groups = rep(as.factor(c("1","2")), each=5)
##' 
##' describeRNA(counts=counts,
##'             biotypes=biotypes,
##'             groups=groups,
##'             filter=2)
#' @export describeRNA
describeRNA = function(counts, biotypes, groups, report = FALSE, verbose = FALSE, filter = 1) 
{
  table = table(biotypes$gene_biotype, exclude = NA)
  target = c("miRNA", "protein_coding", "lincRNA", "pseudogene", 
             "snoRNA", "snRNA", "ribozyme")
  table = data.frame(table)
  index <- table$Var1 %in% target
  barplot = table[index, ]
  if (report) {
    File <- tempfile(fileext = ".pdf")
    warning("\n Temporary report at  ", File, call. = FALSE, 
            immediate. = TRUE)
    dir.create(dirname(File), showWarnings=FALSE)
    pdf(File, width = 15, height = 15)
    oldpar <- par(no.readonly = TRUE)
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
  if (filter == 1) {
    data_filtered = filterByExpr(counts, group = groups)
    data_filtered = counts[data_filtered, ]
  }
  if (filter == 2) {
    data_filtered = filterByExpr(counts, group = groups, 
                                 min.count = 15, min.total.count = 25)
    data_filtered = counts[data_filtered, ]
  }
  if (filter == 3) {
    data_filtered = filterByExpr(counts, group = groups, 
                                 min.count = 25, min.total.count = 40)
    data_filtered = counts[data_filtered, ]
  }
  if (verbose) {
    print(knitr::kable(table))
    cat("\nGENES", sep = "\n")
    cat("Total number of genes:", nrow(counts))
    cat("\nGenes remaining:", dim(data_filtered)[1])
  }
  pos <- 1
  envir = as.environment(pos)
  BioInsight <- assign("BioInsight", data_filtered, envir = envir)
}
