\name{describeRNA}
\alias{describeRNA}
\title{Biotype Summary}
\usage{
describeRNA(counts, biotypes, groups, report=FALSE, verbose=FALSE, filter=1)
}
\description{
This function provide filterByExpr with three custom options. See 'Arguments'.
Different strategies can be useful either if you are trying to compare different approaches of normalization or if you want to analyze a particular biotype where a different variation of expression is expected under certain conditions.
A BioInsight data frame will return with your new count matrix where you can proceed with your Differential Expression Analysis or revealRNA() to produce new plots and tables after the filter.
}
\arguments{
  \item{counts}{Data.frame of counts matrix}
  \item{biotypes}{The count matrix with a gene_biotype column. See 'Examples'}
  \item{groups}{Groups as factor}
  \item{report}{If TRUE, a .pdf file formated in a proportion of 2:2 with a Multidimensional Scalling,barplot, dendrogram and  wordcloud will be generated in a tempdir}
  \item{verbose}{if TRUE, a count table will be printed in the console with the numbers of biotypes}
  \item{filter}{Numeric from 1 to 3. See 'Details'}
}
\details{
1 "DEFAULT" from edgeR: Min.count=10, min.total.count=15;
2 "Slightly" above the DEFAULT: Min.count=15, min.total.count=25;
3 "Restricted" than DEFAULT: Min.count=25, min.total.count=40;
}
\references{
  Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26(1):139-140. doi:10.1093/bioinformatics/btp616
}
\note{
  See trace(describeRNA, edit=T) to modify the values of filterByExpr and/or the biotypes of interest.
}
\examples{
counts = system.file("extdata", "count_matrix.tsv", package = "BioInsight")
counts = read.table(counts[1], header=TRUE, row.names=1)
biotypes = system.file("extdata", "Rattus_Norvegicus_biomart.tsv", package = "BioInsight")
biotypes = read.table(biotypes[1], header=TRUE, row.names=1)

groups = rep(as.factor(c("1","2")), each=5)

describeRNA(counts=counts,
            biotypes=biotypes,
            groups=groups,
            filter=2)

#Annotation - see biomaRt
#mart <- useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl")
#genes = getBM(attributes=c("ensembl_gene_id","external_gene_name",
#                           "gene_biotype", "start_position",
#                           "end_position", "source"),
#              values = row.names(data),
#              mart = mart)
}
