#' Data from sample file
#'
#' The data from Affymetrix microarray of western diet feed mouse (WD)
#' vs. low fat diet feed mouse (LF) were used as an example.
#' 748 differential expressed genes were identified with
#' cutoff FDR<0.05 using limma R package.
#'
#' @format A data frame with 721 rows and 10 variables:
#' \describe{
#'   \item{entrezgene}{entrez gene ID}
#'   \item{mgi_symbol}{gene symbol}
#'   \item{AveExpr}{average expression value of gene}
#'   \item{t}{t value}
#'   \item{P.Value}{p value}
#'   \item{adj.P.Val}{adjusted p value}
#'   \item{B}{B value}
#'   \item{log2FC}{log fold change value}
#'   \item{ensembl_gene_id}{ensembl gene ID}
#'   \item{KEGG_GeneID}{KEGG gene ID}
#' }
"inputSample"


#' Sample in GSEABase class
#'
#' The data from Affymetrix microarray of western diet feed mouse (WD)
#' vs. low fat diet feed mouse (LF)
#' were used as an example.
#' differential expressed genes were identified with
#' cutoff FDR<0.05 using limma R package.
#'
#' @format Sample in GSEABase class:
#' \describe{
#'   \item{geneIds}{entrez gene ID}
#' }
"gs"
