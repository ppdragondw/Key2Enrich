% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{sampleWith3IDs}
\alias{sampleWith3IDs}
\title{Mapping input file with ID to symbol, entrez, ensembl IDs}
\usage{
sampleWith3IDs(filename, IDColumn, logFCColumn, IDType, inputSpecies)
}
\arguments{
\item{filename}{input file in csv format}

\item{IDColumn}{must be numberic,the column number of sample ID}

\item{logFCColumn}{must be numberic,
the column number of sample log fold change}

\item{IDType}{Human: entrezgene, hgnc_symbol, ensembl_gene_id,
Mouse: entrezgene, mgi_symbol, ensembl_gene_id,
Rat: entrezgene, rgd_symbol, ensembl_gene_id}

\item{inputSpecies}{human, mouse, rat}
}
\value{
sample with entrez ID, symbol and ensembl ID
}
\description{
Mapping input file with ID to symbol, entrez, ensembl IDs
}
\examples{
path<-system.file(package = "Key2Enrich")
filePath<-paste(path,"/extdata",sep="")
file<-paste(filePath,sep="","/sampleFile.csv")
\donttest{
thisSampleWith3IDs<-sampleWith3IDs(file,
                                   IDColumn=1,
                                   logFCColumn=2,
                                   IDType="mgi_symbol",
                                   inputSpecies="mouse")
}
}
