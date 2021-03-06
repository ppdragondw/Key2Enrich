% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KEGG.R
\name{get_inputSampleKEGG_m}
\alias{get_inputSampleKEGG_m}
\title{Get sample genes in specific KEGG pathway}
\usage{
get_inputSampleKEGG_m(thisInputSampleKEGG, geneInOnePathwayDF)
}
\arguments{
\item{thisInputSampleKEGG}{input sample with in KEGG gene ID}

\item{geneInOnePathwayDF}{gene list in specific KEGG pathway}
}
\value{
sample genes in specific KEGG pathway
}
\description{
Get sample genes in specific KEGG pathway
}
\examples{
data(inputSample)
inputSample<-as.data.frame(inputSample)
\donttest{
geneInOnePathwayDF<-getGeneInOnePathwayDF("mmu05160","mmu")
get_inputSampleKEGG_m(inputSample,geneInOnePathwayDF)
}
}
