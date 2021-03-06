
#' Get enrichment matrix of enriched Reactome pathway
#'
#' @param geneIDs  gene identifiers
#' @param inputSpecies human, mouse, rat
#' @param adjustMethod p.adjust.methods,
#' c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param filterValue must by numeric, value from 0 to 1
#' @param imgWidth the width of export file
#' @param imgHeight the height of export file
#' @return data in dataframe class with pValue, adjust pValue of significate Reactome pathway
#' @export
#' @examples
#' data(inputSample)
#' inputSample<-as.data.frame(inputSample)
#' head(inputSample)
#' \donttest{
#' reactomePValueMatrix<-Reactome2Enrich(geneIDs=as.character(inputSample$entrezgene),
#'                                       "mouse","fdr",0.05,15,20)
#' }
#' @importFrom  ReactomePA enrichPathway
Reactome2Enrich<-function (geneIDs,inputSpecies,
                           adjustMethod,filterValue,
                           imgWidth,imgHeight) {
reactomePValueMatrix<- enrichPathway(gene=geneIDs,
                                    qvalueCutoff=filterValue,
                                    readable=TRUE,
                                    pAdjustMethod=adjustMethod,
                                    organism = inputSpecies)

reactomePValueDF<-
  data.frame(negativeLog=-log10(reactomePValueMatrix$p.adjust))
reactomePValueDF$pathName<-
  reactomePValueMatrix$Description
reactomePValueDF$GeneRatio<-reactomePValueMatrix$GeneRatio
splitCol<-strsplit(as.character(reactomePValueDF$GeneRatio),'/',fixed=TRUE)
reactomePValueDF<-
  data.frame(reactomePValueDF,do.call(rbind, splitCol))
reactomePValueDF$ratio<-
  as.numeric(as.character(reactomePValueDF$X1))/as.numeric(as.character(reactomePValueDF$X2))
colnames(reactomePValueDF)[c(4,5)]<- c("m","M")

reactomePValueDF$genes<-reactomePValueMatrix$geneID

reactomePValueDF<-reactomePValueDF[order(reactomePValueDF$negativeLog),]

thisKey2EnrichBarplot<-
  Key2EnrichBarplot(reactomePValueDF,"Reactome",imgWidth,imgHeight)

return (reactomePValueDF)
}
