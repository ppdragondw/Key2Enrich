#' Get enrichment matrix of enriched KEGG pathway
#'
#' @param inputSampleKEGG input formatted sample
#' @param inputSpecies human, mouse, rat
#' @param adjustMethod p.adjust.methods, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param filterMethod pValue or other character vector
#' @param filterValue must by numeric, value from 0 to 1
#' @return data in dataframe class with pValue, adjust pValue of significate KEGG pathway
#' @export
#' @examples thisKEGGEnricherMatrix<-KEGGEnricherMatrixRealTime(inputSample,"mouse","fdr","fdr",0.05)

KEGGEnricherMatrixRealTime<-function(inputSampleKEGG,inputSpecies,adjustMethod,filterMethod,filterValue){
  KEGGSpecies<-speciesKEGGConvert(inputSpecies)
  allGeneInPathwayDF<-getAllGeneInPathwayDF(KEGGSpecies)

  N<-getN(allGeneInPathwayDF)
  n<-get_n(inputSampleKEGG,allGeneInPathwayDF)
  pvalueMatrix<-apply(getTotalPathNames(KEGGSpecies),
                    1,
                    getPValue,
                    thisKEGGSpecies=KEGGSpecies,
                    thisInputSampleKEGG=inputSampleKEGG,
                    N=N,
                    n=n)

  pvalueDF<-data.frame(pathID=unlist(sapply(pvalueMatrix,"[",1))
                     ,M=unlist(sapply(pvalueMatrix,"[",2))
                     ,m=unlist(sapply(pvalueMatrix,"[",3))
                     ,N=unlist(sapply(pvalueMatrix,"[",4))
                     ,n=unlist(sapply(pvalueMatrix,"[",5))
                     ,pValue=unlist(sapply(pvalueMatrix,"[",6))
                     ,KEGGGenes=unlist(sapply(pvalueMatrix,"[",7))
                     ,entrezGenes=unlist(sapply(pvalueMatrix,"[",8))
                     ,genes=unlist(sapply(pvalueMatrix,"[",9))
                     ,ensemblGenes=unlist(sapply(pvalueMatrix,"[",10)))

pvalueDF$adjpValue<-p.adjust(pvalueDF$pValue,method = adjustMethod)
filterValueDF<-NULL
if (filterMethod=="pValue"){
  filterValueDF<-pvalueDF[pvalueDF$pValue<filterValue,]
  filterValueDF$negativeLog<- -log10(filterValueDF$pValue)
}
else
{
  filterValueDF<-pvalueDF[pvalueDF$adjpValue<filterValue,]
  filterValueDF$negativeLog<- -log10(filterValueDF$adjpValue)
}
allPathNameAndID<-getAllPathNameAndID(KEGGSpecies)
filterValuePathNameDF<-merge(filterValueDF,allPathNameAndID,by.x="pathID",by.y="pathID")
filterValuePathNameDF <- filterValuePathNameDF[order(filterValuePathNameDF$negativeLog),]
filterValuePathNameDF<- filterValuePathNameDF[filterValuePathNameDF$pathID!="mmu01100",]
return (filterValuePathNameDF)
}

#' Get matrix of enriched KEGG pathway, plot barplot of significantly enriched KEGG pathway, mapping genes on KEGG pathway image
#'
#' @param inputSampleKEGG input formatted sample
#' @param inputSpecies human, mouse, rat
#' @param adjustMethod p.adjust.methods, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param filterMethod pValue or other character vector
#' @param filterValue must by numeric, value from 0 to 1
#' @param ifbarplot barplot or not, TRUE/FALSE
#' @param ifMapOnPath pathway mapping or not, TRUE/FALSE
#' @param mapping gene on pathway or not, TRUE/FALSE
#' @param type KEGG, Reactome, BP,MF,CC
#' @param imgWidth the width of export file
#' @param imgHeight the height of export file
#' @return data in dataframe class with pValue, adjust pValue of significate KEGG pathway
#' @export
#' @examples KEGGSigMx<-Key2KEGGEnrich(inputSampleKEGG,"mouse","fdr","fdr",0.05,TRUE,TRUE,"KEGG",15,20)

Key2KEGGEnrich<-function (inputSampleKEGG,inputSpecies,adjustMethod,filterMethod,filterValue,ifbarplot,ifMapOnPath,type,imgWidth,imgHeight){
  mx<-KEGGEnricherMatrixRealTime(inputSampleKEGG,inputSpecies,adjustMethod,filterMethod,filterValue)
  if (ifbarplot) {Key2EnrichBarplot(mx,type,imgWidth,imgHeight)}
  if (ifMapOnPath) {o<-apply(mx,1,plotSigImg,inputSampleKEGG=inputSampleKEGG,inputSpecies=inputSpecies)}
  return (mx)
}

#' Download KEGG pathway files and mapping genes on pathway
#'
#' @param plotPathID KEGG pathway ID
#' @param inputSampleKEGG input formatted sample
#' @param inputSpecies human, mouse, rat
#' @return NA
#' @export
#' @examples mappingOnSpecifiedKEGGPathway("mmu05160",inputSample,"mouse")

mappingOnSpecifiedKEGGPathway<-function (plotPathID,inputSampleKEGG,inputSpecies){
  #plotPathID:specific KEGG pathway ID
  o<-plotSigImg(plotPathID,inputSampleKEGG,inputSpecies)
}

