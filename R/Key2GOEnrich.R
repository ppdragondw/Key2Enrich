#' Perform Gene Ontology enrichment analysis
#'
#' @param inputSample  formatted input sample
#' @param inputSpecies human, mouse, rat
#' @param adjustMethod p.adjust.methods,
#' c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param filterValue must by numeric, value from 0 to 1
#' @param GOType BP,MF,CC
#' @param imgWidth the width of export file
#' @param imgHeight the height of export file
#' @return data in dataframe class with pValue, adjust pValue of significate Gene Ontology
#' @export
#' @examples
#' egoDF<-GO2Enrich(inputSample=inputSample,
#' inputSpecies="mouse",
#' adjustMethod="fdr",
#' filterValue=0.05,
#' GOType="BP",
#' imgWidth=15,
#' imgHeight=20)
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @importFrom  clusterProfiler enrichGO

GO2Enrich<-function (inputSample,
                     inputSpecies,
                     adjustMethod,
                     filterValue,
                     GOType,
                     imgWidth,
                     imgHeight) {

  Db<-NULL
  attiArray<-NULL
  if (inputSpecies=="human")
  {
    Db<-org.Hs.eg.db
  }
  else if (inputSpecies=="mouse")
  {
    Db<-org.Mm.eg.db
  }
  else if (inputSpecies=="rat")
  {
    Db<-org.Rn.eg.db
  }

  if(inputSpecies=="mouse")
  {
    attiArray=
      c('mgi_symbol', 'ensembl_gene_id','entrezgene')
  }
  else if (inputSpecies=="human")
  {
    attiArray=
      c('hgnc_symbol', 'ensembl_gene_id','entrezgene')
  }
  else if  (inputSpecies=="Rat")
  {
    attiArray=
      c('rgd_symbol', 'ensembl_gene_id','entrezgene')
  }

  biomaRtSpecies<-speciesConvert2Biomart(inputSpecies)
  ensembl=biomaRt::useMart("ensembl")
  ensembl = biomaRt::useDataset(biomaRtSpecies,mart=ensembl)
  IDMapping<-biomaRt::getBM(attributes=attiArray,mart = ensembl)

  bg<-as.character(IDMapping$entrezgene)
  bg<-bg[which(bg!="NA")]
  ego <- enrichGO(gene          = inputSample$entrezgene,
                  universe      = bg,
                  OrgDb         = Db,
                  ont           = GOType,
                  pAdjustMethod = adjustMethod,
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

  egoDF<-data.frame(negativeLog=-log10(ego$p.adjust))
  egoDF$pathName<- ego$Description
  egoDF$GeneRatio<-ego$GeneRatio
  egoDF$Genes<-ego$geneID

  splitCol<-strsplit(as.character(egoDF$GeneRatio),'/',fixed=TRUE)
  egoDF<-data.frame(egoDF,do.call(rbind, splitCol))
  egoDF$ratio<-
    as.numeric(as.character(egoDF$X1))/as.numeric(as.character(egoDF$X2))
  colnames(egoDF)[c(5,6)]<- c("m","M")
  egoDF<-egoDF[order(egoDF$negativeLog),]

  thisKEGGBarplot<-Key2EnrichBarplot(egoDF,GOType,imgWidth,imgHeight)

  return (egoDF)
}
