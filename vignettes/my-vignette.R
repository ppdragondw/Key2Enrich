## ---- include=FALSE------------------------------------------------------
library(Key2Enrich)

## ----echo=TRUE, fig.show='hold', message=FALSE, warning=FALSE------------
# Load sample data.
data(inputSample)
inputSample<-as.data.frame(inputSample)
head(inputSample)

## ----echo=TRUE, fig.show='hold', message=FALSE, warning=FALSE------------

thisKEGGEnricherMatrix<-KEGGEnricherMatrixRealTime(inputSample,"mouse","fdr","fdr",0.05)
colnames(thisKEGGEnricherMatrix)
head(thisKEGGEnricherMatrix[c(1:6,11:13)])

## ----fig.cap="Signicantly enriched KEGG pathways", echo=FALSE, message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/KEGG_barplot.png')

## ----echo=TRUE, fig.show='hold', message=FALSE, warning=FALSE------------

KEGGSigMx<-Key2KEGGEnrich(inputSample=inputSample,inputSpecies="mouse",
                          adjustMethod="fdr",filterMethod="fdr",
                          filterValue=0.05,ifbarplot=TRUE,
                          ifMapOnPath=TRUE,type="KEGG",
                          imgWidth=15,imgHeight=20)
head(KEGGSigMx[c(1:6,11:13)])

## ----echo=FALSE, fig.cap="Signicantly enriched KEGG pathways", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/KEGG_barplot.png')

## ----echo=FALSE, fig.cap="Map differential expressed genes on KEGG pathway", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/KEGG_mapping1.png')
knitr::include_graphics('../inst/extdata/figure/KEGG_mapping2.png')

## ----fig.cap="Reactome enrichment bar plot", echo=TRUE, fig.show='hold', message=FALSE, warning=FALSE, fig.show='hide'----
reactomePValueMatrix<-Reactome2Enrich(geneIDs=as.character(inputSample$entrezgene),"mouse","fdr",0.05)
colnames(reactomePValueMatrix)
head(reactomePValueMatrix[,1:6])

## ----echo=FALSE, fig.cap="Signicantly enriched Reactome pathways", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/reactome_barplot.png')

## ----echo=TRUE, echo=TRUE, message=FALSE, warning=FALSE, out.width='97%',fig.show='hide'----
egoDF<-GO2Enrich(geneIDs=as.character(inputSample$entrezgene),"mouse","fdr",0.05,"CC")
colnames(egoDF)
egoDF[,c(1:3,5:7)]

## ----echo=FALSE, fig.cap="Signicantly enriched Cellular Component terms of gene ontology", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/CC.png')

## ----fig.cap="KEGG enrichment bar plot", echo=TRUE, fig.show='hide', message=FALSE, warning=FALSE----
mappingOnSpecifiedKEGGPathway("mmu04976",inputSample,"mouse")

## ----echo=FALSE, fig.cap="Map differential expressed genes on KEGG pathway", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/KEGG_mapping3.png')

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
mmuInfo<-getAllPathNameAndID("mmu")
head(mmuInfo)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
mmuGeneInfo<-KEGGID2EntrezID("mmu")
head(mmuGeneInfo)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
allMmuPath<-getTotalPathNames("mmu")
head(allMmuPath)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
allGeneInPathwayDF<-getAllGeneInPathwayDF("mmu")
N<-getN(allGeneInPathwayDF)
n<-get_n(inputSample,allGeneInPathwayDF)

pValueMX<-getPValue("mmu00053","mmu",inputSample,N,n)
pValueMX[1:6]

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
pathIDandKEGG_geneID<-getAllGeneInPathwayDF("mmu")
head(pathIDandKEGG_geneID)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
   path<-system.file(package = "Key2Enrich")
   xmlPath<-paste(path,"/extdata",sep="")
   xmlFile<-paste(xmlPath,sep="","/hsa04012.xml")
   r<-parseKGMLFile(xmlFile)
   kegg.nodes <- lapply(r[childIsEntry(r)], parseEntry)
   nodeDF<-data.frame(KEGG_GeneID="NA",type=NA,link=NA,graphicName=NA,
                      graphicType=NA,graphicX=NA,graphicY=NA,
                      graphicWidth=NA,graphicHeight=NA)
    nodeDF<-parseList2Dataframe(kegg.nodes,nodeDF)
    head(nodeDF)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
data(gs)
library(GSEABase)
reactomePValueMatrix<-Reactome2Enrich(geneIds(gs),"mouse","fdr",0.05)

egoDF<-GO2Enrich(geneIds(gs),"mouse","fdr",0.05,"CC")


## ----echo=TRUE-----------------------------------------------------------
sessionInfo()

