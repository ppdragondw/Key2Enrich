## ---- include=FALSE------------------------------------------------------
devtools::load_all()

## ----echo=TRUE, fig.show='hold', message=FALSE, warning=FALSE------------
# Load sample data.
devtools::load_all()
data(inputSample)
inputSample<-as.data.frame(inputSample)
head(inputSample)

## ----fig.cap="Signicantly enriched KEGG pathways", echo=FALSE, message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('./figure/KEGG_barplot.png')

## ----echo=FALSE, fig.cap="Signicantly enriched KEGG pathways", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('./figure/KEGG_barplot.png')

## ----echo=FALSE, fig.cap="Map differential expressed genes on KEGG pathway", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('./figure/KEGG_mapping1.png')
knitr::include_graphics('./figure/KEGG_mapping2.png')

## ----fig.cap="Reactome enrichment bar plot", echo=TRUE, fig.show='hold', message=FALSE, warning=FALSE, fig.show='hide'----
reactomePValueMatrix<-Reactome2Enrich(inputSample,"mouse","fdr",0.05)
colnames(reactomePValueMatrix)
head(reactomePValueMatrix[,1:6])

## ----echo=FALSE, fig.cap="Signicantly enriched Reactome pathways", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('./figure/reactome_barplot.png')

## ----echo=TRUE, echo=TRUE, message=FALSE, warning=FALSE, out.width='97%',fig.show='hide'----
egoDF<-GO2Enrich(inputSample,"mouse","fdr",0.05,"CC")
colnames(egoDF)
egoDF[,c(1:3,5:7)]

## ----echo=FALSE, fig.cap="Signicantly enriched Cellular Component terms of gene ontology", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('./figure/CC.png')

## ----fig.cap="KEGG enrichment bar plot", echo=TRUE, fig.show='hide', message=FALSE, warning=FALSE----
mappingOnSpecifiedKEGGPathway("mmu04976",inputSample,"mouse")

## ----echo=FALSE, fig.cap="Map differential expressed genes on KEGG pathway", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('./figure/KEGG_mapping3.png')

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

## ----echo=TRUE-----------------------------------------------------------
sessionInfo()

