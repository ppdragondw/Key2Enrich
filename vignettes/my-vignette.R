## ---- include=FALSE------------------------------------------------------
library(Key2Enrich)

## ----echo=TRUE, fig.show='hold', message=FALSE, warning=FALSE------------
# Load sample data.
data(inputSample)
inputSample<-as.data.frame(inputSample)
head(inputSample)

## ----fig.cap="Signicantly enriched KEGG pathways", echo=FALSE, message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/KEGG_barplot.png')

## ----echo=FALSE, fig.cap="Signicantly enriched KEGG pathways", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/KEGG_barplot.png')

## ----echo=FALSE, fig.cap="Map differential expressed genes on KEGG pathway", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/KEGG_mapping1.png')
knitr::include_graphics('../inst/extdata/figure/KEGG_mapping2.png')

## ----echo=FALSE, fig.cap="Signicantly enriched Reactome pathways", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/reactome_barplot.png')

## ----echo=FALSE, fig.cap="Signicantly enriched Cellular Component terms of gene ontology", message=FALSE, warning=FALSE, out.width='97%'----
knitr::include_graphics('../inst/extdata/figure/CC.png')

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
data(gs)
library(GSEABase)
reactomePValueMatrix<-Reactome2Enrich(geneIds(gs),"mouse","fdr",0.05)

egoDF<-GO2Enrich(geneIds(gs),"mouse","fdr",0.05,"CC")


## ----echo=TRUE-----------------------------------------------------------
sessionInfo()

