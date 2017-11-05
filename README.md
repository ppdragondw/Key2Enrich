# Key2Enrich
An all-in-one R/Bioconductor package for gene list enrichment analysis and pathway visualization


## Abstract
Enrichment analysis is a popular method for analyzing gene sets generated from transcriptomes and proteomics. Here, we developed Key2Enrich, a novel all-in-one R package for gene set enrichment analysis. It combines the whole pipeline of enrichment analysis in one package, provides enrichment function and visualization methods, while it also maps and renders input sample data on relevant pathway graphs. Key2Enrich automatically downloads the pathway graph data, parses the data file, maps and integrates the genes of input data onto the pathway and renders pathway graphs. Currently, Key2Enrich supports three species, including human, mice and rat. Both enrichment analysis and pathway mapping functions can be integrated with other functional annotation tools for large-scale and fully automated analysis pipelines.


## Get started

Read sample from file, the input file should in csv format. Each line of the file is a data record. Each record consists of one or more fields, separated by commas. The data from Affymetrix microarray of western diet feed mouse (WD) vs. low fat diet feed mouse (LF) were used as an example. 748 differential expressed genes were identified with cutoff FDR<0.05 using limma R package. 
```{r echo=TRUE, fig.show='hold', message=FALSE, warning=FALSE}
  devtools::load_all()
# Read file and add gene symbol, entrez and ensembol IDs
 inputSample<-formatInputSample("../data-raw/sampleMgisymbol.csv", 
                                IDColumn =  1,logFCColumn= 2,
                                IDType="mgi_symbol", inputSpecies = "mouse")
 head(inputSample)
```

KEGG pathway enrichment analysis, this step will take some time, the data is getting from KEGG database by using KEGG API.
```{r, echo=TRUE, message=FALSE, warning=FALSE}
thisKEGGEnricherMatrix<-KEGGEnricherMatrixRealTime(inputSample,"mouse","fdr","fdr",0.05)
colnames(thisKEGGEnricherMatrix)
head(thisKEGGEnricherMatrix[c(1:6,11:13)])

```

Then the result can be show in bar plot.
```{r echo=FALSE, eval=FALSE, message=FALSE, warning=FALSE, include=FALSE}
Key2EnrichBarplot(thisKEGGEnricherMatrix,"KEGG",15,20)
```

![alt text](https://github.com/ppdragondw/dongFile/blob/master/KEGG_barplot.png)


These step can be done by using one function. The figure sizes have been customised so that you can easily put two images side-by-side. The genes can mapped to the nodes of KEGG pathway by using the parameter ifMapOnPath=TRUE


```{r echo=TRUE, message=FALSE, warning=FALSE,fig.show='hide'}
KEGGSigMx<-Key2KEGGEnrich(inputSampleKEGG=inputSample,inputSpecies="mouse",
                          adjustMethod="fdr",filterMethod="fdr",filterValue=0.05,
                          ifbarplot=TRUE,ifMapOnPath=TRUE,type="KEGG",
                          imgWidth=15,imgHeight=20)
head(KEGGSigMx[c(1:6,11:13)])
```
![alt text](https://github.com/ppdragondw/dongFile/blob/master/KEGG_barplot.png)

![alt text](https://github.com/ppdragondw/dongFile/blob/master/KEGG_mapping1.png)
![alt text](https://github.com/ppdragondw/dongFile/blob/master/KEGG_mapping2.png)


Enrichment analysis based on Reactome pathways.
```{r fig.cap="Reactome enrichment bar plot", echo=TRUE, fig.show='hold', message=FALSE, warning=FALSE, fig.show='hide'}
reactomePValueMatrix<-Reactome2Enrich(inputSample,"mouse","fdr",0.05)
colnames(reactomePValueMatrix)
head(reactomePValueMatrix[,1:6])
```

![alt text](https://github.com/ppdragondw/dongFile/blob/master/reactome_barplot.png)


Enrichment analysis based on terms of gene ontology.
```{r echo=TRUE, echo=TRUE, message=FALSE, warning=FALSE, out.width='97%',fig.show='hide'}
egoDF<-GO2Enrich(inputSample,"mouse","fdr",0.05,"CC")
colnames(egoDF)
egoDF[,c(1:3,5:7)]
```

![alt text](https://github.com/ppdragondw/dongFile/blob/master/GO_CC_barplot.png)


## Functions can be called by other packages
The gene list can be mapped to specific KEGG pathway. Key2Enrich can intergate with other genome-wide data analysis tools or packages.
```{r fig.cap="KEGG enrichment bar plot", echo=TRUE, fig.show='hide', message=FALSE, warning=FALSE}
mappingOnSpecifiedKEGGPathway("mmu04976",inputSample,"mouse")
```

Get KEGG pathname and pathway ID of a specific species.
```{r echo=TRUE, message=FALSE, warning=FALSE}
mmuInfo<-getAllPathNameAndID("mmu")
head(mmuInfo)
```
Get all KEGG gene ID and entrezgene ID of a specific species.
```{r echo=TRUE, message=FALSE, warning=FALSE}
mmuGeneInfo<-KEGGID2EntrezID("mmu")
head(mmuGeneInfo)
```

Get all KEGG pathway name of a specific species.
```{r echo=TRUE, message=FALSE, warning=FALSE}
allMmuPath<-getTotalPathNames("mmu")
head(allMmuPath)
```

Calculate p-value and adjust p-value for a specific KEGG pathway on input sample.
```{r echo=TRUE, message=FALSE, warning=FALSE}
allGeneInPathwayDF<-getAllGeneInPathwayDF("mmu")
N<-getN(allGeneInPathwayDF)
n<-get_n(inputSample,allGeneInPathwayDF)

pValueMX<-getPValue("mmu00053","mmu",inputSample,N,n)
pValueMX[1:6]
```
Get all KEGG pathway ID and its corresponding gene list in KEGG gene ID format.
```{r echo=TRUE, message=FALSE, warning=FALSE}
pathIDandKEGG_geneID<-getAllGeneInPathwayDF("mmu")
head(pathIDandKEGG_geneID)
```

