#' Read samples files with different annotation ID
#'
#' @param fileName input file in csv format
#' @param IDColumn must be numberic,the column number of sample ID
#' @param logFCColumn must be numberic,the column number of sample log fold change
#' @param IDType Human: entrezgene, hgnc_symbol, ensembl_gene_id,
#' Mouse: entrezgene, mgi_symbol, ensembl_gene_id,
#' Rat: entrezgene, rgd_symbol, ensembl_gene_id
#' @return Sample file in data.frame format
#' @examples #entrezgene_sample<-readInputFile("sampleFile_entrezgene.csv",1,2,"entrezgene")
#'
readInputFile<-function (fileName,IDColumn,logFCColumn,IDType){
inputSample<-read.csv(fileName)
sample<-subset(inputSample,select=-c(IDColumn,logFCColumn))

sample$ID<-inputSample[,IDColumn]
sample$log2FC<-inputSample[,logFCColumn]

colnames(sample)[which(names(sample)=="ID")]<-IDType

return (sample)
}

#' Format input file to dataframe class with entrez ID, symbol and ensembl ID
#'
#' @param filename input file in csv format
#' @param IDColumn must be numberic,the column number of sample ID
#' @param logFCColumn must be numberic,the column number of sample log fold change
#' @param IDType Human: entrezgene, hgnc_symbol, ensembl_gene_id,
#' Mouse: entrezgene, mgi_symbol, ensembl_gene_id,
#' Rat: entrezgene, rgd_symbol, ensembl_gene_id
#' @param inputSpecies human, mouse, rat
#' @return input sample in dataframe class with entrez ID, symbol and ensembl ID
#' @export
#' @examples #inputSample<-formatInputSample("./inst/extdata/sampleFile.csv",1,2,"entrezgene","mouse")

formatInputSample<-function (filename,IDColumn,logFCColumn,IDType,inputSpecies){
  sample<-readInputFile(filename,IDColumn,logFCColumn,IDType)
  IDColValue<-sample[,(names(sample)==IDType)]
  IDMapping<-getAllTypesID(inputSpecies,IDType,IDColValue)
  inputSampleWith3IDs<-merge(sample,
                             IDMapping,
                             by.x=IDType,
                             by.y=IDType)
  KEGGSpecies<-speciesKEGGConvert(inputSpecies)
  inputSampleKEGG<-merge(inputSampleWith3IDs,
                         KEGGID2EntrezID(KEGGSpecies),
                         by.x="entrezgene",
                         by.y="entrezgene")
  return (inputSampleKEGG)
}

#' Convert input species to biomart species
#'
#' @param inputSpecies human, mouse, rat
#' @return species of biomart format,
#' human:hsapiens_gene_ensembl,
#' mouse:mmusculus_gene_ensembl,
#' rat:rnorvegicus_gene_ensembl
#' @examples #biomartSpecies<-speciesConvert2Biomart("mouse")
#'
speciesConvert2Biomart <- function(inputSpecies) {
  switch(inputSpecies,
         mouse = "mmusculus_gene_ensembl",
         human = "hsapiens_gene_ensembl",
         rat = "rnorvegicus_gene_ensembl")
}

#' Convert input species to KEGG species
#'
#' @param inputSpecies human, mouse, rat
#' @return species of KEGG format, human:hsa, mouse:mmu, rat:rno
#' @examples #KEGGSpecies<-speciesKEGGConvert("mouse")
#'
speciesKEGGConvert <- function(inputSpecies) {
  switch(inputSpecies,
         mouse = "mmu",
         human = "hsa",
         rat = "rno")
}
#' Get gene symbol from different species
#'
#' @param thisKEGGSpecies species in KEGG format
#' @return name of gene symbol
#' @examples #geneSymbol<-getGeneSymbol("mouse")
#'
getGeneSymbol <- function(thisKEGGSpecies) {
    switch(thisKEGGSpecies,
    mmu = "mgi_symbol",
    hsa = "hgnc_symbol",
    rno = "rgd_symbol")
}

#' Convert input species to KEGG species with flag
#'
#' @param inputSpecies human, mouse, rat
#' @return species of KEGG format with flag,
#' human:- Homo sapiens,
#' mouse:- Mus musculus,
#' rat:- Rattus norvegicus
#' @examples #KEGGSpeciesFlag<-speciesKEGGFlagConvert("mouse")
#'
speciesKEGGFlagConvert	 <- function(inputSpecies){
  #mouse:"- Mus musculus"
  #human:"- Homo sapiens"
  #rat:"- Rattus norvegicus"
  switch(inputSpecies,
         mmu = "- Mus musculus",
         hsa = "- Homo sapiens",
         rno = "- Rattus norvegicus")
}

#' Mapping input file with ID to symbol, entrez, ensembl IDs
#'
#' @param filename input file in csv format
#' @param IDColumn must be numberic,the column number of sample ID
#' @param logFCColumn must be numberic,the column number of sample log fold change
#' @param IDType Human: entrezgene, hgnc_symbol, ensembl_gene_id,
#' Mouse: entrezgene, mgi_symbol, ensembl_gene_id,
#' Rat: entrezgene, rgd_symbol, ensembl_gene_id
#' @param inputSpecies human, mouse, rat
#' @return sample with entrez ID, symbol and ensembl ID
#' @examples #thisSampleWith3IDs<-sampleWith3IDs("sampleFile_entrezgene.csv",1,2,"entrezgene","mouse")
#'
sampleWith3IDs<-function (filename,IDColumn,logFCColumn,IDType,inputSpecies){
#inputSpecies mouse,human,rat
  sample<-readInputFile(filename,IDColumn,logFCColumn,IDType)
  #KEGGSpecies<-speciesKEGGConvert(inputSpecies)

  IDColValue<-sample[,(names(sample)==IDType)]

  IDMapping<-getAllTypesID(inputSpecies,IDType,IDColValue)
  inputSampleWith3IDs<-merge(sample,IDMapping,by.x=IDType,by.y=IDType)

  return (inputSampleWith3IDs)
}
