#' Convert input file ID to dataframe class with entrez ID, symbol and ensembl ID
#'
#' @param inputSpecies human, mouse, rat
#' @param inputIDType Human: entrezgene, hgnc_symbol, ensembl_gene_id,
#' Mouse: entrezgene, mgi_symbol, ensembl_gene_id,
#' Rat: entrezgene, rgd_symbol, ensembl_gene_id
#' @param geneIDCol the input ID
#' @return input sample in dataframe class with entrez ID, symbol and ensembl ID
#' @import biomaRt

getAllTypesID<-function(inputSpecies,inputIDType,geneIDCol){
  biomaRtSpecies<-speciesConvert2Biomart(inputSpecies)
  ensembl= biomaRt:: useMart("ensembl")
  ensembl = biomaRt:: useDataset(biomaRtSpecies,mart=ensembl)
  queryData<-as.character(geneIDCol)

   attiArray<-NULL
   if(inputSpecies=="mouse")
   {
     attiArray=
       c('mgi_symbol', 'ensembl_gene_id','entrezgene')}
        else if  (inputSpecies=="human")
        {
          attiArray=
            c('hgnc_symbol', 'ensembl_gene_id','entrezgene')
        }
        else if  (inputSpecies=="Rat"){
           attiArray=
             c('rgd_symbol', 'ensembl_gene_id','entrezgene')
        }

  IDMapping<- biomaRt::getBM(attributes=attiArray,
                    filters = inputIDType,
                    values = queryData,
                    mart = ensembl)

  return (IDMapping)
}

