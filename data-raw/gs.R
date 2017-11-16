# Read sample from file, the input file should in csv format.
# Each line of the file is a data record.
# Each record consists of one or more fields, separated by commas.
# The data from Affymetrix microarray of western diet feed mouse (WD)
# vs. low fat diet feed mouse (LF) were used as an example.
# 748 differential expressed genes were identified with cutoff FDR<0.05 using limma R package.
# gs is a sample in GSEABase class

library(GSEABase)

inputSample<-formatInputSample("./inst/extdata/sampleFile.csv",
                               IDColumn =  1,
                               logFCColumn= 2,
                               IDType="mgi_symbol",
                               inputSpecies = "mouse")
noDuplicatedEntrezgene<-as.character(inputSample[!duplicated(inputSample$entrezgene),]$entrezgene)
gs <- GeneSet(geneIds=noDuplicatedEntrezgene,
              setName="GSEABaseSample",
              geneIdType=EntrezIdentifier(),
              collectionType=ExpressionSetCollection())
devtools::use_data(gs,overwrite = TRUE)


