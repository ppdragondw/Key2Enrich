# Read sample from file, the input file should in csv format.
# Each line of the file is a data record.
# Each record consists of one or more fields, separated by commas.
# The data from Affymetrix microarray of western diet feed mouse (WD)
# vs. low fat diet feed mouse (LF) were used as an example.
# 748 differential expressed genes were identified with cutoff FDR<0.05 using limma R package.

inputSample<-formatInputSample("./inst/extdata/sampleFile.csv",
                               IDColumn =  1,
                               logFCColumn= 2,
                               IDType="mgi_symbol",
                               inputSpecies = "mouse")

devtools::use_data(inputSample,overwrite = TRUE)
