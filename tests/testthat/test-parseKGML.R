file<-system.file("extdata", "hsa04012.xml", package = "Key2Enrich")

test_that("Output of Reactome pathway enrichment analysis", {
  r<-parseKGMLFile(file)
  kegg.nodes <- lapply(r[childIsEntry(r)], parseEntry)
  nodeDF<-data.frame(KEGG_GeneID="NA",
                     type=NA,
                     link=NA,
                     graphicName=NA,
                     graphicType=NA,
                     graphicX=NA,
                     graphicY=NA,
                     graphicWidth=NA,
                     graphicHeight=NA)
  nodeDF<-parseList2Dataframe(kegg.nodes,nodeDF)
  expect_true(nrow(nodeDF)!=0)
})
