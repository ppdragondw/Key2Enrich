

context("enrichment analysis output")

data(inputSample)
inputSample<-as.data.frame(inputSample)

test_that("Output of Reactome pathway enrichment analysis", {
  reactomePValueMatrix<-Reactome2Enrich(inputSample$entrezgene,"mouse","fdr",0.05,15,20)
  expect_is (reactomePValueMatrix$ratio,"numeric")
})

test_that("Output of GO enrichment analysis", {
  egoDF<-GO2Enrich(inputSample$entrezgene,
                   inputSpecies="mouse",adjustMethod="fdr",
                   filterValue=0.05,GOType="BP",
                   imgWidth=15,imgHeight=20)
  expect_is (egoDF$ratio,"numeric")
})
