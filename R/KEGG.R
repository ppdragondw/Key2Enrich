#' Get all KEGG pathway name and ID of species
#'
#' @param KEGGSpecies species in KEGG format
#' @return KEGG pathway ID of input species
#' @export
#' @import KEGGREST
#' @examples
#' library(Key2Enrich)
#' getAllPathNameAndID("mmu")
getAllPathNameAndID<-function (KEGGSpecies){
  pathways <- keggList("pathway",KEGGSpecies)
  pathways<-as.data.frame(pathways)
  pathways$pathName<-pathways$pathways
  pathways$pathID<-rownames(pathways)
  row.names(pathways)<-NULL
  pathways$pathways<-NULL

  splitCol1<-strsplit(as.character(pathways$pathID),':',fixed=TRUE)
  pathways<-data.frame(pathways,do.call(rbind, splitCol1))

  #splitCol2<-strsplit(as.character(pathways$pathName),"- Mus musculus",fixed=TRUE).
  splitCol2<-strsplit(as.character(pathways$pathName),
                      speciesKEGGFlagConvert(KEGGSpecies),fixed=TRUE)

  pathways<-data.frame(pathways,do.call(rbind, splitCol2))

  pathways<-pathways[c(5,4)]
  colnames(pathways)<-c("pathName","pathID")

  return (pathways)
}

#' Get all KEGG gene ID and entrezgene ID in a specific species
#'
#' @param KEGGSpecies species in KEGG format
#' @return all KEGG Gene ID and its entrezgene ID of input KEGG species
#' @export
#' @examples
#' library(Key2Enrich)
#' KEGGID2EntrezID("mmu")

KEGGID2EntrezID<-function(KEGGSpecies){
  all_KEGG_GeneIDtoGeneID<-keggConv("ncbi-geneid", KEGGSpecies)
  all_KEGG_GeneIDtoGeneIDDF<-data.frame(names(all_KEGG_GeneIDtoGeneID),all_KEGG_GeneIDtoGeneID)
  splitCol<-strsplit(as.character(all_KEGG_GeneIDtoGeneIDDF$names.all_KEGG_GeneIDtoGeneID.),':',fixed=TRUE)
  all_KEGG_GeneIDtoGeneIDDF<-data.frame(all_KEGG_GeneIDtoGeneIDDF,do.call(rbind, splitCol))
  all_KEGG_GeneIDtoGeneIDDF<-all_KEGG_GeneIDtoGeneIDDF[,-c(2,3)]
  colnames(all_KEGG_GeneIDtoGeneIDDF)<-c("KEGG_GeneID","entrezgene")
  rownames(all_KEGG_GeneIDtoGeneIDDF)<- NULL
  return (all_KEGG_GeneIDtoGeneIDDF)
}

#' Get all KEGG pathway name for species
#'
#' @param KEGGSpecies species in KEGG format
#' @return all KEGG pathway name of input KEGG species
#' @export
#' @examples
#' library(Key2Enrich)
#' getTotalPathNames("mmu")
#'
getTotalPathNames<-function (KEGGSpecies){
  pathways <- keggList("pathway",KEGGSpecies)
  KEGGPathwayIDs<-data.frame(names(pathways))
  splitCol<-strsplit(as.character(KEGGPathwayIDs$names.pathways.),':',fixed=TRUE)
  KEGGPathwayIDs<-data.frame(KEGGPathwayIDs,do.call(rbind, splitCol))
  KEGGPathwayIDs<-data.frame(KEGGPathwayIDs[,c(3)])
  colnames(KEGGPathwayIDs)<-"pathwayName"
  return (KEGGPathwayIDs)
}

#' Calculate p-value and adjust p-value of KEGG pathway on input sample
#'
#' @param thispathwayID_DF KEGG pathway ID
#' @param thisKEGGSpecies species in KEGG format
#' @param thisInputSampleKEGG input sample with in KEGG gene ID
#' @param N the number of all genes with KEGG gene ID
#' @param n the number of sample genes with KEGG gene ID
#' @return vector in dataframe class containing p-value, adjust p-value,
#' N:the number of all genes with KEGG gene ID,
#' n the number of sample genes with KEGG gene ID,
#' M:the number of all genes in specific KEGG pathway,
#' m: the number of sample genes in specific KEGG pathway
#' @export
#' @examples
#' library(Key2Enrich)
#' data(inputSample)
#' inputSample<-as.data.frame(inputSample)
#' allGeneInPathwayDF<-getAllGeneInPathwayDF("mmu")
#' N<-getN(allGeneInPathwayDF)
#' n<-get_n(inputSample,allGeneInPathwayDF)
#' getPValue("mmu00053","mmu",inputSample,N,n)
getPValue <- function(thispathwayID_DF,thisKEGGSpecies,thisInputSampleKEGG,N,n) {
  thisSymbol<-getGeneSymbol(thisKEGGSpecies)
  #print (paste("Working on ",thispathwayID_DF[1],sep=""))
  #cat("\nCaculate M\n")
  geneInOnePathway<- keggLink(thisKEGGSpecies,thispathwayID_DF[1]) #keggLink("mmu","mmu05140")
  geneInOnePathwayDF<-data.frame(names(geneInOnePathway),geneInOnePathway)
  colnames(geneInOnePathwayDF)<-c("pathwayID","KEGG_GeneID")
  #M
  M<-length(unique(geneInOnePathwayDF$KEGG_GeneID))
  #print (M)

  #cat("Caculate m\n")
  inputSampleKEGG_M<-merge(thisInputSampleKEGG,
                           geneInOnePathwayDF,
                           by.x="KEGG_GeneID",
                           by.y="KEGG_GeneID")
  inputSampleKEGG_m<-unique(inputSampleKEGG_M$KEGG_GeneID)
  m<-length(inputSampleKEGG_m)
  #print (m)
  #print (inputSampleKEGG_m)

  #cat("N:\n")
  #print (N)
  #cat("n:\n")
  #print (n)

  #cat("Caculate p value\n")
  p.value <-  stats::phyper(q=m-1, m=M, n=(N-M), k=n, lower.tail=FALSE)
  #print (p.value)
  #print("###################################################")

 genes<- unique(as.character(inputSampleKEGG_M[,which(colnames(inputSampleKEGG_M)== thisSymbol)]))
 KEGGGenes=unique(as.character(inputSampleKEGG_M$KEGG_GeneID))
 entrezGenes=unique(as.character(inputSampleKEGG_M$entrezgene))
 ensemblGenes=unique(as.character(as.character(inputSampleKEGG_M$ensembl_gene_id)))

  return (data.frame(pathName=thispathwayID_DF[1],
  M=M,m=m,N=N,n=n,
  pValue=p.value,
  KEGGGenes=paste("\"",KEGGGenes,"\"",collapse=", ",sep=""),
  entrezGenes=paste("\"",entrezGenes,"\"",collapse=", ",sep=""),
  genes=paste("\"",genes,"\"",collapse=", ",sep=""),
  ensemblGenes=paste("\"",ensemblGenes,"\"",collapse=", ",sep="")))
}

#' Get all KEGG pathway ID and its corresponding gene list in KEGG gene ID format
#'
#' @param KEGGSpecies species in KEGG format
#' @return KEGG pathway ID, and its gene list in KEGG Gene ID format
#' @export
#' @examples
#' library(Key2Enrich)
#' getAllGeneInPathwayDF("mmu")

getAllGeneInPathwayDF<-function(KEGGSpecies){
  allGeneInPathway<-keggLink(KEGGSpecies, "pathway")
  allGeneInPathwayDF<-data.frame(names(allGeneInPathway),allGeneInPathway)
  colnames(allGeneInPathwayDF)<-c("pathID","KEGG_GeneID")
  return (allGeneInPathwayDF)
}

#' Get the number of all genes with KEGG gene ID
#'
#' @param allGeneInPathwayDF all KEGG pathway ID and its gene list in KEGG gene ID format
#' @return N:the number of all genes with KEGG gene ID
#' @export
#' @examples
#' library(Key2Enrich)
#' allGeneInPathwayDF<-getAllGeneInPathwayDF("mmu")
#' print(N<-getN(allGeneInPathwayDF))
#'
getN<-function (allGeneInPathwayDF){
  N<-length(unique(allGeneInPathwayDF$KEGG_GeneID))
  return (N)
}

#' Get the number of sample genes with KEGG gene ID
#'
#' @param inputSampleKEGG input sample with in KEGG gene ID
#' @param allGeneInPathwayDF all KEGG pathway ID and its gene list in KEGG gene ID format
#' @return n: the number of sample genes with KEGG gene ID
#' @export
#' @examples
#' library(Key2Enrich)
#' data(inputSample)
#' inputSample<-as.data.frame(inputSample)
#' allGeneInPathwayDF<-getAllGeneInPathwayDF("mmu")
#' print(n<-get_n(inputSample,allGeneInPathwayDF))
#'
get_n<-function (inputSampleKEGG,allGeneInPathwayDF){
  inputSampleKEGG_N<-merge(inputSampleKEGG,allGeneInPathwayDF,by.x="KEGG_GeneID",by.y="KEGG_GeneID")
  nrow(inputSampleKEGG_N)
  n<-length(unique(inputSampleKEGG_N$KEGG_GeneID))
  return (n)
}

#' Get gene list in specific KEGG pathway
#'
#' @param thispathID KEGG pathway ID
#' @param KEGGSpecies species in KEGG format
#' @return gene list in specific KEGG pathway
#' @export
#' @examples
#' library(Key2Enrich)
#' getGeneInOnePathwayDF("mmu05160","mmu")

getGeneInOnePathwayDF<-function(thispathID,KEGGSpecies){ #("mmu05160","mmu")
  geneInOnePathway<-keggLink(KEGGSpecies,thispathID) #keggLink("mmu","mmu05160")
  geneInOnePathwayDF<-data.frame(names(geneInOnePathway),geneInOnePathway)
  colnames(geneInOnePathwayDF)<-c("pathID","KEGG_GeneID")
  return (geneInOnePathwayDF)
}

#' Get sample genes in specific KEGG pathway
#'
#' @param thisInputSampleKEGG input sample with in KEGG gene ID
#' @param geneInOnePathwayDF gene list in specific KEGG pathway
#' @return sample genes in specific KEGG pathway
#' @export
#' @examples
#' library(Key2Enrich)
#' geneInOnePathwayDF<-getGeneInOnePathwayDF("mmu05160","mmu")
#' data(inputSample)
#' inputSample<-as.data.frame(inputSample)
#' get_inputSampleKEGG_m(inputSample,geneInOnePathwayDF)
#'
get_inputSampleKEGG_m<-function (thisInputSampleKEGG,geneInOnePathwayDF){
  inputSampleKEGG_M<-merge(thisInputSampleKEGG,
                           geneInOnePathwayDF,
                           by.x="KEGG_GeneID",
                           by.y="KEGG_GeneID")
  inputSampleKEGG_m<-unique(inputSampleKEGG_M$KEGG_GeneID)
  return (inputSampleKEGG_M[inputSampleKEGG_M$KEGG_GeneID %in% inputSampleKEGG_m,])
}

#' Plot heat map rectangle on the KEGG pathway node
#'
#' @param filterPValueName p-value and adjust-pvalue matrix of sample on all KEGG pathways
#' @param inputSampleKEGG input sample with in KEGG gene ID
#' @param inputSpecies KEGGSpecies species in KEGG format
#' @return NA
#' @import graphics
#' @import grDevices
#' @import utils
#' @import png

plotSigImg<-function(filterPValueName,inputSampleKEGG,inputSpecies){

  plotPathID<-as.character(filterPValueName[1])
  #plotPathID<-"mmu05160"
  #print (plotPathID)

  KGMLName<-paste(plotPathID,"xml",sep='.')
  path<-paste("http://rest.kegg.jp/get/","/kgml",sep=plotPathID)
  download.file(path,KGMLName,mode="wb",quiet=TRUE)
  #print ("downloading")

  imgName<-paste(plotPathID,"png",sep='.')
  path<-paste("http://rest.kegg.jp/get/","/image",sep=plotPathID)
  download.file(path,imgName,mode="wb",quiet=TRUE)

  r<-parseKGMLFile(KGMLName)
  kegg.nodes <- lapply(r[childIsEntry(r)], parseEntry)
  nodeDF<-data.frame(KEGG_GeneID="NA",type=NA,
                     link=NA,graphicName=NA,
                     graphicType=NA,graphicX=NA,
                     graphicY=NA,graphicWidth=NA,
                     graphicHeight=NA)
  nodeDF<-parseList2Dataframe(kegg.nodes,nodeDF)

  img <- readPNG(imgName)

  width<-dim(img)[2]
  height<-dim(img)[1]
  KEGGSpecies<-speciesKEGGConvert(inputSpecies)
  geneInOnePathwayDF<-getGeneInOnePathwayDF(plotPathID,KEGGSpecies)
  inputSampleKEGG_m<-get_inputSampleKEGG_m(inputSampleKEGG,geneInOnePathwayDF)
  node2Plot<-merge(inputSampleKEGG_m,nodeDF,by.x="KEGG_GeneID",by.y="KEGG_GeneID")
  node2Plot<-node2Plot[!(duplicated(node2Plot$graphicX) & duplicated(node2Plot$graphicY)),]

  rgbImg<-rgb(img[,,1], img[,,2], img[,,3])
  dim(rgbImg) <- dim(img)[1:2]

  #print (node2Plot)
  imgNamePlot<- paste("K2E", paste(plotPathID,"png",sep='.'),sep='_')

  png(imgNamePlot, width = width, height = height, res=300)
  op=par(mar = c(0, 0, 0, 0))
  plot(c(0, width), c(0, height), type = "n", xlab = "", ylab = "",xaxs = "i",yaxs = "i")
  rasterImage(img, 0, 0, width, height, interpolate = F)

  for (ii in 1:nrow(node2Plot)){

    if(!is.na(node2Plot$graphicWidth[ii])){
    geneWidth<-as.numeric(as.character(node2Plot$graphicWidth[ii]))
    geneHeight<-as.numeric(as.character(node2Plot$graphicHeight[ii]))
    geneX<-as.numeric(as.character(node2Plot$graphicX[ii]))
    geneY<-as.numeric(as.character(node2Plot$graphicY[ii]))

    rectX<- geneX-0.5*geneWidth+1
    rectX2<-geneX-0.5*geneWidth+1+geneWidth
    rectY<-geneY-0.5*geneHeight+1
    rectY2<-geneY-0.5*geneHeight+1+geneHeight

    # print (paste("rectX",rectX))
    # print (paste("rectX2",rectX2))
    # print (paste("rectY",rectY))
    # print (paste("rectY2",rectY2))

    plotArea<-rgbImg[rectY:rectY2,rectX:rectX2]

    thisLogFC<-as.numeric(node2Plot$log2FC[ii])

    plotArea[plotArea=="#BFFFBF"]<-getColor(as.numeric(thisLogFC))

    rasterImage(plotArea, geneX-0.5*geneWidth,
                height-geneY-0.5*geneHeight,
                geneX-0.5*geneWidth+geneWidth,
                height-geneY-0.5*geneHeight+geneHeight,
                interpolate=FALSE)
    }
  }
  plotGradient(width,height)
  textStamp(width)
  dev.off()
}
