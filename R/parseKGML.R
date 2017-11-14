#' Read XML file
#'
#' @param file file name
#' @return xmlRoot of file
#' @export
#' @examples
#' path<-system.file(package = "Key2Enrich")
#' xmlPath<-paste(path,"/extdata",sep="")
#' xmlFile<-paste(xmlPath,sep="","/hsa04012.xml")
#' parseKGMLFile(xmlFile)
#' @importFrom XML xmlTreeParse xmlAttrs xmlChildren xmlRoot xmlName xmlErrorCumulator

parseKGMLFile <- function(file) {
    tryCatch(
        doc <-  xmlTreeParse(file, getDTD=FALSE,
                            error=XML::xmlErrorCumulator(immediate=FALSE)),
        error = function(e) {
            fileSize <- file.info(file)$size[1]
            bytes <- sprintf("%d byte%s",
                             fileSize, ifelse(fileSize>1, "s", ""))
            msg <- paste("The file",
                         file,
                         "seems not to be a valid KGML file\n")
            if(fileSize<100L)
                msg <- paste(msg,
                             "[WARNING] File size (",
                             bytes,
                             ") is unsually small; please check.\n", sep="")
            msg <- paste(msg,
                         "\nDetailed error messages from",
                         "XML::xmlTreeParse:\n", sep="")
           cat(msg)
           stop(e)
        })
  r <- xmlRoot(doc)
  return (r)
  }

#' Get KEGG pathway info from xmlRoot
#'
#' @param r xmlRoot
#' @return list class with entryID,name,type,link, X of graphic,Y of graphic,width of graphic,height of graphic
#' @export
#' @examples #getPathInfo(r)
  getPathInfo<-function(r){
  ## parse them
    attrs <- xmlAttrs(r)
  ## required: name, org, number
  pathName <- attrs[["name"]]
  pathOrg <- attrs[["org"]]
  pathNumber <- attrs[["number"]]
  ## implied: title, image, link
  pathTitle <- attrs[["title"]]
  pathImage <-   attrs[["image"]]
  pathLnk <-    attrs[["link"]]
  return (data.frame(pathName=pathName,
                     pathOrg=pathOrg,
                     pathNumber=pathNumber,
                     pathTitle=pathTitle,
                     pathImage=pathImage,
                     pathLnk=pathLnk))
  }

  childIsEntry<-function(r){
  ## possible elements: entry, relation and reaction
  childnames <- sapply(xmlChildren(r), xmlName)
  isEntry <- childnames == "entry"
  return (isEntry)
  }

  parseEntry <- function(entry) {
  attrs <- xmlAttrs(entry)

  entryID <- attrs[["id"]]
  name <- unname(unlist(strsplit(attrs["name"]," ")))
  type <- attrs[["type"]]
  #link <- attrs[["link"]]
  link<-""

  graphics <- xmlChildren(entry)$graphics
  graphicsAttrs <- xmlAttrs(graphics)

  graphicName<-unname(graphicsAttrs["name"])
  graphicType<-unname(graphicsAttrs["type"])
  graphicX<-unname(graphicsAttrs["x"])
  graphicY<-unname(graphicsAttrs["y"])
  graphicWidth<-unname(graphicsAttrs["width"])
  graphicHeight<-unname(graphicsAttrs["height"])

  return (list(entryID=entryID,
               name=name,
               type=type,
               link=link,
  graphicName=graphicName,
  graphicType=graphicType,
  graphicX=graphicX,
  graphicY=graphicY,
  graphicWidth=graphicWidth,
  graphicHeight=graphicHeight))
  }

  #' Get info from list class of graphic info
  #'
  #' @param thisGeneKEGGID gene KEGG ID
  #' @return character vector with name,type, X of graphic,Y of graphic,width of graphic,height of graphic
  #' @export
  #' @examples #parseGraphic(thisGeneKEGGID)
   parseGraphic <- function(thisGeneKEGGID) {
   graphicName<-unname(thisGeneKEGGID["name"])
   graphicType<-unname(thisGeneKEGGID["type"])
   graphicX<-unname(thisGeneKEGGID["x"])
   graphicY<-unname(thisGeneKEGGID["y"])
   graphicWidth<-unname(thisGeneKEGGID["width"])
   graphicHeight<-unname(thisGeneKEGGID["height"])

   return (c(graphicName,graphicType,graphicX,graphicY,graphicWidth,graphicHeight))
}

   #' Convert list class of node info to dataframe
   #'
   #' @param nodeList list class of node info
   #' @param nodeDF output dataframe
   #' @export
   #' @return dataframe class of node info
   #' @examples #parseList2Dataframe(nodeList,nodeDF)
parseList2Dataframe<-function(nodeList,nodeDF){


	for(i in 1:length(nodeList)){
	node<-nodeList[[i]]
	if(node$type=="gene"){
	thisNodeGene<-node$name
	thisNodeType<-node$type
	thisNodeLink<-node$link
	thisNodeGraphicName<-node$graphicName
	thisNodeGraphicType<-node$graphicType
	thisNodeGraphicX<-node$graphicX
	thisNodeGraphicY<-node$graphicY
	thisNodeGraphicWidth<-node$graphicWidth
	thisNodeGraphicHeight<-node$graphicHeight

	thisNodeGeneNum<-length(thisNodeGene)

   for (m in 1:thisNodeGeneNum)
   {
   newNodeDF<-data.frame(KEGG_GeneID=thisNodeGene[m],
   type=thisNodeType,
   link=thisNodeLink,
	graphicName=thisNodeGraphicName,
	graphicType=thisNodeGraphicType,
	graphicX=thisNodeGraphicX,
	graphicY=thisNodeGraphicY,
	graphicWidth=thisNodeGraphicWidth,
	graphicHeight=thisNodeGraphicHeight)
	nodeDF<-rbind(newNodeDF,nodeDF)
   }
}
}
nodeDF<-nodeDF[-c(nrow(nodeDF)),]
return (nodeDF)
}
