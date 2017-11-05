#' Get color for node from log fold change value
#'
#' @param thisValue log fold change value of input
#' @return color for node
#' @examples getColor(1.6)

getColor<-function (thisValue){
  thisCOlor<-"#BEBEBE"
  thisArray<-seq(0,2,length.out = 10)
  colfunc <- colorRampPalette(c("white", "red"))
  colArray<-colfunc(10)

  thisArray2<-seq(-2,0,length.out = 10)
  colfunc2 <- colorRampPalette(c("blueviolet", "white"))
  colArray2<-colfunc2(10)

  if (thisValue>=2) {
    thisCOlor<-colArray[10]
  }
  else if (thisValue>0 & thisValue<2){
    thisCOlor<-colArray[getIndex(thisArray,thisValue)]
  }
  else if (thisValue<=-2) {
    thisCOlor<-colArray2[1]
  }
  else if (thisValue<0 & thisValue>-2) {
    thisCOlor<-colArray2[getIndex2(thisArray2,thisValue)]
  }else if (thisValue==0){
    thisCOlor<-"#BEBEBE"
  }

  return (thisCOlor)
}

#' Get index of value in array for right color key
#'
#' @param thisArray color array
#' @param thisValue log fold change value
#' @return index
#' @examples  thisArray<-seq(0,2,length.out = 10)
#' @examples  getIndex(thisArray,1.6)

getIndex<-function (thisArray,thisValue){
  thisIndex<-1
  l<-length(thisArray)-1
  for (i in 1:l){
    if(thisValue>=thisArray[i] & thisValue<thisArray[i+1]) {
      thisIndex<-i
      break
    }
  }
  return (thisIndex)
}

#' Get index of value in array for left color key
#'
#' @param thisArray color array
#' @param thisValue log fold change value
#' @return index
#' @examples  thisArray2<-seq(-2,0,length.out = 10)
#' @examples  getIndex2(thisArray2,1.6)

getIndex2<-function (thisArray,thisValue){
  thisIndex<-10
  l<-length(thisArray)-1
  for (i in 1:l){
    if(thisValue>thisArray[i] & thisValue<=thisArray[i+1]) {
      thisIndex<-i+1
      break
    }
  }
  return (thisIndex)
}

#' Plot heatmap color key for node
#'
#' @param thisWidth the width of color key
#' @param thisHeight the height of color key
#' @return NA
#' @examples  plotGradient(120,5)
#'
plotGradient<- function (thisWidth,thisHeight){
  gradientWidth<-7
  gradientHeight<-20
  x<-6
  y<-20

  colfunc <- colorRampPalette(c("white", "red"))
  colArray<-colfunc(10)

  colfunc2 <- colorRampPalette(c("blueviolet", "white"))
  colArray2<-colfunc2(10)
  for (i in 10:1){
  rect(thisWidth-y-gradientWidth*i,thisHeight-x-gradientHeight,thisWidth-y-(i-1),thisHeight-x,col=colArray[11-i],border = "black",lwd=0.2)
  rect(thisWidth-y-gradientWidth*i-gradientWidth*10,thisHeight-x-gradientHeight,thisWidth-y-(i-1)-gradientWidth*10,thisHeight-x,col=colArray2[11-i],border = "black",lwd=0.2)
  }
  text(thisWidth-y,thisHeight-x-gradientHeight-5,"1",cex=0.2)
  text(thisWidth-y-gradientWidth*20,thisHeight-x-gradientHeight-5,"-1",cex=0.2)
}

#' Add text stamp on image
#'
#' @param thisWidth the width of text stamp
#' @return NA
#' @examples  textStamp(100)
#'
textStamp<- function (thisWidth){
  text(thisWidth-136,11,"Rendered by Key2Enrich on KEGG",cex=0.3,font=2)
}
