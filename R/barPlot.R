#' Plot barplot
#'
#' @param filterValuePathNameDF  matrix of enriched dataset
#' @param type the type of mapping database, KEGG, Reactome, BP, MF, CC
#' @param imgWidth the width of export file
#' @param imgHeight the height of export file
#' @return data in dataframe class with pValue, adjust pValue
#' @export
#' @examples thisKEGGBarplot<-Key2EnrichBarplot(filterValuePathNameDF,"KEGG",15,20)

Key2EnrichBarplot<-function (filterValuePathNameDF,type,imgWidth,imgHeight)
{
  thisColor="black"
  fileName="pathway"
  if (type=="KEGG") {
    thisColor="cyan4"
    fileName<- paste(type,"pathway enrichment.pdf")
  }
  else if (type=="Reactome") {
    thisColor="springgreen4"
    fileName<- paste(type,"pathway enrichment.pdf")
    }
  else if (type=="BP") {
    thisColor="deeppink4"
    fileName<- "GO BP enrichment.pdf"
    }
  else if (type=="MF") {
    thisColor="navyblue"
    fileName<- "GO MF enrichment.pdf"
  }
  else if (type=="CC") {
    thisColor="tan4"
    fileName<- "GO CC enrichment.pdf"
  }

  pdf(file=fileName,width=imgWidth,height=imgHeight)

    filterValuePathNameDF$ratio<- as.numeric(as.character( filterValuePathNameDF$m ))/as.numeric(as.character( filterValuePathNameDF$M))

    p1<-ggplot(data= filterValuePathNameDF, aes(x=pathName, y=negativeLog)) +
      geom_bar(stat="identity",color=thisColor, fill=thisColor,width=0.5)+
      scale_x_discrete(limits= filterValuePathNameDF$pathName) +
      coord_flip()+ geom_hline(yintercept=1.3, linetype="dashed", color = "red")+
      labs(title="",x ="", y = "-log(q-value)")+
      theme_bw()+theme(axis.text=element_text(size=14),axis.title.x=element_text(size=14))

    p2<-ggplot(data= filterValuePathNameDF, aes(x=pathName, y=ratio)) +
      geom_line(data= filterValuePathNameDF, aes(x=pathName,y = filterValuePathNameDF$ratio,color="orange",group=1),show.legend = F,size=1)+
      geom_point(data= filterValuePathNameDF, mapping = aes(x=pathName,y = filterValuePathNameDF$ratio,group=1),size =2, shape = 21, fill = "yellow") +
      scale_x_discrete(limits= filterValuePathNameDF$pathName) +
      labs(title="",x ="", y = "Gene ratio")+
      coord_flip()+ theme_bw() %+replace%
      theme(panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),axis.title.y=element_blank(),
            axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.ticks.length=unit(0.18, "cm"))
    multiplot(p1,p2)

  dev.off()

}

#' Plot two image on same panel
#'
#' @param p1 ggplot
#' @param p2 ggplot
#' @return plot
#' @examples multiplot(p1,p2)

multiplot <- function(p1,p2) {
  grid.newpage()

  #extract gtable
  g1<-ggplot_gtable(ggplot_build(p1))
  g2<-ggplot_gtable(ggplot_build(p2))

  #overlap the panel of the 2nd plot on that of the 1st plot

  pp<-c(subset(g1$layout, name=="panel", se=t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]], pp$t, pp$l, pp$b,
                       pp$l)

  ## steal axis from second plot and modify
  ia <- which(g2$layout$name == "axis-b")
  ga <- g2$grobs[[ia]]
  ax <- ga$children[[2]]

  ## switch position of ticks and labels
  ax$heights <- rev(ax$heights)
  ax$grobs <- rev(ax$grobs)
  ax$grobs[[2]]$y <- ax$grobs[[2]]$y - unit(1, "npc") + unit(0.15, "cm")

  ## modify existing row to be tall enough for axis
  g$heights[[2]] <- g$heights[g2$layout[ia,]$t]

  ## add new axis
  g <- gtable_add_grob(g, ax, 2, 4, 2, 4)

  ## add new row for upper axis label
  g <- gtable_add_rows(g, g2$heights[1], 1)

  # draw it
  grid.draw(g)

}
