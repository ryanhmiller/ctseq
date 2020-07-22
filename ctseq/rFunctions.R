# make sure a list of names is in the correct format
formatNames=function(listOfNames){
  newListOfNames=c()
  for(i in seq(1,length(listOfNames),1)){
    currName=listOfNames[i]
    firstChar=substr(currName,1,1)
    
    if(firstChar=="0" || firstChar=="1" || firstChar=="2" || firstChar=="3" || firstChar=="4" || firstChar=="5" || firstChar=="6" || firstChar=="7" || firstChar=="8" || firstChar=="9"){
      currName=paste("X",currName,sep = "")
    }
    
    currName=gsub('-','.',currName) # replace '-' with '.'
    newListOfNames=append(newListOfNames,currName)
  }
  return(newListOfNames)
}


# merge two data frames
mergeDF=function(df1,df2){
  finalDF=merge(df1,df2,by="row.names",all=T)
  rownames(finalDF)=finalDF$Row.names
  finalDF$Row.names=NULL
  return(finalDF)
}


# load and assemble all the data necessary to plot
loadFiles=function(fragInfoFileName,sampleStatsFileName,totalMolFileName,methMolFileName,methRatioFileName,sampleInfoFileName){
  # info about fragments
  fragInfo=read.table(fragInfoFileName, header = T, row.names = 1, sep = '\t')
  
  # get order of frags for molecule depth plot
  fragOrder_molDepthPlot=fragInfo[with(fragInfo, order(fragOrder)),]
  fragOrder_molDepthPlot=row.names(fragOrder_molDepthPlot)
  
  fragOrder_heatmaps=rownames(fragInfo)
  
  # info about each sample
  sampleStats=read.table(file = sampleStatsFileName, header = T, row.names = 1, sep = "\t")
  
  
  # R is weird adds "X" to column names beginning with an integer
  # rename sample names in 'sampleStats' if first character of sample name is an integer
  # replace '-' with '.'
  newSampleStatsRowNames=formatNames(rownames(sampleStats))
  rownames(sampleStats)=newSampleStatsRowNames
  
  ##
  
  totalMolFile=read.table(totalMolFileName, header = T, row.names = 1, sep = "\t")
  totalMolFile=totalMolFile[fragOrder_heatmaps,] # only grab frags in info file - we only want to plot these
  
  methMolFile=read.table(methMolFileName,header = T, row.names = 1, sep = "\t")
  methMolFile=methMolFile[fragOrder_heatmaps,]
  
  methRatioFile=read.table(ratioFileName,header = T, row.names = 1, sep = "\t")
  methRatioFile=methRatioFile[fragOrder_heatmaps,]
  

  if(sampleInfoFileName!="NOSAMPLEINFO"){
    sampleInfoFile=read.table(sampleInfoFileName, header = T, row.names = 1, sep = "\t")

    newSampleInfoRowNames=formatNames(rownames(sampleInfoFile))
    rownames(sampleInfoFile)=newSampleInfoRowNames

    # put data file cols in same order at sampleInfoFile rows
    totalMolFile=totalMolFile[,rownames(sampleInfoFile)]
    colnames(totalMolFile)=paste(rownames(sampleInfoFile),as.character(sampleInfoFile[,1]),sep = '.')
    colnames(totalMolFile)=formatNames(colnames(totalMolFile))

    methMolFile=methMolFile[,rownames(sampleInfoFile)]
    colnames(methMolFile)=paste(rownames(sampleInfoFile),as.character(sampleInfoFile[,1]),sep = '.')
    colnames(methMolFile)=formatNames(colnames(methMolFile))

    methRatioFile=methRatioFile[,rownames(sampleInfoFile)]
    colnames(methRatioFile)=paste(rownames(sampleInfoFile),as.character(sampleInfoFile[,1]),sep = '.')
    colnames(methRatioFile)=formatNames(colnames(methRatioFile))


    sampleStats=sampleStats[rownames(sampleInfoFile),]
    rownames(sampleStats)=paste(rownames(sampleInfoFile),as.character(sampleInfoFile[,1]),sep = '.')
    rownames(sampleStats)=formatNames(rownames(sampleStats))

  }
  
  returnList=list(fragOrder_molDepthPlot=fragOrder_molDepthPlot,
                  totalMolFile=totalMolFile,
                  methMolFile=methMolFile,
                  methRatioFile=methRatioFile,
                  sampleStats=sampleStats)
  
  return(returnList)
}


# molecule depth plot and heatmaps
plotData=function(runName,fragInfoFileName,sampleStats,totalMolFile,methMolFile,methRatioFile){
  library(ggplot2)
  library(reshape)
  library(pheatmap)
  
  # info about fragments
  fragInfo=read.table(fragInfoFileName, header = T, row.names = 1, sep = '\t')
  
  # get order of frags for molecule depth plot
  fragOrder_molDepthPlot=fragInfo[with(fragInfo, order(fragOrder)),]
  fragOrder_molDepthPlot=row.names(fragOrder_molDepthPlot)
  
  fragOrder_heatmaps=rownames(fragInfo)
  
  ##### create color palette for all annotations we will use in the heatmaps #####
  annotationFeatures=c(colnames(fragInfo),colnames(sampleStats))
  
  # colors to pick from
  colorList=c("aquamarine4","blue4","brown4","burlywood4","cadetblue4","chartreuse4","chocolate4","coral4","darkgoldenrod4","darkgreen",
              "darkmagenta","darkorange3","darkorchid4","darkred","darkslateblue","darkslategray","deeppink4","deepskyblue4","dodgerblue4", "firebrick3",
              "gray9","indianred3","hotpink3","purple4","navy","olivedrab4","magenta3","sienna","wheat4","turquoise4",
              "palevioletred3","salmon","cornflowerblue","palegreen4","maroon","green4","royalblue3","rosybrown4","orange","purple")
  
  set.seed(1) # set seed to keep order of colors the same
  colorList=sample(colorList) # shuffle colors
  
  my_ann_colors = list()
  
  for(i in seq(1,length(annotationFeatures),1)){
    feature=annotationFeatures[i]
    color=colorList[i]
    
    my_ann_colors[[feature]]=c("grey95",color)
    
  }
  
  
  ####### total molecules plot #######
  # load total molecules file
  
  # totalMolFile=read.table(totalMolFileName, header = T, row.names = 1, sep = "\t")
  # totalMolFile=totalMolFile[fragOrder_heatmaps,] # only grab frags in info file - we only want to plot these
  colNamesTotalMol=colnames(totalMolFile)
  colNamesTotalMol=c("Locus",colNamesTotalMol)
  
  totalMolFile$Locus=rownames(totalMolFile)
  totalMolFile=totalMolFile[,colNamesTotalMol]
  rownames(totalMolFile)=NULL
  
  meltedMol=melt(totalMolFile) 
  colnames(meltedMol)=c("Fragment","Sample","Molecules")
  
  # change all NA's to 1 because will be taking log
  meltedMol[is.na(meltedMol)]=1
  
  totalMolPlotName=paste(runName,"_totalMoleculesPlot.pdf",sep="")
  
  moleculePlot=ggplot(meltedMol, aes(x=factor(Fragment), y=Molecules, fill=factor(Sample))) +
    ggtitle(paste(runName,"Unique Molecules")) +
    xlab("Fragment") +
    scale_fill_discrete(name = "Sample") +
    geom_dotplot(binaxis = "y", dotsize=0.4) +
    scale_x_discrete(limits=fragOrder_molDepthPlot) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 0.01, color = "white")) +
    scale_y_log10(limits = c(1,10000))
  
  pdf(totalMolPlotName)
  print(moleculePlot)
  dev.off()
  
  
  
  
  ####### HEATMAPS ########
  rowFontSize=3
  colFontSize=6
  fontSize=6
  
  if(ncol(methMolFile)>25 && ncol(methMolFile)<40){
    colFontSize=4
  }
  else if(ncol(methMolFile)>=40){
    colFontSize=3
  }
  
  if(nrow(fragInfo)>160){
    rowFontSize=2
  }
  
  ### total molecule heatmap ###
  rownames(totalMolFile)=totalMolFile$Locus
  totalMolFile$Locus=NULL
  
  fragInfo$fragOrder=NULL
  
  totalMolHeatmapName=paste(runName,"_totalMoleculesHeatmap.pdf",sep = "")
  
  pdf(totalMolHeatmapName)
  pheatmap(log2(totalMolFile[fragOrder_heatmaps,]+1), main=paste(runName,"Aligned Molecules (log2)"),  fontsize_row=rowFontSize, fontsize_col=colFontSize, fontsize = fontSize,
           show_rownames = T, cluster_cols=F, cluster_rows=F, color=colorRampPalette(c("navy","yellow","forestgreen","darkgreen"))(50),
           annotation_col = sampleStats, annotation_row = fragInfo, na_col = "snow4", annotation_colors = my_ann_colors)
  dev.off()
  
  ### meth mol heatmap ###
  # methMolFile=read.table(methMolFileName,header = T, row.names = 1, sep = "\t")
  
  methMolHeatmapName=paste(runName,"_methylatedMoleculesHeatmap.pdf",sep = "")
  
  pdf(methMolHeatmapName)
  pheatmap(log2(methMolFile[fragOrder_heatmaps,]+1), main=paste(runName,"Methylated Molecules (log2)"), fontsize_row=rowFontSize, fontsize_col=colFontSize, fontsize = fontSize,
           show_rownames = T, cluster_cols=F, cluster_rows=F, color=colorRampPalette(c("navy","white","firebrick3","firebrick4"))(50),
           annotation_col = sampleStats, annotation_row = fragInfo, na_col = "snow4", annotation_colors = my_ann_colors)
  dev.off()
  
  
  ### meth ratio heatmap ###
  # methRatioFile=read.table(ratioFileName,header = T, row.names = 1, sep = "\t")
  
  methRatioHeatmapName=paste(runName,"_methylationRatioHeatmap.pdf",sep = "")
  
  pdf(methRatioHeatmapName)
  pheatmap(methRatioFile[fragOrder_heatmaps,], main=paste(runName,"Methylation Ratio"), fontsize_row=rowFontSize, fontsize_col=colFontSize, fontsize = fontSize,
           show_rownames = T, cluster_cols=F, cluster_rows=F, color=colorRampPalette(c("navy","white","sienna1","sienna2","sienna3","sienna","sienna4"))(50),
           annotation_col = sampleStats, annotation_row = fragInfo, na_col = "snow4", annotation_colors = my_ann_colors)
  dev.off()
  
  
}