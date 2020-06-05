
### arguments ###
args=commandArgs(trailingOnly = TRUE)

runName=args[1]
totalMolFileNameExt=args[2]
sampleStatsFileNameExt=args[3]
ratioFileNameExt=args[4]
methMolFileNameExt=args[5]
sampleInfoFileNameExt=args[6]
fragInfoFileName=args[7]
outputDir=args[8]
directoriesFileName=args[9]
rHelperFxnsPath=args[10]

### run ###

setwd(outputDir)

source(rHelperFxnsPath)

directories=read.table(directoriesFileName)
directories=as.character(directories$V1)

aggregatedTotalMol=''
aggregatedSampleStats=''
aggregatedRatios=''
aggregatedMethMol=''
aggregatedSampleInfo=''

for(i in seq(1,length(directories),1)){
  currDir=directories[i]
  
  setwd(currDir)
  
  totalMolFileName=list.files(pattern = paste("\\",totalMolFileNameExt,"$",sep = ""))[1]
  sampleStatsFileName=list.files(pattern = paste("\\",sampleStatsFileNameExt,"$",sep = ""))[1]
  ratioFileName=list.files(pattern = paste("\\",ratioFileNameExt,"$",sep = ""))[1]
  methMolFileName=list.files(pattern = paste("\\",methMolFileNameExt,"$",sep = ""))[1]
  
  sampleInfoFileName=list.files(pattern = paste("\\",sampleInfoFileNameExt,"$",sep = ""))
  
  if(length(sampleInfoFileName)<1){
    sampleInfoFileName="NOSAMPLEINFO"
  }
  else{
    sampleInfoFileName=sampleInfoFileName[1]
  }
  
  outputFiles=loadFiles(fragInfoFileName,sampleStatsFileName,totalMolFileName,methMolFileName,methRatioFileName,sampleInfoFileName)

  if(i==1){
    aggregatedMethMol=outputFiles$methMolFile
    aggregatedRatios=outputFiles$methRatioFile
    aggregatedTotalMol=outputFiles$totalMolFile
    
    aggregatedSampleStats=outputFiles$sampleStats
  }
  else{
    aggregatedMethMol=mergeDF(aggregatedMethMol,outputFiles$methMolFile)
    aggregatedRatios=mergeDF(aggregatedRatios,outputFiles$methRatioFile)
    aggregatedTotalMol=mergeDF(aggregatedTotalMol,outputFiles$totalMolFile)
    
    aggregatedSampleStats=rbind(aggregatedSampleStats,outputFiles$sampleStats)
  }
}

setwd(outputDir)

plotData(runName,fragInfoFileName,aggregatedSampleStats,aggregatedTotalMol,aggregatedMethMol,aggregatedRatios)
  
