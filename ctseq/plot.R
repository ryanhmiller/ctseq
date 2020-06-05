
### arguments ###
args=commandArgs(trailingOnly = TRUE)

rHelperFxnsPath=args[1]
workingDir=args[2]
runName=args[3]
totalMolFileName=args[4]
methMolFileName=args[5]
ratioFileName=args[6]
sampleStatsFileName=args[7]
fragInfoFileName=args[8]
sampleInfoFileName=args[9]

### run ###

setwd(workingDir)

source(rHelperFxnsPath)

outputFiles=loadFiles(fragInfoFileName,sampleStatsFileName,totalMolFileName,methMolFileName,methRatioFileName,sampleInfoFileName)

plotData(runName,fragInfoFileName,outputFiles$sampleStats,outputFiles$totalMolFile,outputFiles$methMolFile,outputFiles$methRatioFile)
