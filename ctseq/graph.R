
args=commandArgs(trailingOnly = TRUE)

workingDir=args[1]
totalMoleculeExtension=args[2]
fragOrderFileName=args[3]

setwd(workingDir)


library(ggplot2)
library(reshape)

moleculeFileName=list.files(pattern = paste("\\",totalMoleculeExtension,"$",sep=""))[1]

moleculeFile=read.table(moleculeFileName, header = T, row.names = 1, sep = "\t")

fragOrder=read.table(fragOrderFileName)
fragOrder=as.character(fragOrder$V1)


moleculeFileName=list.files(pattern = paste("\\",totalMoleculeExtension,"$",sep=""))[1]
runName=strsplit(moleculeFileName,split = "_")[[1]][1]

moleculeFile=read.table(moleculeFileName, header = T, sep = "\t")

meltedMol=melt(moleculeFile)
colnames(meltedMol)=c("Fragment","Sample","Molecules")

meltedMol[is.na(meltedMol)]=1


pdf(paste(runName,"_totalMolecules.pdf",sep=""))
ggplot(meltedMol, aes(x=factor(Fragment), y=Molecules, fill=factor(Sample))) +
  ggtitle(paste(runName,"unique molecules")) +
  xlab("Fragment") +
  scale_fill_discrete(name = "Sample") +
  geom_dotplot(binaxis = "y", dotsize=0.4) +
  scale_x_discrete(limits=fragOrder) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 0.01, color = "white")) +
  scale_y_log10(limits = c(1,10000))
dev.off()
