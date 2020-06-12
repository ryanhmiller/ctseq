from multiprocessing import Pool
import sys
import re
from simplesam import Reader, Writer
import os
from functools import partial
from umi_tools import UMIClusterer
import gc
from . import utilities

### I've had to delete some objects after use because the memory was ballooning on very large data sets

##### LOAD REFERENCE SEQUENCES INTO DICT #####
def getRefSeq(refFileName):
    refSeq={}
    refLociNames=[]
    with open(refFileName, "r") as refFile:
        lines=refFile.readlines()
        for i in range(0,len(lines),2):
            fragName=lines[i].strip("\n")[1:]
            seq=lines[i+1].strip("\n")
            refSeq[fragName]=seq.upper()

            refLociNames.append(fragName)

    return(refSeq,refLociNames)


def getRefCGindices(refSeqDict):
    refCGindices={}
    for frag in refSeqDict:
        seq=refSeqDict[frag]
        cgIndices=[]

        # iterate through ref seq to find indices of CGs
        for i in range(0,(len(seq)-1),1):
            kmer=seq[i:(i+2)].upper()

            if kmer=="CG":
                cgIndices.append(i)

        refCGindices[frag]=cgIndices

    return(refCGindices)

def splitSamFile(samFileName,numProcesses):
    # find the number of header lines
    numHeaderLines=0
    headerLines=[]
    with open(samFileName,"r") as infile:
        # figure out how many header lines there are
        for line in infile:
            if line[0] == "@":
                numHeaderLines+=1
                headerLines.append(line)
            else:
                break

    intermedFileNames=[]
    with open(samFileName,"r") as infile:
        sampleName=samFileName.split('.')[0]

        print('reading in all reads for ',samFileName,utilities.getDate())
        allLines=infile.readlines()
        allLines=allLines[numHeaderLines:] # remove headerlines
        print('done reading in all reads',utilities.getDate())

        sizeChunks=len(allLines)//numProcesses
        # if odd, need to add one to the number to make sure it's even
        if sizeChunks%2 != 0:
            sizeChunks+=1

        myIndex=0

        print('splitting reads into '+str(numProcesses)+' temp files',utilities.getDate())
        for i in range(0,numProcesses,1):
            currIntermedFile=sampleName+"_TEMP_"+str(i)+".sam"
            intermedFileNames.append(currIntermedFile)

            with open(currIntermedFile,"w") as intermedFile:
                if i != (numProcesses-1):
                    intermedFile.write("".join(headerLines+allLines[myIndex:myIndex+sizeChunks]))
                    myIndex+=sizeChunks
                else:
                    intermedFile.write("".join(headerLines+allLines[myIndex:])) # the last temp file contains all the remaining lines

        # clear up RAM
        del allLines
        print("done splitting reads into temp files",utilities.getDate())

    return(intermedFileNames)


class Read:
  def __init__(self, myRead, myRefSeqDict, myRefCGindices):
    self.umi=myRead.qname.split("+")[1]
    self.locus=myRead.rname
    self.cigar=myRead.cigar
    alignedRef="".join(myRead.parse_md())
    self.alignedRef=alignedRef
    alignedCoords=myRead.coords
    self.startCoord=alignedCoords[0]-1 # WE WANT THE 0-BASED COORDINATES!!!
    self.endCoord=self.startCoord+len(alignedRef)-1

    self.refCGindices=myRefCGindices[self.locus]

    relativeRefPos=0
    seqPos=0

    self.gappedmethIndexString=myRead.gapped('XM')
    self.methIndexString=[]


    self.methIndices=[]
    self.unmethIndices=[]
    self.methString=""

    for realIndex in self.refCGindices:
        relativeIndex=realIndex-self.startCoord
        if realIndex >= self.startCoord and realIndex <= self.endCoord: #is ref CG index in this read?
            gappedChar=self.gappedmethIndexString[relativeIndex]
            if gappedChar == "Z":
                self.methIndices.append(realIndex)
            else:
                self.unmethIndices.append(realIndex)

            self.methIndexString.append(str(realIndex)+gappedChar)
            self.methString+=gappedChar

    self.methIndexString=":".join(self.methIndexString)



# class to act as a wrapper around results from analyzing each temp/intermediate file
class ReadsChunkResults:
    def __init__(self, myIntermedFileName, myReadDepthDict, myLocusUMIreadsDict):
        self.intermedFileName=myIntermedFileName
        self.readDepthDict=myReadDepthDict
        self.locusUMIreadsDict=myLocusUMIreadsDict


# function to analyze the reads from a temp/intermediate file
def analyzeReads(myIntermedFileName,myRefSeq,myRefCGindices):
    locusUMIreadsDict={} # locus : UMI : allReadsforthatUMI

    readDict={} # fragName:totalReads
    # initialize readDict
    for frag in myRefSeq:
        readDict[frag]=0


    infile = open(myIntermedFileName, 'r')
    samFile = Reader(infile)


    for forwardRead in samFile:
        if forwardRead.mapped and forwardRead.paired: # makes sure reads map to a reference seq and are paired
            reverseRead=samFile.next()

            #### UMI...
            if forwardRead.rname == reverseRead.rname: ### this is ADDED 11/20/19 to make sure both the forward and the reverse read map to the same fragment

                # make sure read name, UMI are same
                if forwardRead.qname != reverseRead.qname:
                    print("READS NOT PAIRED, out of sync")
                    print("forwardRead: ",forwardRead.qname,"reverseRead: ",reverseRead.qname)
                    print("EXITING...")
                    sys.exit()

                fUMI=forwardRead.qname.split("+")[1]
                rUMI=reverseRead.qname.split("+")[1]

                # if "N" not in fUMI and "N" not in rUMI: ### ADDED 11/21/19 - PREVENTS US FROM ANALYZING READS THAT HAVE "N"s IN THE UMIs

                forwardMethylation=Read(forwardRead, myRefSeq, myRefCGindices)

                reverseMethylation=Read(reverseRead, myRefSeq, myRefCGindices)

                myLocus=forwardMethylation.locus

                #ADDED - keeping track of total number of reads - PUT THIS AT THE END
                readDict[myLocus]+=1

                locusCGindices=myRefCGindices[myLocus]

                myUMI=forwardMethylation.umi

                consensusIndexString=""
                consensusMethString=""
                consensusMethIndices=[]
                consensusUnmethIndices=[]


                for index in locusCGindices:
                    if index >= forwardMethylation.startCoord and index < reverseMethylation.startCoord and index in forwardMethylation.methIndices:
                        consensusIndexString+=(str(index)+"Z")
                        consensusMethString+="Z"
                        consensusMethIndices.append(index)
                    elif index > forwardMethylation.endCoord and index <= reverseMethylation.endCoord and index in reverseMethylation.methIndices:
                        consensusIndexString+=(str(index)+"Z")
                        consensusMethString+="Z"
                        consensusMethIndices.append(index)
                    elif index >= forwardMethylation.startCoord and index >= reverseMethylation.startCoord and index <= forwardMethylation.endCoord and index <= reverseMethylation.endCoord and index in forwardMethylation.methIndices and index in reverseMethylation.methIndices: #index >= forwardMethylation.startCoord was in here twice. changed second one to index >= reverseMethylation.startCoord
                        consensusIndexString+=(str(index)+"Z")
                        consensusMethString+="Z"
                        consensusMethIndices.append(index)
                    else:
                        consensusIndexString+=(str(index)+"z")
                        consensusMethString+="z"
                        consensusUnmethIndices.append(index)

                if myLocus in locusUMIreadsDict:
                    if myUMI in locusUMIreadsDict[myLocus]:
                        locusUMIreadsDict[myLocus][myUMI].append(consensusMethString)
                    else:
                        locusUMIreadsDict[myLocus][myUMI]=[consensusMethString]
                else:
                    readList=[consensusMethString]
                    myUMIdict={myUMI:readList}
                    locusUMIreadsDict[myLocus]=myUMIdict

    myResults=ReadsChunkResults(myIntermedFileName,readDict,locusUMIreadsDict)
    return myResults

### combine results from first parallel analysis
# function to merge the locusUMIreadDicts
def mergeLocusUMIreadDicts(dict1,dict2):
    # merging things into dict1
    for locus in set(list(dict1.keys())+list(dict2.keys())):
        if locus in dict1 and locus in dict2:
            for UMI in set(list(dict1[locus].keys())+list(dict2[locus].keys())):
                if UMI in dict1[locus] and UMI in dict2[locus]:
                    dict1[locus][UMI]+=dict2[locus][UMI]

                elif UMI not in dict1[locus] and UMI in dict2[locus]:
                    dict1[locus][UMI]=dict2[locus][UMI]


        elif locus not in dict1 and locus in dict2:
            dict1[locus]=dict2[locus]

# function to merge the readDepthDict
def mergeReadDepthDicts(dict1,dict2):
    # merging things into dict1
    for locus in set(list(dict1.keys())+list(dict2.keys())):
        if locus in dict1 and locus in dict2:
            dict1[locus]+=dict2[locus]

        elif locus not in dict1 and locus in dict2:
            dict1[locus]=dict2[locus]


# split up combined locusUMIreadsDict into separate dicts to analyze in parallel
class MethylationChunkInput:
    def __init__(self, myLocusUMIreadsDict, myRefCGindices, myConsensusCutoff,mySampleName, myTempNumber, myOutputDir, myUMIcollapseThreshold, myUMIcollapseAlg):
        self.locusUMIreadsDict=myLocusUMIreadsDict
        self.refCGindices=myRefCGindices
        self.myConsensusCutoff=myConsensusCutoff
        self.mySampleName=mySampleName
        self.myTempNumber=myTempNumber
        self.myOutputDir=myOutputDir
        self.myUMIcollapseThreshold=myUMIcollapseThreshold
        self.myUMIcollapseAlg=myUMIcollapseAlg

#### counting molecules and methlyation concordance for a subset of loci
def callMolecules(myMethylationInput):
    mySampleName=myMethylationInput.mySampleName
    tempFileNumber=myMethylationInput.myTempNumber
    myOutputDir=myMethylationInput.myOutputDir
    collapseThreshold=myMethylationInput.myUMIcollapseThreshold
    collapseAlg=myMethylationInput.myUMIcollapseAlg


    tempFileName=myOutputDir+mySampleName+"_TEMP_ALLMOLECULES_"+str(tempFileNumber)+".txt"
    tempFile=open(tempFileName,"w")



    for locus in myMethylationInput.locusUMIreadsDict:

        #####################
        ### collapse UMIs ###
        #####################

        myLocusUMIs={}   # UMI (in byte form, b'ATGTC') : len(number reads with UMI)

        for UMI in myMethylationInput.locusUMIreadsDict[locus]:
            myLocusUMIs[UMI.encode()]=len(myMethylationInput.locusUMIreadsDict[locus][UMI])


        # initialize UMI clusterer
        clusterer = UMIClusterer(cluster_method=collapseAlg)
        clustered_umis = clusterer(myLocusUMIs, threshold=collapseThreshold)


        for cluster in clustered_umis:
            if len(cluster) > 1:
                headUMI="" # we'll just combine all the reads for all the UMIs in this cluster into the dictionary entry for the first UMI (called "headUMI" here)
                for i in range(1,len(cluster),1): # since combining everything into first UMI's entry, start on second entry
                    headUMI=cluster[0].decode()
                    currUMI=cluster[i]
                    currUMI=currUMI.decode()
                    myMethylationInput.locusUMIreadsDict[locus][headUMI]+=myMethylationInput.locusUMIreadsDict[locus][currUMI] # combine reads together into one UMI entry
                    del myMethylationInput.locusUMIreadsDict[locus][currUMI] # remove extra UMI entry

        # print("locus: ",locus)
        # print(len(myLocusUMIs))
        # print(len(clustered_umis))
        # print(len(list(myMethylationInput.locusUMIreadsDict[locus].keys())))
        # print(clustered_umis[:5],"\n")

        #################################
        #################################

        for UMI in myMethylationInput.locusUMIreadsDict[locus]:
            consensusString=""
            consensusMeth=0
            consensusUnmeth=0

            for i in range(0,len(myMethylationInput.refCGindices[locus]),1):
                methReads=0
                totalReads=0

                for copyTemplate in myMethylationInput.locusUMIreadsDict[locus][UMI]:
                    if copyTemplate[i]=="Z":
                        methReads+=1
                    totalReads+=1

                if float(methReads)/totalReads >= myMethylationInput.myConsensusCutoff:
                    consensusString+="Z"
                    consensusMeth+=1
                else:
                    consensusString+="z"
                    consensusUnmeth+=1

            consensusLength=consensusMeth+consensusUnmeth

            myOutputLine=[mySampleName,locus,UMI,consensusString,str(consensusMeth),str(consensusLength),str(consensusMeth/consensusLength),str(len(myMethylationInput.locusUMIreadsDict[locus][UMI]))]

            tempFile.write("\t".join(myOutputLine)+"\n")

    tempFile.close()
    del myMethylationInput # we need this otherwise the memory balloons with large files

def run(args):
    refDir=args.refDir
    samDir=args.dir
    consensus=float(args.consensus)
    processes=int(args.processes)
    umiCollapseThreshold=int(args.umiThreshold)
    umiCollapseAlg=args.umiCollapseAlg


    #############
    # arg check #
    #############
    # check that paths are valid
    refDir=utilities.validDir(refDir)
    samDir=utilities.validDir(samDir)

    # check that samDir/fastqDir actually has sam/reference files in them
    utilities.fileCheck(refDir,".fa")
    utilities.fileCheck(samDir,".sam")

    # make sure sam file dir has '/' on end
    if samDir[-1]!='/':
        samDir+='/'

    ####################
    ### SET UP STUFF ###
    ####################
    # get reference .fa file
    refFilePath=utilities.getReferenceFile(refDir)

    ### get reference sequences from reference file
    (refSeq,refLociNames)=getRefSeq(refFilePath)

    ### get indices of cgs in reference sequences
    refCGindices=getRefCGindices(refSeq)
    ####################


    os.chdir(samDir)
    listOfSamFiles=utilities.getFiles(samDir,'.sam')

    ## loop through sam files
    for samFileName in listOfSamFiles:
        # call_molecules(args,samFileName)
        if "TEMP" not in samFileName: # so don't start analyzing old temp files
            sampleName=samFileName.split(".")[0]
            print('\n**************')
            print('Calling molecules for',sampleName,utilities.getDate())
            print('**************\n')

            #### split sam file into smaller files (n=numProcesses)
            os.chdir(samDir)
            intermedFileNames=splitSamFile(samFileName,processes)
            # gc.collect() # free up unreferenced memory

            ### analyze the reads from the temp .sam files for this sample in parallel
            chunksAnalyzingReads=[]
            print('analyzing the reads for ',sampleName,'in parallel',utilities.getDate())
            p = Pool(processes)

            analyzeReads_multipleArgs=partial(analyzeReads, myRefSeq=refSeq, myRefCGindices=refCGindices) # see 'Example 2' - http://python.omics.wiki/multiprocessing_map/multiprocessing_partial_function_multiple_arguments, this is to use 1 dynamic fxn arguments and these 2 static arguments for multiprocessing

            gc.collect()
            chunksAnalyzingReads=p.map(analyzeReads_multipleArgs, intermedFileNames)
            p.close()
            p.join()

            print("exiting pool",utilities.getDate()) # ADDED

            #### merge all the locusUMIreadDicts and readDepthDicts
            for i in range(1,len(chunksAnalyzingReads),1): # we are starting at second element in list; we are appending everything into the first dictionary in this list
                mergeLocusUMIreadDicts(chunksAnalyzingReads[0].locusUMIreadsDict,chunksAnalyzingReads[i].locusUMIreadsDict)
                del chunksAnalyzingReads[i].locusUMIreadsDict # free up some RAM

                mergeReadDepthDicts(chunksAnalyzingReads[0].readDepthDict,chunksAnalyzingReads[i].readDepthDict)
                del chunksAnalyzingReads[i].readDepthDict # free up some RAM

            # figure out how many loci to analyze in each process
            numLociPerChunk=len(refLociNames)//processes

            splitUpMethInput=[]
            splitUpLociNames=[]
            readDict=chunksAnalyzingReads[0].readDepthDict

            # free up some RAM
            del chunksAnalyzingReads[0].readDepthDict

            myIndex=0
            for i in range(0,processes,1):
                tempDict={}
                tempLoci=[]

                if i != (processes-1):
                    tempLoci=refLociNames[myIndex:myIndex+numLociPerChunk]
                    myIndex+=numLociPerChunk
                else:
                    tempLoci=refLociNames[myIndex:] # last chunk of data will contain data for remaining loci

                for locus in tempLoci:
                    if locus in chunksAnalyzingReads[0].locusUMIreadsDict:
                        tempDict[locus]=chunksAnalyzingReads[0].locusUMIreadsDict[locus]
                    #else:

                splitUpLociNames.append(tempLoci)
                splitUpMethInput.append(MethylationChunkInput(myLocusUMIreadsDict=tempDict, myRefCGindices=refCGindices, myConsensusCutoff=consensus, mySampleName=sampleName, myTempNumber=i, myOutputDir=samDir, myUMIcollapseThreshold=umiCollapseThreshold, myUMIcollapseAlg=umiCollapseAlg))


            # clear up some RAM
            del chunksAnalyzingReads

            # finalResults=[]
            print("entering pool to analyze molecules for sample",sampleName,utilities.getDate())
            p = Pool(processes)

            # finalResults=p.map(analyzeMethylation, splitUpMethInput)

            gc.collect()
            p.map(callMolecules, splitUpMethInput)
            p.close()
            p.join()                        
            print("exiting pool",utilities.getDate())

            # combining temp output files
            print("combining temp output files (ALLMOLECULES)",sampleName,utilities.getDate())
            # outputFileName=samDir+sampleName+"_allMolecules_"+str(consensus)+"consensus.txt"
            outputFileName=samDir+sampleName+"_allMolecules.txt"

            with open(outputFileName,"w") as outputFile:
                myHeaderline="sample\tlocus\tUMI\tconsensusZstring_{0}\tmethCGs\ttotalCGs\tmethRatio\tnumReads\n".format(str(consensus))
                outputFile.write(myHeaderline)

            print("combining temp files to create", outputFileName.split('/')[-1], utilities.getDate())
            catTempFilesCommand="cat "+samDir+sampleName+"_TEMP_ALLMOLECULES_* >> "+outputFileName
            os.system(catTempFilesCommand)

            print("removing temp files for",sampleName,utilities.getDate())
            rmTempFilesCommand="rm "+samDir+sampleName+"_TEMP*"
            os.system(rmTempFilesCommand)

            print("compressing ",samFileName,utilities.getDate())
            if os.path.exists(samFileName+'.gz'): # remove sampleName.sam.gz if it already exists (means you are running this again). Need to remove so pigz will work properly without needing user input
                os.system('rm '+samFileName+'.gz')
            compressSamFileCommand="pigz "+samFileName
            os.system(compressSamFileCommand)

            print('\n**************')
            print('Done calling molecules for',sampleName,utilities.getDate())
            print('**************\n')
