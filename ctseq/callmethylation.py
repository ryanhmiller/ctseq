import os
import glob
from multiprocessing import Pool
from functools import partial
import subprocess
# from collections import defaultdict
from . import utilities

class Locus:
  def __init__(self):
      self.name=""
      self.methMolecules=0
      self.unmethMolecules=0
      self.totalMolecules=0
      self.totalReads=0      

def callMethylation(fileName,myBlacklistUMI,myCisCGcutoff,myMoleculeThreshold):
    locusDict={} # fragName:locusObject
    sampleName=fileName.split("_")[0]
    
    with open(fileName,"r") as infile:
        # infile.readline()

        # methylationDict={} # frag:%methylated
        # readDepthDict={} # frag:readDepth
        # umiDict={} # frag:molDepth
        # methMolDict={} # frag:#methMol

        for line in infile:

            line=line.strip("\n").split("\t")            
            # print(line)
            if line[0] != "sample" and line[1] != "locus": # making sure we're not reading the header line
                umi=line[2]
                zString=line[3]                

                if zString!="NA" and umi != myBlacklistUMI:
                    try:
                        locus=line[1]       
                        # totalCGs=int(line[5])
                        methRatio=float(line[6])
                        #molecules=int(line[5])
                        #molDepth=int(line[6])
                        readDepth=int(line[7])
                        # totalReadsPerSample[sampleName]+=readDepth
                    except:
                        print(line)



                    ##### NEED TO COUNT UMI DEPTH + READ DEPTH #####

                    # if methRatio >= cgMethylationCutoff:
                    #
                    #     print locus, str(molecules),str(methRatio), str(cgMethylationCutoff)


                    if locus not in locusDict:
                        newLocus=Locus()
                        newLocus.name=locus
                        # newLocus.totalCGs=totalCGs
                        #newLocus.molDepth=molDepth
                        #newLocus.totalReads=readDepth
                        locusDict[locus]=newLocus

                    locusDict[locus].totalReads+=readDepth


                    #### updated 4/15/20 - different read depth cutoffs for meth and unmeth mols
                    # if methRatio >= methCutoff and readDepth >= methCoverageCutoff: # this is a methylated molecule and passes read depth cutoff for meth molecules
                    #     locusDict[locus].molDepth+=1
                    #     locusDict[locus].methMolecules+=1
                    # elif methRatio < methCutoff and readDepth >= unmethCoverageCutoff: # this is a methylated molecule and passes read depth cutoff for meth molecules
                    #     locusDict[locus].molDepth+=1
                    #     locusDict[locus].unmethMolecules+=1

                    if readDepth >= myMoleculeThreshold: # checking if we can count this as a unique molecule
                        locusDict[locus].totalMolecules+=1
                        if methRatio >= myCisCGcutoff:
                            locusDict[locus].methMolecules+=1
                        else:
                            locusDict[locus].unmethMolecules+=1

    return(locusDict)


def combineLocusData(myListOfLociDicts):
    combinedLociDict={}

    for tempLocusDict in myListOfLociDicts:
        for myLocus in tempLocusDict:
            if myLocus not in combinedLociDict:             
                newLocus=Locus()
                newLocus.name=myLocus
                combinedLociDict[myLocus]=newLocus

            combinedLociDict[myLocus].methMolecules+=tempLocusDict[myLocus].methMolecules
            combinedLociDict[myLocus].unmethMolecules+=tempLocusDict[myLocus].unmethMolecules
            combinedLociDict[myLocus].totalMolecules+=tempLocusDict[myLocus].totalMolecules
            combinedLociDict[myLocus].totalReads+=tempLocusDict[myLocus].totalReads   

    return(combinedLociDict)


def run(args):
    refDir=args.refDir
    fileDir=args.dir
    processes=int(args.processes)
    cisCGcutoff=float(args.cisCG)
    moleculeThreshold=int(args.moleculeThreshold)
    fileExt='_allMolecules.txt'
    tempFileBase='_allMolecules_TEMP_'

    blacklistUMI=''

    #############
    # arg check #
    #############
    # check that paths are valid
    utilities.validDir(refDir)
    utilities.validDir(fileDir)

    # check that refDir/fileDir actually have the necessary files in them
    utilities.fileCheck(refDir,".fa")
    utilities.fileCheck(fileDir,fileExt)

    # make sure file dir has '/' on end
    if fileDir[-1]!='/':
        fileDir+='/'
    ####################
    ### SET UP STUFF ###
    ####################
    # get reference .fa file
    refFileName=utilities.getReferenceFile(refDir)


    # get list of reference fragment names
    refFrags=[]

    with open(refFileName, "r") as refFile:
        lines=refFile.readlines()
        for i in range(0,len(lines),2):
            fragName=lines[i].strip("\n")[1:]
            # seq=lines[i+1].strip("\n")
            # seq=seq.upper()
            refFrags.append(fragName)


    os.chdir(fileDir)

    listOfFiles=utilities.getFiles(fileDir,fileExt)    
    
    # ####################
    mySamples=[] # builds out list of all samples in directory

    lociDictsAllSamples={}

    for fileName in listOfFiles:
        # split input file into N (# processes) 
        sampleName=fileName.split("_")[0]
        mySamples.append(sampleName)

        print('\n**************')
        print('Calling methylation for',sampleName,utilities.getDate())
        print('**************\n')

        # grab number of lines in file
        totalLines = int(subprocess.getoutput('wc -l '+fileName).split(' ')[0])
        linesPerTempFile=(totalLines//processes)+1
        
        # split main file into temp files (N=number of processes)        
        print('splitting',fileName,'into',str(processes),'temp files',utilities.getDate())
        splitFileCmd='split -l '+str(linesPerTempFile)+' -d '+sampleName+fileExt+' --additional-suffix=.txt '+sampleName+tempFileBase
        os.system(splitFileCmd)


        tempFiles=glob.glob(sampleName+'_allMolecules_TEMP_*')
        # print(tempFiles)
        print("entering pool to analyze methylation for sample",sampleName,utilities.getDate())

        p = Pool(processes)
        callMethylation_multipleArgs=partial(callMethylation, myBlacklistUMI=blacklistUMI,myCisCGcutoff=cisCGcutoff,myMoleculeThreshold=moleculeThreshold) # see 'Example 2' - http://python.omics.wiki/multiprocessing_map/multiprocessing_partial_function_multiple_arguments, this is to use 1 dynamic fxn arguments and these 2 static arguments for multiprocessing
        methylResults=p.map(callMethylation_multipleArgs, tempFiles)
        print("exiting pool",utilities.getDate())

        combinedLociDict=combineLocusData(myListOfLociDicts=methylResults)

        lociDictsAllSamples[sampleName]=combinedLociDict
        # break
    print(list(lociDictsAllSamples.keys()))



        # splitForwardFileCmd=['split','-l',str(linesPerFile),'-d',forwardFile,'--additional-suffix=.fastq',sampleName+'_forward_TEMP_']


        # with open(fileName,"r") as infile:
        #     infile.readline()

        #     sampleName=fileName.split("_")[0]
        #     mySamples.append(sampleName)

        #     locusDict={} # fragName:locusObject

        #     methylationDict={} # frag:%methylated
        #     readDepthDict={} # frag:readDepth
        #     umiDict={} # frag:molDepth
        #     methMolDict={} # frag:#methMol

        #     for line in infile:

        #         line=line.strip("\n").split("\t")

        #         umi=line[2]
        #         zString=line[3]                

        #         if zString!="NA" and umi != blacklistUMI:
        #             locus=line[1]       
        #             totalCGs=int(line[5])
        #             methRatio=float(line[6])
        #             #molecules=int(line[5])
        #             #molDepth=int(line[6])
        #             readDepth=int(line[7])
        #             totalReadsPerSample[sampleName]+=readDepth


        #             ##### NEED TO COUNT UMI DEPTH + READ DEPTH #####

        #             # if methRatio >= cgMethylationCutoff:
        #             #
        #             #     print locus, str(molecules),str(methRatio), str(cgMethylationCutoff)


        #             if locus not in locusDict:
        #                 newLocus=Locus()
        #                 newLocus.name=locus
        #                 newLocus.totalCGs=totalCGs
        #                 #newLocus.molDepth=molDepth
        #                 #newLocus.totalReads=readDepth
        #                 locusDict[locus]=newLocus

        #             locusDict[locus].totalReads+=readDepth


        #             #### updated 4/15/20 - different read depth cutoffs for meth and unmeth mols
        #             # if methRatio >= methCutoff and readDepth >= methCoverageCutoff: # this is a methylated molecule and passes read depth cutoff for meth molecules
        #             #     locusDict[locus].molDepth+=1
        #             #     locusDict[locus].methMolecules+=1
        #             # elif methRatio < methCutoff and readDepth >= unmethCoverageCutoff: # this is a methylated molecule and passes read depth cutoff for meth molecules
        #             #     locusDict[locus].molDepth+=1
        #             #     locusDict[locus].unmethMolecules+=1

        #             if readDepth >= cisCGcutoff: # checking if we can count this as a unique molecule
        #                 locusDict[locus].molDepth+=1
        #                 if methRatio >= methCutoff:
        #                     locusDict[locus].methMolecules+=1
        #                 else:
        #                     locusDict[locus].unmethMolecules+=1






