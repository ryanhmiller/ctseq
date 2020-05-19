import os
import glob
from multiprocessing import Pool
from functools import partial
import subprocess
# from collections import defaultdict
from . import utilities
import sys

class Locus:
  def __init__(self):
      self.name=''
      self.methMolecules=0
      self.unmethMolecules=0
      self.totalMolecules=0
      self.totalReads=0

def callMethylation(fileName,myBlacklistUMI,myCisCGcutoff,myMoleculeThreshold):
    locusDict={} # fragName:locusObject
    sampleName=fileName.split('_')[0]

    with open(fileName,'r') as infile:

        for line in infile:

            line=line.strip('\n').split('\t')
            # print(line)
            if line[0] != 'sample' and line[1] != 'locus': # making sure we're not reading the header line
                umi=line[2]
                zString=line[3]

                if zString!='NA' and umi != myBlacklistUMI:
                    try:
                        locus=line[1]
                        methRatio=float(line[6])
                        readDepth=int(line[7])
                        # totalReadsPerSample[sampleName]+=readDepth
                    except:
                        print(line)

                    if locus not in locusDict:
                        newLocus=Locus()
                        newLocus.name=locus
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

def writeReport(myListOfLociDicts,reportType,myRefFrags,mySampleNames,myRunName):
    #write header line
    reportName=''
    if reportType == 'methMolecules':
        reportName='methylatedMolecules'
    elif reportType == 'totalMolecules':
        reportName='totalMolecules'
    elif reportType == 'totalReads':
        reportName='totalReads'
    elif reportType == 'methRatio':
        reportName='methylationRatio'

    outputFileName=myRunName+'_'+reportName+'.txt'

    with open(outputFileName, 'w') as outputFile:
        # write header line
        headerLine=['Locus']+mySampleNames
        outputFile.write('\t'.join(headerLine)+'\n')

        # write out results, locus by locus
        for locus in myRefFrags:
            outputLine=[locus]
            for sample in mySampleNames:
                elementToAdd=''
                if locus in myListOfLociDicts[sample]:
                    if reportType == 'methRatio':
                        methMol=myListOfLociDicts[sample][locus].methMolecules
                        totalMol=myListOfLociDicts[sample][locus].totalMolecules
                        if totalMol == 0:
                            elementToAdd='NA'
                        else:
                            elementToAdd=str(methMol/totalMol)
                    elif reportType == 'totalMolecules' and myListOfLociDicts[sample][locus].totalMolecules==0:
                        elementToAdd='NA'
                    elif reportType == 'methMolecules' and myListOfLociDicts[sample][locus].totalMolecules==0:
                        elementToAdd='NA'
                    else:
                        elementToAdd=str(getattr(myListOfLociDicts[sample][locus],reportType)) # for more info: https://stackoverflow.com/questions/2157035/accessing-an-attribute-using-a-variable-in-python
                else:
                    elementToAdd='NA' # if locus not in dictionary, there were no molecules for that locus in that sample

                outputLine.append(elementToAdd)
            outputFile.write('\t'.join(outputLine)+'\n')



class sampleStats:
    def __init__(self):
        self.bismarkAlignedReads=0
        self.bismarkTotalReads=0

        self.methCpG=0
        self.unmethCpG=0

        self.methCHG=0
        self.unmethCHG=0

        self.methCHH=0
        self.unmethCHH=0

        self.methUnknown=0
        self.unmethUnknown=0

def writeRunStatsReport(myListOfLociDicts,myRefFrags,mySampleNames,myRunName):

    outputFileName=myRunName+'_runStatistics.txt'

    # grab stats from bismark report files
    bismarkStatsContainer={}
    bismarkReportFiles=utilities.getFiles(path=os.getcwd(),fileExt='report.txt')

    for reportFileName in bismarkReportFiles:
       # print(reportFileName)
       with open(reportFileName,'r') as reportFile:
            currSample=""
            currAlignedReads=0
            currTotalReads=0

            currMethCpG=0
            currUnmethCpG=0

            currMethCHG=0
            currUnmethCHG=0

            currMethCHH=0
            currUnmethCHH=0

            currMethUnknown=0
            currUnmethUnknown=0

            # percent=""
            for line in reportFile:
                line=line.strip('\n').split(':')
                if line[0]=='Bismark report for':
                    # currSample=line[1].split('_')[0].split(' ')[1]
                    currSample=line[1].split(' ')[1].split('/')[-1].split('_')[0]
                    #print(currSample,'\n')

                elif line[0]=='Number of paired-end alignments with a unique best hit':
                    currAlignedReads=int(line[1].split('\t')[1])
                elif line[0]=='Sequence pairs analysed in total':
                    currTotalReads=int(line[1].split('\t')[1])


                elif line[0]=="Total methylated C's in CpG context":
                    currMethCpG=int(line[1].split('\t')[1])
                elif line[0]=="Total unmethylated C's in CpG context":
                    currUnmethCpG=int(line[1].split('\t')[1])


                elif line[0]=="Total methylated C's in CHG context":
                    currMethCHG=int(line[1].split('\t')[1])
                elif line[0]=="Total unmethylated C's in CHG context":
                    currUnmethCHG=int(line[1].split('\t')[1])


                elif line[0]=="Total methylated C's in CHH context":
                    currMethCHH=int(line[1].split('\t')[1])
                elif line[0]=="Total unmethylated C's in CHH context":
                    currUnmethCHH=int(line[1].split('\t')[1])


                elif line[0]=="Total methylated C's in Unknown context":
                    currMethUnknown=int(line[1].split('\t')[1])
                elif line[0]=="Total unmethylated C's in Unknown context":
                    currUnmethUnknown=int(line[1].split('\t')[1])


                    if currSample not in bismarkStatsContainer:
                        currObj=sampleStats()
                        bismarkStatsContainer[currSample]=currObj

                    bismarkStatsContainer[currSample].bismarkAlignedReads+=currAlignedReads
                    bismarkStatsContainer[currSample].bismarkTotalReads+=currTotalReads

                    bismarkStatsContainer[currSample].methCpG+=currMethCpG
                    bismarkStatsContainer[currSample].unmethCpG+=currUnmethCpG

                    bismarkStatsContainer[currSample].methCHG+=currMethCHG
                    bismarkStatsContainer[currSample].unmethCHG+=currUnmethCHG

                    bismarkStatsContainer[currSample].methCHH+=currMethCHH
                    bismarkStatsContainer[currSample].unmethCHH+=currUnmethCHH

                    bismarkStatsContainer[currSample].methUnknown+=currMethUnknown
                    bismarkStatsContainer[currSample].unmethUnknown+=currUnmethUnknown


                    currSample=""
                    currAlignedReads=0
                    currTotalReads=0

                    currMethCpG=0
                    currUnmethCpG=0

                    currMethCHG=0
                    currUnmethCHG=0

                    currMethCHH=0
                    currUnmethCHH=0

                    currMethUnknown=0
                    currUnmethUnknown=0


    # iterate through myListOfLociDicts and calc # aligned reads/molecules for each sample
    #alignedReadsMolDict={} # sampleName : [alignedReads, alignedMolecules]

    # write out stats
    with open(outputFileName, 'w') as outputFile:
        # write header line
        headerLine=['sample','percentAlignedReads','alignedReads','alignedMol','methCpG','methCHG','methCHH','methUnknownCNorCHN']
        outputFile.write('\t'.join(headerLine)+'\n')

        # loop through and grab aligned reads/molecules stats for each sample and then write out all stats for that sample
        for sample in mySampleNames:
            sampleAlignedReads=0
            sampleAlignedMolecules=0

            for frag in myRefFrags:
                if frag in myListOfLociDicts[sample]:
                    sampleAlignedReads+=myListOfLociDicts[sample][frag].totalReads
                    sampleAlignedMolecules+=myListOfLociDicts[sample][frag].totalMolecules


            percentReadsAligned='0'
            methCpG='0'
            methCHG='0'
            methUnknown='0'

            if bismarkStatsContainer[sample].bismarkTotalReads != 0:
                percentReadsAligned=str(round((bismarkStatsContainer[sample].bismarkAlignedReads/bismarkStatsContainer[sample].bismarkTotalReads)*100,1))

            if bismarkStatsContainer[sample].methCpG !=0 or bismarkStatsContainer[sample].unmethCpG !=0:
                methCpG=str(round(((bismarkStatsContainer[sample].methCpG/(bismarkStatsContainer[sample].methCpG+bismarkStatsContainer[sample].unmethCpG))*100),1))

            if bismarkStatsContainer[sample].methCHG != 0 or bismarkStatsContainer[sample].unmethCHG != 0:
                methCHG=str(round(((bismarkStatsContainer[sample].methCHG/(bismarkStatsContainer[sample].methCHG+bismarkStatsContainer[sample].unmethCHG))*100),1))

            if bismarkStatsContainer[sample].methCHH != 0 or bismarkStatsContainer[sample].unmethCHH != 0:
                methCHH=str(round(((bismarkStatsContainer[sample].methCHH/(bismarkStatsContainer[sample].methCHH+bismarkStatsContainer[sample].unmethCHH))*100),1))

            if bismarkStatsContainer[sample].methUnknown != 0 or bismarkStatsContainer[sample].unmethUnknown != 0:
                methUnknown=str(round(((bismarkStatsContainer[sample].methUnknown/(bismarkStatsContainer[sample].methUnknown+bismarkStatsContainer[sample].unmethUnknown))*100),1))

            outputLine=[sample,percentReadsAligned,str(sampleAlignedReads),str(sampleAlignedMolecules),methCpG,methCHG,methCHH,methUnknown]

            outputFile.write('\t'.join(outputLine)+'\n')


def run(args):
    refDir=args.refDir
    fileDir=args.dir
    nameRun=args.nameRun
    processes=int(args.processes)
    cisCGcutoff=float(args.cisCG)
    moleculeThreshold=int(args.moleculeThreshold)
    fileExt='_allMolecules.txt'
    tempFileBase='_allMolecules_TEMP_'

    blacklistUMI='CGTCTTAGAGAA'

    #############
    # arg check #
    #############
    # check that paths are valid
    refDir=utilities.validDir(refDir)
    fileDir=utilities.validDir(fileDir)

    # check that refDir/fileDir actually have the necessary files in them
    utilities.fileCheck(refDir,'.fa')
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

    with open(refFileName, 'r') as refFile:
        lines=refFile.readlines()
        for i in range(0,len(lines),2):
            fragName=lines[i].strip('\n')[1:]
            # seq=lines[i+1].strip('\n')
            # seq=seq.upper()
            refFrags.append(fragName)


    os.chdir(fileDir)

    listOfFiles=utilities.getFiles(fileDir,fileExt)

    # ####################
    allSamples=[] # builds out list of all samples in directory

    lociDictsAllSamples={}

    print('\n**************')
    print('CALLING METHYLATION FOR RUN',nameRun,utilities.getDate())
    print('**************\n')

    for fileName in listOfFiles:
        # split input file into N (# processes)
        sampleName=fileName.split('_')[0]
        allSamples.append(sampleName)

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

        tempFiles=utilities.naturalSort(glob.glob(sampleName+'_allMolecules_TEMP_*'))
        # print(tempFiles)

        # entering pool to analyze the methylation
        print('entering pool to analyze methylation for sample',sampleName,utilities.getDate())
        p = Pool(processes)
        callMethylation_multipleArgs=partial(callMethylation, myBlacklistUMI=blacklistUMI,myCisCGcutoff=cisCGcutoff,myMoleculeThreshold=moleculeThreshold) # see 'Example 2' - http://python.omics.wiki/multiprocessing_map/multiprocessing_partial_function_multiple_arguments, this is to use 1 dynamic fxn arguments and these 2 static arguments for multiprocessing
        methylResults=p.map(callMethylation_multipleArgs, tempFiles)
        print('exiting pool',utilities.getDate())

        ## remove temp files
        print('removing temp files for',sampleName,utilities.getDate())
        rmTempFilesCmd='rm '+sampleName+'_allMolecules_TEMP_*'+'.txt'
        os.system(rmTempFilesCmd)

        combinedLociDict=combineLocusData(myListOfLociDicts=methylResults)

        lociDictsAllSamples[sampleName]=combinedLociDict
        # break


    print('writing reports for ',nameRun,utilities.getDate())
    writeReport(reportType='totalMolecules',myListOfLociDicts=lociDictsAllSamples,myRefFrags=refFrags,mySampleNames=allSamples,myRunName=nameRun)
    writeReport(reportType='methMolecules',myListOfLociDicts=lociDictsAllSamples,myRefFrags=refFrags,mySampleNames=allSamples,myRunName=nameRun)
    writeReport(reportType='totalReads',myListOfLociDicts=lociDictsAllSamples,myRefFrags=refFrags,mySampleNames=allSamples,myRunName=nameRun)
    writeReport(reportType='methRatio',myListOfLociDicts=lociDictsAllSamples,myRefFrags=refFrags,mySampleNames=allSamples,myRunName=nameRun)

    writeRunStatsReport(myListOfLociDicts=lociDictsAllSamples,myRefFrags=refFrags,mySampleNames=allSamples,myRunName=nameRun)
    print('finished writing reports for ',nameRun,utilities.getDate())


    print('\n**************')
    print('DONE CALLING METHYLATION FOR RUN',nameRun,utilities.getDate())
    print('**************\n')
