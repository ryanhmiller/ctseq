import os
import sys
import glob
from . import utilities


def run(args):
    dir=args.dir
    runType=args.umiType
    umiLength=int(args.umiLength)
    forwardExt=args.forwardExt
    reverseExt=args.reverseExt
    umiExt=args.umiExt

    # if dir[-1]!='/':
    #     dir+='/'

    ##############################################
    # arg check - make sure this is a valid path #
    ##############################################
    print('\n**************')
    print('ADDING UMIs',utilities.getDate())
    print('**************\n')

    dir=utilities.validDir(dir)
    utilities.fileCheck(dir,forwardExt)
    utilities.fileCheck(dir,reverseExt)

    os.chdir(dir)

    if runType=='separate' and umiExt=='NOTSPECIFIED':
        print('\n**ERROR**')
        print('Since your UMIs are in separate fastq files, please specify the unique file extension of the UMI fastq files for the \'--umiExt\' flag (e.g. R2_001.fastq OR R2_001.fastq.gz)')
        print('**EXITING**')
        sys.exit()
    elif runType=='inline' and umiExt!='NOTSPECIFIED':
        print('\n**ERROR**')
        print('You specified that the UMIs are \'inline\' for the \'--type\' flag yet you also used the \'--umiExt\' flag implying you have separate file with UMIs in it. Did you mean to use \'separate\' with the \'type\' flag instead?')
        print('**EXITING**')
        sys.exit()
    elif runType=='separate' and umiExt!='NOTSPECIFIED':
        utilities.fileCheck(dir,umiExt)


    #############################
    # continue with adding umis #
    #############################
    allForwardFiles=utilities.naturalSort(glob.glob('*'+forwardExt))

    for forwardFileName in allForwardFiles:
        sampleName=forwardFileName.split("_")[0]

        reverseFileName=glob.glob(sampleName+'_*'+reverseExt)[0]

        # unzip forward/reverse files if needed
        forwardFileName=decompressFile(forwardFileName)
        reverseFileName=decompressFile(reverseFileName)

        print('****')
        if runType=='separate':
            umiFileName=glob.glob(sampleName+'_*'+umiExt)[0]
            umiFileName=decompressFile(umiFileName)

            # make dictionary of umis
            umiDict={}
            with open(umiFileName,"r") as umiFile:
                umiData=umiFile.readlines()

                for i in range(0,len(umiData),1):
                    if i % 4 == 0:
                        line=umiData[i]
                        readName=line.strip("\n").split(" ")[0]
                        umi=umiData[i+1].strip("\n")[:umiLength] # grab the specified length of umi
                        umiDict[readName]=umi

            # add UMIs to forward and reverse file
            print('adding UMIs to fastq file with forward reads for',forwardFileName,utilities.getDate())
            addUMIs_separate(myUMIdict=umiDict,myFastqFileName=forwardFileName,myBaseName='forward')
            print('adding UMIs to fastq file with reverse reads for',reverseFileName,utilities.getDate())
            addUMIs_separate(myUMIdict=umiDict,myFastqFileName=reverseFileName,myBaseName='reverse')

            print('compressing UMI file for sample',sampleName,utilities.getDate())
            os.system('pigz '+umiFileName)

        elif runType=='inline':
            print('adding inline UMIs to forward reads for',forwardFileName,utilities.getDate())
            addUMIs_inline(myFastqFileName=forwardFileName,myBaseName='forward',myUMIlength=umiLength)
            print('adding inline UMIs to reverse reads for',reverseFileName,utilities.getDate())
            addUMIs_inline(myFastqFileName=reverseFileName,myBaseName='reverse',myUMIlength=umiLength)

        # compress input files
        print('compressing forward file for sample',sampleName,utilities.getDate())
        os.system('pigz '+forwardFileName)
        print('compressing reverse file for sample',sampleName,utilities.getDate())
        os.system('pigz '+reverseFileName)
        print('')


    print('\n**************')
    print('DONE ADDING UMIs',utilities.getDate())
    print('**************\n')


def addUMIs_separate(myUMIdict,myFastqFileName,myBaseName):
    sampleName=myFastqFileName.split('_')[0]
    outputFileName=sampleName+'_'+myBaseName+'ReadsWithUMIs.fastq'
    # add umis to data file
    with open(myFastqFileName, 'r') as fastqFile:

        with open(outputFileName,'w') as outFile:
            lineNumber=0

            for line in fastqFile:
                if lineNumber % 4 ==0:
                    line = line.strip('\n').split(' ')
                    myReadName=line[0]
                    sampleBarcode=line[1].split(':')[-1] # just grabbing the sample umi 3:N:0:CGTGTAAT
                    line=myReadName+'-'+sampleBarcode+'+'+myUMIdict[myReadName]+'\n'

                outFile.write(line)
                lineNumber+=1


def addUMIs_inline(myFastqFileName,myBaseName,myUMIlength):
    sampleName=myFastqFileName.split('_')[0]
    outputFileName=sampleName+'_'+myBaseName+'ReadsWithUMIs.fastq'

    with open(myFastqFileName, 'r') as dataFile:
        with open(outputFileName,'w') as outFile:
            lineNumber=0

            for line in dataFile:
                if lineNumber % 4 ==0:
                    line = line.strip('\n').split(' ')
                    umi=line[0].split(':')[-1] # grabbing the UMI
                    umi=umi[:myUMIlength] # make sure UMI is correct length
                    myReadName=':'.join(line[0].split(':')[:-1]) # grabbing read name minus the UMI
                    sampleBarcode=line[1].split(':')[-1] # just grabbing the sample barcode 3:N:0:CGTGTAAT
                    line=myReadName+'-'+sampleBarcode+'+'+umi+'\n'

                outFile.write(line)
                lineNumber+=1

# this fxn decompresses file if file ends in '.gz'
def decompressFile(fileName):
    if fileName[-3:]=='.gz':
        os.system('pigz -d '+fileName)
        fileName=fileName[:-3]

    return(fileName)
