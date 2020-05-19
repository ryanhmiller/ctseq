import os
import sys
import subprocess
import glob
from . import utilities

#### static variables (make this more flexible later on)
# args.forwardAdapter='AGTGTGGGAGGGTAGTTGGTGTT' #rcHA2 AGTGTGGGAGGGTAGTTGGTGT*T, forwardRead_3_primer_adapter
# args.reverseAdapter='ACTCCCCACCTTCCTCATTCTCTAAGACGGTGT' #rcHA1, reverseRead_3_primer_adapter
# cutadaptCores=str(18)
# bismarkCores=str(6)

# readsPerFile=1500000


def cutAdapters(args,fExt,rExt,newForExt,newRevExt):
    print('\n**************')
    print('Trimming adapters with Cutadapt '+utilities.getDate())
    print('**************\n')

    cutadaptCores=str(args.cutadaptCores)
    # loop through fastq files
    allForwardFiles=utilities.naturalSort(glob.glob('*'+fExt))

    for forwardFile in allForwardFiles:
        sampleName=forwardFile.split("_")[0]

        reverseFile=glob.glob(sampleName+'*'+rExt)[0]

        forwardTrimmedFile=sampleName+newForExt
        reverseTrimmedFile=sampleName+newRevExt

        # print(forwardTrimmedFile,forwardFile)
        # print(reverseTrimmedFile,reverseFile)
        # print(args.forwardAdapter,args.reverseAdapter)
        cutadaptCmd=['cutadapt','-a',args.forwardAdapter,'-A',args.reverseAdapter,'--cores',cutadaptCores,'--minimum-length',str(10),'-o',forwardTrimmedFile,'-p',reverseTrimmedFile,forwardFile,reverseFile]


        externalProcess = subprocess.Popen(cutadaptCmd)
        exitCode = externalProcess.wait()
        error = externalProcess.stderr

        if exitCode != 0:
            print('\n**ERROR**')
            print('Something went wrong with cutadapt. Exiting...')
            sys.exit()

        print('Removing ',forwardFile,'and',reverseFile,utilities.getDate())
        removeFiles([forwardFile,reverseFile]) # removing input files because no longer need them

    print('\n**************')
    print('Done trimming adapters with Cutadapt '+utilities.getDate())
    print('**************\n')


def catFiles(listOfFiles,finalFileName):

    catCmd=['cat']+listOfFiles

    with open(finalFileName,'w') as outputFile:
        catProcess = subprocess.Popen(catCmd,stdout=outputFile)
        catExitCode = catProcess.wait()
        catError = catProcess.stderr

        if catExitCode != 0:
            print('\n**ERROR**')
            print('Something went wrong with cat-ing '+finalFileName+' files. Exiting...')
            sys.exit()

def removeFiles(listOfFiles):
    rmCmd=['rm']+listOfFiles
    rmProcess = subprocess.Popen(rmCmd)
    rmExitCode = rmProcess.wait()
    rmError = rmProcess.stderr

    if rmExitCode != 0:
        print('\n**ERROR**')
        print('Something went wrong with rm-ing these files. Exiting...')
        print(listOfFiles)
        print('\n')
        sys.exit()


def align(args,myRefDir,newForExt,newRevExt):
    print('\n**************')
    print('Aligning reads with Bismark '+utilities.getDate())
    print('**************\n')

    sampleNames=[]

    ############################
    # loop through fastq files #
    ############################
    allForwardFiles=utilities.naturalSort(glob.glob('*'+newForExt))

    for forwardFile in allForwardFiles:
        sampleName=forwardFile.split("_")[0]
        sampleNames.append(sampleName)

        reverseFile=glob.glob(sampleName+'*'+newRevExt)[0]

        ###############################################################
        ## split fastq files up here if needed to run through Bismark #
        ###############################################################

        linesPerFile=args.readsPerFile*4

        splitForwardFileCmd=['split','-l',str(linesPerFile),'-d',forwardFile,'--additional-suffix=.fastq',sampleName+'_forward_TEMP_']

        splitReverseFileCmd=['split','-l',str(linesPerFile),'-d',reverseFile,'--additional-suffix=.fastq',sampleName+'_reverse_TEMP_']

        # splitProcessForward = subprocess.Popen(splitForwardFileCmd)
        # forwardSplitExitCode = splitProcessForward.wait()
        # forwardSplitError = splitProcessForward.stderr

        os.system(' '.join(splitForwardFileCmd))

        # splitProcessReverse = subprocess.Popen(splitReverseFileCmd)
        # reverseSplitExitCode = splitProcessReverse.wait()
        # reverseSplitError = splitProcessReverse.stderr
        #

        os.system(' '.join(splitReverseFileCmd))

        # if forwardSplitExitCode != 0 or reverseSplitExitCode != 0:
        #     print('\n**ERROR**')
        #     print('Something went wrong with splitting up files before running Bismark. Exiting...')
        #     sys.exit()

        ################################################
        ## run split up fastqs and run through Bismark #
        ################################################

        for forwardSplitFile in utilities.naturalSort(glob.glob(sampleName+'_forward_TEMP_*.fastq')):
            subsetNum=forwardSplitFile.split(".")[0].split("_")[3]
            reverseSplitFile=forwardSplitFile.replace("forward","reverse")

            bismarkCmd=['bismark','--bowtie2','--multicore',str(args.bismarkCores),'--ambiguous','--unmapped','--genome_folder',
                        myRefDir,'-1',forwardSplitFile,'-2',reverseSplitFile]

            print('Aligning '+forwardSplitFile+' with Bismark '+utilities.getDate())

            os.system(' '.join(bismarkCmd))
            # externalProcess = subprocess.Popen(bismarkCmd)
            # exitCode = externalProcess.wait()
            # error = externalProcess.stderr
            #
            # if exitCode != 0:
            #     print('\n**ERROR**')
            #     print('Something went wrong with Bismark. Exiting...')
            #     sys.exit()

        #####################################################
        # combine bams (if only one bam, still renames bam) #
        #####################################################
        thisSampleBams=utilities.naturalSort(glob.glob(sampleName+'_forward_TEMP_*_bismark_bt2_pe.bam'))
        combinedBamFile=sampleName+'.bam'

        combineBamsCmd=['samtools','cat','-o',combinedBamFile]+thisSampleBams

        print('Combining '+sampleName+' bam files '+utilities.getDate())

        os.system(' '.join(combineBamsCmd))
        # combineBamsProcess = subprocess.Popen(combineBamsCmd)
        # combineBamsExitCode = combineBamsProcess.wait()
        # combineBamsError = combineBamsProcess.stderr
        #
        # if combineBamsExitCode != 0:
        #     print('\n**ERROR**')
        #     print('Something went wrong combining bam files. Exiting...')
        #     sys.exit()

        ###########################
        # convert bam to sam file #
        ###########################

        combinedSamFile=sampleName+'.sam'

        bam2samCmd=['samtools','view','-@',str(args.cutadaptCores),'-h',combinedBamFile]

        print('Converting '+sampleName+' bam file to sam format '+utilities.getDate())

        # output of converting bam to sam goes to stdout so we need to open a file and push all that data there
        with open(combinedSamFile,'w') as mysamfile:
            bam2samProcess = subprocess.Popen(bam2samCmd,stdout=mysamfile)
            bam2samExitCode = bam2samProcess.wait()
            bam2samError = bam2samProcess.stderr

            if bam2samExitCode != 0:
                print('\n**ERROR**')
                print('Something went wrong with converting .bam to .sam after Bismark for '+sampleName+'. Exiting...')
                sys.exit()

        #############################################
        # combine ambiguous and unmapped reads file #
        #############################################
        # ambig files
        print('Combining '+sampleName+' ambiguous reads files '+utilities.getDate())
        forwardAmbigFiles=utilities.naturalSort(glob.glob(sampleName+'_forward_TEMP_*.fastq_ambiguous_reads_1.fq.gz'))
        catFiles(forwardAmbigFiles,(sampleName+'_forward.fastq_ambiguous_reads_1.fq.gz'))

        reverseAmbigFiles=utilities.naturalSort(glob.glob(sampleName+'_reverse_TEMP_*.fastq_ambiguous_reads_2.fq.gz'))
        catFiles(reverseAmbigFiles,(sampleName+'_reverse.fastq_ambiguous_reads_2.fq.gz'))

        # unmapped files
        print('Combining '+sampleName+' unmapped reads files '+utilities.getDate())
        forwardUnmappedFiles=utilities.naturalSort(glob.glob(sampleName+'_forward_TEMP_*.fastq_unmapped_reads_1.fq.gz'))
        catFiles(forwardUnmappedFiles,(sampleName+'_forward.fastq_unmapped_reads_1.fq.gz'))

        reverseUnmappedFiles=utilities.naturalSort(glob.glob(sampleName+'_reverse_TEMP_*.fastq_unmapped_reads_2.fq.gz'))
        catFiles(reverseUnmappedFiles,(sampleName+'_reverse.fastq_unmapped_reads_2.fq.gz'))


        # rm unnecessary files for this sample
        print('Removing temp files for',sampleName,utilities.getDate())
        removeFiles(glob.glob(sampleName+'_*TEMP*.bam'))
        removeFiles(glob.glob(sampleName+'_*TEMP*.fastq'))
        removeFiles(glob.glob(sampleName+'_*TEMP*.fq.gz'))
        removeFiles(glob.glob(sampleName+'_*ReadsWithUMIsTRIMMED.fastq'))
        removeFiles([sampleName+".bam"])

        # print('Renaming bismark report file',sampleName,utilities.getDate())
        # bismarkReportFileName=glob.glob(sampleName+'_*bismark_bt2_PE_report.txt')[0]
        # newBismarkReportFileName=sampleName+'_bismarkReport.txt'
        # renameBismarkReportCmd='mv '+bismarkReportFileName+' '+newBismarkReportFileName
        # os.system(renameBismarkReportCmd)

        print('\nDone aligning sample '+sampleName+' '+utilities.getDate()+'\n')


    print('\n**************')
    print('Done aligning reads with Bismark '+utilities.getDate())
    print('**************\n')


def run(args):
    refDir=args.refDir
    dir=args.dir
    cutadaptCores=int(args.cutadaptCores)
    bismarkCores=int(args.bismarkCores)
    readsPerFile=int(args.readsPerFile)

    fExt='_forwardReadsWithUMIs.fastq'
    rExt='_reverseReadsWithUMIs.fastq'

    newForExt='_forwardReadsWithUMIsTRIMMED.fastq'
    newRevExt='_reverseReadsWithUMIsTRIMMED.fastq'

    #############
    # arg check #
    #############
    # check that paths are valid
    refDir=utilities.validDir(refDir)
    dir=utilities.validDir(dir)

    # check that dir actually has fastq files in it
    utilities.fileCheck(dir,fExt)
    utilities.fileCheck(dir,rExt)

    #########################
    # continue with command #
    #########################

    # make methylation reference genome with bismark
    # try:
    os.chdir(dir)

    cutAdapters(args,fExt=fExt,rExt=rExt,newForExt=newForExt,newRevExt=newRevExt)

    align(args,myRefDir=refDir,newForExt=newForExt,newRevExt=newRevExt)
