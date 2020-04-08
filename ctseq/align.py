import os
import sys
import subprocess
import glob
from . import utilities

#### static variables (make this more flexible later on)
forwardReadAdapter='AGTGTGGGAGGGTAGTTGGTGTT' #rcHA2 AGTGTGGGAGGGTAGTTGGTGT*T, forwardRead_3_primer_adapter
reverseReadAdapter='ACTCCCCACCTTCCTCATTCTCTAAGACGGTGT' #rcHA1, reverseRead_3_primer_adapter
cutadaptCores=str(18)
bismarkCores=str(6)
trimmedFileExt='.TRIMMED.fastq'
readsPerFile=1500000


def cutAdapters(args):
    print('\n**************')
    print('Trimming adapters with Cutadapt '+utilities.getDate())
    print('**************\n')

    # loop through fastq files
    allForwardFiles=glob.glob('*'+args.fExt)

    for forwardFile in allForwardFiles:
        sampleName=forwardFile.split("_")[0]

        reverseFile=glob.glob(sampleName+'*'+args.rExt)[0]

        forwardTrimmedFile=forwardFile+trimmedFileExt
        reverseTrimmedFile=reverseFile+trimmedFileExt

        # print(forwardTrimmedFile,forwardFile)
        # print(reverseTrimmedFile,reverseFile)
        # print(forwardReadAdapter,reverseReadAdapter)
        cutadaptCmd=['cutadapt','-a',forwardReadAdapter,'-A',reverseReadAdapter,'--cores',cutadaptCores,'--minimum-length',str(10),'-o',forwardTrimmedFile,'-p',reverseTrimmedFile,forwardFile,reverseFile]


        externalProcess = subprocess.Popen(cutadaptCmd)
        exitCode = externalProcess.wait()
        error = externalProcess.stderr

        if exitCode != 0:
            print('\n**ERROR**')
            print('Something went wrong with cutadapt. Exiting...')
            sys.exit()

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


def align(args):
    print('\n**************')
    print('Aligning reads with Bismark '+utilities.getDate())
    print('**************\n')

    sampleNames=[]

    ############################
    # loop through fastq files #
    ############################
    allForwardFiles=glob.glob('*'+args.fExt+trimmedFileExt)

    for forwardFile in allForwardFiles:
        sampleName=forwardFile.split("_")[0]
        sampleNames.append(sampleName)

        reverseFile=glob.glob(sampleName+'*'+args.rExt+trimmedFileExt)[0]

        ###############################################################
        ## split fastq files up here if needed to run through Bismark #
        ###############################################################

        linesPerFile=readsPerFile*4

        splitForwardFileCmd=['split','-l',str(linesPerFile),'-d',forwardFile,'--additional-suffix=.fastq',sampleName+'_forward_TEMP_']

        splitReverseFileCmd=['split','-l',str(linesPerFile),'-d',reverseFile,'--additional-suffix=.fastq',sampleName+'_reverse_TEMP_']

        splitProcessForward = subprocess.Popen(splitForwardFileCmd)
        forwardSplitExitCode = splitProcessForward.wait()
        forwardSplitError = splitProcessForward.stderr

        splitProcessReverse = subprocess.Popen(splitReverseFileCmd)
        reverseSplitExitCode = splitProcessReverse.wait()
        reverseSplitError = splitProcessReverse.stderr

        if forwardSplitExitCode != 0 or reverseSplitExitCode != 0:
            print('\n**ERROR**')
            print('Something went wrong with splitting up files before running Bismark. Exiting...')
            sys.exit()

        ################################################
        ## run split up fastqs and run through Bismark #
        ################################################

        for forwardSplitFile in utilities.naturalSort(glob.glob(sampleName+'_forward_TEMP_*.fastq')):
            subsetNum=forwardSplitFile.split(".")[0].split("_")[3]
            reverseSplitFile=forwardSplitFile.replace("forward","reverse")

            bismarkCmd=['bismark','--bowtie2','--multicore',str(bismarkCores),'--ambiguous','--unmapped','--genome_folder',
                        args.refDir,'-1',forwardSplitFile,'-2',reverseSplitFile]

            print('Aligning '+forwardSplitFile+' with Bismark '+utilities.getDate())
            externalProcess = subprocess.Popen(bismarkCmd)
            exitCode = externalProcess.wait()
            error = externalProcess.stderr

            if exitCode != 0:
                print('\n**ERROR**')
                print('Something went wrong with Bismark. Exiting...')
                sys.exit()

        #####################################################
        # combine bams (if only one bam, still renames bam) #
        #####################################################
        thisSampleBams=utilities.naturalSort(glob.glob(sampleName+'_forward_TEMP_*_bismark_bt2_pe.bam'))
        combinedBamFile=sampleName+'.bam'

        combineBamsCmd=['samtools','cat','-o',combinedBamFile]+thisSampleBams

        print('Combining '+sampleName+' bam files '+utilities.getDate())
        combineBamsProcess = subprocess.Popen(combineBamsCmd)
        combineBamsExitCode = combineBamsProcess.wait()
        combineBamsError = combineBamsProcess.stderr

        if combineBamsExitCode != 0:
            print('\n**ERROR**')
            print('Something went wrong combining bam files. Exiting...')
            sys.exit()

        ###########################
        # convert bam to sam file #
        ###########################

        combinedSamFile=sampleName+'.sam'

        bam2samCmd=['samtools','view','-@',cutadaptCores,'-h',combinedBamFile]

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
        print('Removing temp files for '+sampleName+' '+utilities.getDate())
        removeFiles(glob.glob(sampleName+'_*TEMP*'))
        removeFiles(glob.glob(sampleName+'_*TRIMMED*'))
        removeFiles([sampleName+".bam"])


        print('\nDone aligning sample '+sampleName+' '+utilities.getDate()+"\n")


        break

    print('\n**************')
    print('Done aligning reads with Bismark '+utilities.getDate())
    print('**************\n')


def run(args):
    refDir=args.refDir
    fqDir=args.fqDir
    fExt=args.fExt
    rExt=args.rExt

    #############
    # arg check #
    #############
    # check that paths are valid
    utilities.validDir(refDir)
    utilities.validDir(fqDir)

    # check that fqDir actually has fastq files in it
    utilities.fileCheck(fqDir,fExt)
    utilities.fileCheck(fqDir,rExt)

    #########################
    # continue with command #
    #########################

    # make methylation reference genome with bismark
    # try:
    os.chdir(fqDir)


    cutAdapters(args)

    align(args)
