import os
import inspect
import sys
from . import utilities

def checkInputIsUnique(fileList,fileExt,myFileDir):
    if len(fileList) > 1:
        print('ERROR: it looks like you have more than one *'+fileExt+' file at '+myFileDir)
        print('Please make sure there is only one of these files at this location')
        print('Exiting...')
        sys.exit()
    else:
        return(fileList[0])

def run(args):
    fileDir=args.dir
    fragInfoFile=args.fragInfo

    totalMolExt='_totalMolecules.txt'
    methMolExt='_methylatedMolecules.txt'
    methRatioExt='_methylationRatio.txt'
    runStatsExt='_runStatistics.txt'
    sampleInfoExt='_sampleInfo.txt'

    rscriptPlotSingleFileName='plotSingle.R'
    rHelperFxnsFileName='rFunctions.R'

    #############
    # arg check #
    #############
    # check that paths are valid
    fileDir=utilities.validDir(fileDir)

    os.chdir(fileDir)

    # make sure file dir has '/' on end
    if fileDir[-1]!='/':
        fileDir+='/'

    # check fragInfo file
    if os.path.isfile(fragInfoFile) == False:
        print('ERROR: Your --fragInfo file does not exist')
        print(fragInfoFile)
        print('Exiting...')
        sys.exit()

    totalMolFile=utilities.getFiles(fileDir,totalMolExt)
    totalMolFile=checkInputIsUnique(totalMolFile,totalMolExt,fileDir)

    methMolFile=utilities.getFiles(fileDir,methMolExt)
    methMolFile=checkInputIsUnique(methMolFile,methMolExt,fileDir)

    methRatioFile=utilities.getFiles(fileDir,methRatioExt)
    methRatioFile=checkInputIsUnique(methRatioFile,methRatioExt,fileDir)

    runStatsFile=utilities.getFiles(fileDir,runStatsExt)
    runStatsFile=checkInputIsUnique(runStatsFile,runStatsExt,fileDir)

    runName=totalMolFile.split('_')[0]

    rscriptsPath = '/'.join(inspect.getfile(utilities).split('/')[:-1])
    # path = os.path.abspath(__file__)

    rscriptPlotSinglePath=rscriptsPath+'/'+rscriptPlotSingleFileName
    rHelperScriptPath=rscriptsPath+'/'+rHelperFxnsFileName

    sampleInfoFile=utilities.getFiles(fileDir,sampleInfoExt)

    if len(sampleInfoFile)==0:
        sampleInfoFile='NOSAMPLEINFO'
    elif len(sampleInfoFile)==1:
        sampleInfoFile=sampleInfoFile[0]
    else:
        print('ERROR: it looks like you have more than one sample info file at',fileDir)
        print('Exiting...')
        sys.exit()

    rscriptCmd=['Rscript',rscriptPlotSinglePath,rHelperScriptPath,fileDir,runName,totalMolFile,methMolFile,methRatioFile,runStatsFile,fragInfoFile,sampleInfoFile]

    os.system(' '.join(rscriptCmd))
