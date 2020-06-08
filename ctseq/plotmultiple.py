import os
import inspect
import sys
from . import utilities


def checkFile(filePath,fileType):
    if os.path.isfile(filePath) == False:
        print('ERROR: Your',fileType,'file does not exist')
        print(filePath)
        print('Exiting...')
        sys.exit()


def run(args):
    fileDir=args.dir
    fragInfoFile=args.fragInfo
    runName=args.name

    totalMolExt='_totalMolecules.txt'
    methMolExt='_methylatedMolecules.txt'
    methRatioExt='_methylationRatio.txt'
    runStatsExt='_runStatistics.txt'
    sampleInfoExt='_sampleInfo.txt'
    directoriesFileExt='_directories.txt'

    rscriptPlotMultipleFileName='plotMultiple.R'
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
    checkFile(fragInfoFile,'fragInfo')

    # check directories file
    directoriesFileName=utilities.getFiles(fileDir,directoriesFileExt)
    directoriesFileName=utilities.checkInputFileIsUnique(directoriesFileName,directoriesFileExt,fileDir)

    # loop through paths in dir and check make sure necessary files are present
    with open(directoriesFileName,'r') as dirFile:
        for line in dirFile:
            currDir=line.strip('\n')
            currDir=utilities.validDir(currDir)

            totalMolFile=utilities.getFiles(currDir,totalMolExt)
            totalMolFile=utilities.checkInputFileIsUnique(totalMolFile,totalMolExt,currDir)

            methMolFile=utilities.getFiles(currDir,methMolExt)
            methMolFile=utilities.checkInputFileIsUnique(methMolFile,methMolExt,currDir)

            methRatioFile=utilities.getFiles(currDir,methRatioExt)
            methRatioFile=utilities.checkInputFileIsUnique(methRatioFile,methRatioExt,currDir)

            runStatsFile=utilities.getFiles(currDir,runStatsExt)
            runStatsFile=utilities.checkInputFileIsUnique(runStatsFile,runStatsExt,currDir)


    rscriptsPath = '/'.join(inspect.getfile(utilities).split('/')[:-1])
    # path = os.path.abspath(__file__)

    rscriptPlotMultiplePath=rscriptsPath+'/'+rscriptPlotMultipleFileName
    rHelperScriptPath=rscriptsPath+'/'+rHelperFxnsFileName

    rscriptCmd=['Rscript',rscriptPlotMultiplePath,runName,totalMolExt,runStatsExt,methRatioExt,methMolExt,sampleInfoExt,fragInfoFile,fileDir,directoriesFileName,rHelperScriptPath]

    os.system(' '.join(rscriptCmd))
