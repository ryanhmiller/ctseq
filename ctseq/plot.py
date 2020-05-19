import os
import inspect
import sys
from . import utilities



def run(args):
    fileDir=args.dir
    fragOrderFile=args.molDepthOrder

    totalMolExt='_totalMolecules.txt'
    rscriptName='plot.R'
    #############
    # arg check #
    #############
    # check that paths are valid
    fileDir=utilities.validDir(fileDir)

    # make sure file dir has '/' on end
    if fileDir[-1]!='/':
        fileDir+='/'

    # check totalmolecules file
    if os.path.isfile(fileDir+fragOrderFile) == False:
        print('ERROR: Your --molDepthOrder order file does not exist')
        print(fragOrderFile)
        print('Exiting...')
        sys.exit()

    totalMolFile=utilities.getFiles(fileDir,totalMolExt)

    if len(totalMolFile) > 1:
        print('ERROR: it looks like you have more than one *'+totalMolExt+' file at '+fileDir)
        print('Please make sure there is only one of these files at this location')
        print('Exiting...')
        sys.exit()
    else:
        totalMolFile=totalMolFile[0]


    rscriptPath = '/'.join(inspect.getfile(utilities).split('/')[:-1])
    # path = os.path.abspath(__file__)

    rscriptPath+='/'+rscriptName

    rscriptCmd=['Rscript',rscriptPath,fileDir,totalMolFile,fragOrderFile]

    os.system(' '.join(rscriptCmd))
