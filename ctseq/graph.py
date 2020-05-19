import os
import inspect
from . import utilities





def run(args):
    fileDir=args.dir
    fragOrderFile=args.molDepthOrder

    totalMolExt='_totalMolecules.txt'
    rscriptName='graph.R'
    #############
    # arg check #
    #############
    # check that paths are valid
    fileDir=utilities.validDir(fileDir)

    # check that refDir/fileDir actually have the necessary files in them
    utilities.fileCheck(fileDir,totalMolExt)
    utilities.fileCheck(fileDir,fragOrderFile)

    # make sure file dir has '/' on end
    if fileDir[-1]!='/':
        fileDir+='/'


    rscriptPath = '/'.join(inspect.getfile(utilities).split('/')[:-1])
    # path = os.path.abspath(__file__)

    rscriptPath+='/'+rscriptName


    rscriptCmd=['Rscript',rscriptPath,fileDir,totalMolExt,fragOrderFile]

    os.system(' '.join(rscriptCmd))
