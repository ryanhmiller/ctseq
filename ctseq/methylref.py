import os
import sys
import subprocess
from . import utilities


def run(args):
    refDir=args.refDir

    # print(__name__)
    ##############################################
    # arg check - make sure this is a valid path #
    ##############################################
    utilities.validDir(refDir)

    #########################
    # continue with command #
    #########################

    # make methylation reference genome with bismark
    try:
        os.chdir(refDir)

        print('\n**************')
        print('Making methylation reference with Bismark '+utilities.getDate())
        print('**************\n')

        makeMethylRefCmd=['bismark_genome_preparation','--bowtie2',refDir]
        externalProcess = subprocess.Popen(makeMethylRefCmd, stdout=subprocess.PIPE)
        exitCode = externalProcess.wait()
        error = externalProcess.stderr


        if exitCode==0:
            print('\n**************')
            print('Done making methylation reference with Bismark '+utilities.getDate())
            print('**************\n')
        else:
            raise Exception

    except:
        e = sys.exc_info()[0]
        print('\n**ERROR**\n(make_methyl_ref)\n')
        print('There was an error making the methylation reference using bismark.')
        print('Bismark exit code:')
        print(str(exitCode)+'\n*********')
