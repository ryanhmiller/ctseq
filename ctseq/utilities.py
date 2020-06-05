from datetime import datetime
import sys
import os
import glob
import re
from .__main__ import defaultDir

def getDate():
    return(str(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

def validDir(path):
    if path==defaultDir: # if no path given, return current working directory
        return(os.getcwd())

    elif os.path.isdir(path): # if path is valid, return path
        return(path)

    else: # did not provide valid path
    #if not os.path.isdir(path):
        print('\n**ERROR**')
        print('The provided path is not a valid path. Exiting...')
        print(path)
        print('\n')
        sys.exit()


    # else:
    #     print("passdir")

def fileCheck(path,fileExt):
    os.chdir(path)

    myFiles=glob.glob('*'+fileExt)

    if len(myFiles) == 0:
        print('\n**ERROR**')
        print('There are no ('+fileExt+') files at:\n'+path)
        print('Exiting...\n')
        sys.exit()
    # else:
    #     print("passfile")


def getFiles(path,fileExt):
    os.chdir(path)

    myFiles=glob.glob('*'+fileExt)

    return(naturalSort(myFiles))

def checkInputFileIsUnique(fileList,fileExt,myoutputDir):
    if len(fileList) > 1:
        print('ERROR: it looks like you have more than one *'+fileExt+' file at '+myoutputDir)
        print('Please make sure there is only one of these files at this location')
        print('Exiting...')
        sys.exit()
    elif len(fileList)==0:
        print('ERROR: it looks like you are missing your *'+fileExt+' file at '+myoutputDir)
        print('Exiting...')
        sys.exit()
    else:
        return(fileList[0])

def naturalSort(mylist):
    #got code from https://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(mylist, key = alphanum_key)


def getReferenceFile(refDir):
    os.chdir(refDir)

    refFiles=getFiles(refDir,'.fa')

    refFilePath=refDir

    if refFilePath[-1]!='/':
        refFilePath+='/'

    if len(refFiles) > 1:
        print('\n**ERROR**')
        print('You have more than 1 \'.fa\' reference file at: '+refDir)
        print('\n**Exiting**')
        sys.exit()
    else:
        refFilePath+=refFiles[0]

    return(refFilePath)
