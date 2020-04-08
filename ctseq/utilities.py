from datetime import datetime
import sys
import os
import glob
import re

def getDate():
    return(str(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

def validDir(path):
    if not os.path.isdir(path):

        print('\n**ERROR**')
        print('The provided path is not a valid path. Exiting...')
        print(path)
        print('\n')
        sys.exit()
    else:
        print("passdir")

def fileCheck(path,fileExt):
    os.chdir(path)

    myFiles=glob.glob('*'+fileExt)

    if len(myFiles) == 0:
        print('\n**ERROR**')
        print('There are no ('+fileExt+') files at:\n'+path)
        print('Exiting...\n')
        sys.exit()
    else:
        print("passfile")


def getFiles(path,fileExt):
    os.chdir(path)

    myFiles=glob.glob('*'+fileExt)

    return(naturalSort(myFiles))


def naturalSort(mylist):
    #got code from https://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(mylist, key = alphanum_key)
