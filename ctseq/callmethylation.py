import os
from . import utilities




def run(args):
    refDir=args.refDir
    fileDir=args.dir
    processes=int(args.processes)
    cisCGcutoff=float(args.cisCG)
    umiThreshold=int(args.umiThreshold)
    umiCollapseAlg=args.umiCollapseAlg
    fileExt='_allMolecules.txt'

    blacklistUMI=''

    #############
    # arg check #
    #############
    # check that paths are valid
    utilities.validDir(refDir)
    utilities.validDir(fileDir)

    # check that refDir/fileDir actually have the necessary files in them
    utilities.fileCheck(refDir,".fa")
    utilities.fileCheck(fileDir,fileExt)

    # make sure sam file dir has '/' on end
    if fileDir[-1]!='/':
        fileDir+='/'
    ####################
    ### SET UP STUFF ###
    ####################
    # get reference .fa file
    refFileName=utilities.getReferenceFile(refDir)


    # get list of reference fragment names
    refFrags=[]

    with open(refFileName, "r") as refFile:
        lines=refFile.readlines()
        for i in range(0,len(lines),2):
            fragName=lines[i].strip("\n")[1:]
            seq=lines[i+1].strip("\n")
            seq=seq.upper()

            refFrags.append(fragName)


    os.chdir(fileDir)

    listOfFiles=utilities.getFiles(fileDir,fileExt)

    
    # ####################
    mySamples=[] # builds out list of all samples in directory

    for fileName in listOfFiles:
        with open(fileName,"r") as infile:
            infile.readline()

            sampleName=fileName.split("_")[0]
            mySamples.append(sampleName)

            locusDict={} # fragName:locusObject

            methylationDict={} # frag:%methylated
            readDepthDict={} # frag:readDepth
            umiDict={} # frag:umiDepth
            methMolDict={} # frag:#methMol

            for line in infile:

                line=line.strip("\n").split("\t")

                umi=line[2]
                zString=line[3]                

                if zString!="NA" and umi != blacklistUMI:
                    locus=line[1]       
                    totalCGs=int(line[5])
                    methRatio=float(line[6])
                    #molecules=int(line[5])
                    #umiDepth=int(line[6])
                    readDepth=int(line[7])
                    totalReadsPerSample[sampleName]+=readDepth


                    ##### NEED TO COUNT UMI DEPTH + READ DEPTH #####

                    # if methRatio >= cgMethylationCutoff:
                    #
                    #     print locus, str(molecules),str(methRatio), str(cgMethylationCutoff)


                    if locus not in locusDict:
                        newLocus=Locus()
                        newLocus.name=locus
                        newLocus.totalCGs=totalCGs
                        #newLocus.umiDepth=umiDepth
                        #newLocus.totalReads=readDepth
                        locusDict[locus]=newLocus

                    locusDict[locus].totalReads+=readDepth








