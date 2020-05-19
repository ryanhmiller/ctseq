from .version import __version__

import argparse
import sys

def run_subcommand(args):
    if args.subcommand=='make_methyl_ref':
        from .methylref import run

    elif args.subcommand=='add_umis':
        from .addumis import run

    elif args.subcommand=='align':
        from .align import run

    elif args.subcommand=='call_molecules':
        from .callmolecules import run

    elif args.subcommand=='call_methylation':
        from .callmethylation import run

    elif args.subcommand=='analyze':
        from .analyze import run

    elif args.subcommand=='graph':
        from .graph import run

    # run the chosen command
    run(args)


defaultDir='CurrDir'

def main():

    ####################
    # top level parser #
    ####################
    parser = argparse.ArgumentParser(prog='ctseq',description='ctseq analyzes patch pcr data (only methylated patch pcr data at this point)')
    parser.add_argument('-v', '--version', help="show the installed ctseq version",
                        action="version",
                        version="%(prog)s installed version: " + str(__version__))

    subparsers = parser.add_subparsers(title='[Subcommands]',
                                        description='valid subcommands',
                                        help='additional help',
                                        dest='subcommand')

    ###########################
    # individual tool parsers #
    ###########################

    ###################
    # make_methyl_ref #
    ###################
    parser_make_methyl_ref = subparsers.add_parser('make_methyl_ref', help='make the necessary reference files to use Bismark to align & call methylation')
    parser_make_methyl_ref.add_argument('-r','--refDir', help='Path to the directory where you want to build your reference methylation genome. Must contain a reference file for your intended targets with extension (.fa). If this flag is not used, it is assumed that your current working directory contains your reference file and that you want to build your reference in your current working directory.', default=defaultDir)
    parser_make_methyl_ref.set_defaults(func=run_subcommand)

    ############
    # add_umis #
    ############
    parser_add_umis = subparsers.add_parser('add_umis', help='properly format and add unique molecular identifiers to your fastq files')
    parser_add_umis.add_argument('--umiType', choices=['separate', 'inline'], help='Choose \'separate\' if the UMIs for the reads are contained in a separate fastq file where the line after the read name is the UMI. Choose \'inline\' if the UMIs are already included in the forward/reverse read fastq files in the following format: \'@M01806:488:000000000-J36GT:1:1101:15963:1363:GTAGGTAAAGTG 1:N:0:CGAGTAAT\' where \'GTAGGTAAAGTG\' is the UMI', required=True)
    parser_add_umis.add_argument('-l','--umiLength', help='Length of UMI sequence, e.g. 12 (required)', type=int, required=True)
    parser_add_umis.add_argument('-d','--dir', help='Path to directory containing fastq files; forward/reverse reads and umi files. If no \'--dir\' is specified, ctseq will look in your current directory.', default=defaultDir)
    parser_add_umis.add_argument('--forwardExt', help='Unique extension of fastq files containing FORWARD reads. Make sure to include \'.gz\' if your files are compressed (required)', required=True)
    parser_add_umis.add_argument('--reverseExt', help='Unique extension of fastq files containing REVERSE reads. Make sure to include \'.gz\' if your files are compressed (required)', required=True)
    parser_add_umis.add_argument('--umiExt', help='Unique extension of fastq files containing the UMIs (This flag is REQUIRED if UMIs are contained in separate fastq file). Make sure to include \'.gz\' if your files are compressed.', default='NOTSPECIFIED')
    parser_add_umis.set_defaults(func=run_subcommand)

    #########
    # align #
    #########
    parser_align = subparsers.add_parser('align', help='trims adapters and aligns reads with Bismark')
    parser_align.add_argument('-r','--refDir', help='Full path to directory where you have already built your methylation reference files. If no \'--refDir\' is specified, ctseq will look in your current directory.', default=defaultDir)
    parser_align.add_argument('-d','--dir', help='Path to directory where you have your fastq files. If no \'--dir\' is specified, ctseq will look in your current directory.', default=defaultDir)
    parser_align.add_argument('--forwardAdapter', help='adapter sequence to remove from FORWARD reads (default=AGTGTGGGAGGGTAGTTGGTGTT)', default='AGTGTGGGAGGGTAGTTGGTGTT')
    parser_align.add_argument('--reverseAdapter', help='adapter sequence to remove from REVERSE reads (default=ACTCCCCACCTTCCTCATTCTCTAAGACGGTGT)', default='ACTCCCCACCTTCCTCATTCTCTAAGACGGTGT')
    parser_align.add_argument('-c','--cutadaptCores', help='number of cores to use with Cutadapt. Default=1. Highly recommended to run with more than 1 core, try starting with 18 cores', default=1, type=int)
    parser_align.add_argument('-b','--bismarkCores', help='number of cores to use to align with Bismark. Default=1. Highly recommended to run with more than 1 core, try starting with 6 cores', default=1, type=int)
    parser_align.add_argument('--readsPerFile', help='number of reads to analyze per fastq file (should only adjust this if you think you are too big of a file through bismark). Default=5000000 (5 million)', default=5000000, type=int)
    parser_align.set_defaults(func=run_subcommand)

    ##################
    # call_molecules #
    ##################
    parser_callmolecules = subparsers.add_parser('call_molecules', help='call molecules from the aligned reads from Bismark')
    parser_callmolecules.add_argument('-r','--refDir', help='Full path to directory where you have already built your methylation reference files. If no \'--refDir\' is specified, ctseq will look in your current directory.', default=defaultDir)
    parser_callmolecules.add_argument('-d','--dir', help='Full path to directory where your .sam files are located. If no \'--dir\' is specified, ctseq will look in your current directory.', default=defaultDir)
    parser_callmolecules.add_argument('-c','--consensus', help='consensus threshold to make consensus methylation call from all the reads with the same UMI (default=0.9)', default=0.9, type=float)
    parser_callmolecules.add_argument('-p','--processes', help='number of processes (default=1; default settings could take a long time to run)', default=1, type=int)
    parser_callmolecules.add_argument('-u','--umiThreshold', help='UMIs with this edit distance will be collapsed together, default=0 (don\'t collapse)', default=0, type=int)
    parser_callmolecules.add_argument('-a','--umiCollapseAlg', help='algorithm used to collapse UMIs, options: default=directional', default='directional')
    parser_callmolecules.set_defaults(func=run_subcommand)

    ####################
    # call_methylation #
    ####################
    parser_callmethylation = subparsers.add_parser('call_methylation', help='call methylation from the \'*allMolecules.txt\' file')
    parser_callmethylation.add_argument('-r','--refDir', help='Full path to directory where you have already built your methylation reference files. If no \'--refDir\' is specified, ctseq will look in your current directory.', default=defaultDir)
    parser_callmethylation.add_argument('-d','--dir', help='Full path to directory where your \'*allMolecules.txt\' files are located. If no \'--dir\' is specified, ctseq will look in your current directory.', default=defaultDir)
    parser_callmethylation.add_argument('-n','--nameRun', help='number of reads needed to be counted as a unique molecule (required)', required=True)
    parser_callmethylation.add_argument('-p','--processes', help='number of processes (default=1)', default=1, type=int)
    parser_callmethylation.add_argument('-c','--cisCG', help='cis-CG threshold to determine if a molecule is methylated (default=0.75)', default=0.75, type=float)
    parser_callmethylation.add_argument('-m','--moleculeThreshold', help='number of reads needed to be counted as a unique molecule (default=5)', default=5, type=int)
    parser_callmethylation.set_defaults(func=run_subcommand)

    ########################################################
    # analyze - combines add_umis through call_methylation #
    ########################################################
    parser_analyze = subparsers.add_parser('analyze', help='runs all the steps in the pipeline')
    parser_analyze.add_argument('--umiType', choices=['separate', 'inline'], help='Choose \'separate\' if the UMIs for the reads are contained in a separate fastq file where the line after the read name is the UMI. Choose \'inline\' if the UMIs are already included in the forward/reverse read fastq files in the following format: \'@M01806:488:000000000-J36GT:1:1101:15963:1363:GTAGGTAAAGTG 1:N:0:CGAGTAAT\' where \'GTAGGTAAAGTG\' is the UMI', required=True)
    parser_analyze.add_argument('--umiLength', help='Length of UMI sequence, e.g. 12 (required)', type=int, required=True)
    parser_analyze.add_argument('-d','--dir', help='Path to directory where you have your fastq files. If no \'--dir\' is specified, ctseq will look in your current directory.', default=defaultDir)
    parser_analyze.add_argument('-r','--refDir', help='Full path to directory where you have already built your methylation reference files. If no \'--refDir\' is specified, ctseq will look in your current directory.', default=defaultDir)
    parser_analyze.add_argument('--forwardExt', help='Unique extension of fastq files containing FORWARD reads. Make sure to include \'.gz\' if your files are compressed (required)', required=True)
    parser_analyze.add_argument('--reverseExt', help='Unique extension of fastq files containing REVERSE reads. Make sure to include \'.gz\' if your files are compressed (required)', required=True)
    parser_analyze.add_argument('--umiExt', help='Unique extension of fastq files containing the UMIs (This flag is REQUIRED if UMIs are contained in separate fastq file). Make sure to include \'.gz\' if your files are compressed.', default='NOTSPECIFIED')
    parser_analyze.add_argument('--forwardAdapter', help='adapter sequence to remove from FORWARD reads (default=AGTGTGGGAGGGTAGTTGGTGTT)', default='AGTGTGGGAGGGTAGTTGGTGTT')
    parser_analyze.add_argument('--reverseAdapter', help='adapter sequence to remove from REVERSE reads (default=ACTCCCCACCTTCCTCATTCTCTAAGACGGTGT)', default='ACTCCCCACCTTCCTCATTCTCTAAGACGGTGT')
    parser_analyze.add_argument('--cutadaptCores', help='number of cores to use with Cutadapt. Default=1. Highly recommended to run with more than 1 core, try starting with 18 cores', default=1, type=int)
    parser_analyze.add_argument('--bismarkCores', help='number of cores to use to align with Bismark. Default=1. Highly recommended to run with more than 1 core, try starting with 6 cores', default=1, type=int)
    parser_analyze.add_argument('--readsPerFile', help='number of reads to analyze per fastq file (should only adjust this if you think you are too big of a file through bismark). Default=5000000 (5 million)', default=5000000, type=int)


    parser_analyze.add_argument('--consensus', help='consensus threshold to make consensus methylation call from all the reads with the same UMI (default=0.9)', default=0.9, type=float)
    parser_analyze.add_argument('-p','--processes', help='number of processes (default=1; default settings could take a long time to run)', default=1, type=int)
    parser_analyze.add_argument('--umiThreshold', help='UMIs with this edit distance will be collapsed together, default=0 (don\'t collapse)', default=0, type=int)
    parser_analyze.add_argument('--umiCollapseAlg', help='algorithm used to collapse UMIs, options: default=directional', default='directional')

    parser_analyze.add_argument('-n','--nameRun', help='number of reads needed to be counted as a unique molecule (required)', required=True)
    parser_analyze.add_argument('--cisCG', help='cis-CG threshold to determine if a molecule is methylated (default=0.75)', default=0.75, type=float)
    parser_analyze.add_argument('--moleculeThreshold', help='number of reads needed to be counted as a unique molecule (default=5)', default=5, type=int)
    parser_analyze.set_defaults(func=run_subcommand)


    #########
    # graph #
    #########
    parser_graph = subparsers.add_parser('graph', help='graph output from ctseq')
    parser_graph.add_argument('--dir', help='Path to directory where you have your graph input files. If no \'--dir\' is specified, ctseq will look in your current directory.', default=defaultDir)
    parser_graph.add_argument('--molDepthOrder', help='Name of file containing order of your fragments to be displayed on the x-axis of the molecule depth plot', required=True)

    parser_graph.set_defaults(func=run_subcommand)



    #######################################
    # parse args and call proper function #
    #######################################

    # args = parser.parse_args()

    if sys.argv[1]=='graph': # we can run 'graph' without any args
        args = parser.parse_args()
        args.func(args)
    elif sys.argv[1]!='graph' and len(sys.argv) > 1: # can't run any other subcommands without any args
        args = parser.parse_args()
        args.func(args)
    else:
        # if no subcommand is specified, print help menu
        parser.print_help()




if __name__ == "__main__":
     main()
