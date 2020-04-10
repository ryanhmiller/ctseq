from .version import __version__

import argparse
import sys

def run_subcommand(args):
    if args.subcommand=='make_methyl_ref':
        from .methylref import run

    elif args.subcommand=='add_umis':
        print("adding umis")

    elif args.subcommand=='align':
        from .align import run

    elif args.subcommand=='call_molecules':
        from .callmolecules import run


    # run the chosen command
    run(args)




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

    ############
    # add_umis #
    ############
    parser_add_umis = subparsers.add_parser('add_umis', help='properly format and add unique molecular identifiers to your fastq files')
    parser_add_umis.add_argument('-r','--refdir', help='full path to directory where you want to build your reference methylation genome\nmust contain a reference file for your intended targets with extension (.fa)', required=True)
    parser_add_umis.set_defaults(func=run_subcommand)

    ###################
    # make_methyl_ref #
    ###################
    parser_make_methyl_ref = subparsers.add_parser('make_methyl_ref', help='make the necessary reference files to use Bismark to align & call methylation')
    parser_make_methyl_ref.add_argument('-r','--refDir', help='Full path to directory where you want to build your reference methylation genome. Must contain a reference file for your intended targets with extension (.fa)', required=True)
    parser_make_methyl_ref.set_defaults(func=run_subcommand)

    #########
    # align #
    #########
    parser_align = subparsers.add_parser('align', help='trims adapters and aligns reads with Bismark')
    parser_align.add_argument('-r','--refDir', help='Full path to directory where you have already built your methylation reference files', required=True)
    parser_align.add_argument('-f','--fqDir', help='Full path to directory where you have your fastq files', required=True)
    parser_align.add_argument('--fExt', help='unique extension of fastq files containing forward reads', required=True)
    parser_align.add_argument('--rExt', help='unique extension of fastq files containing forward reads', required=True)
    parser_align.set_defaults(func=run_subcommand)

    #########
    # align #
    #########
    parser_align = subparsers.add_parser('call_molecules', help='call molecules from the aligned reads from Bismark')
    parser_align.add_argument('-r','--refDir', help='Full path to directory where you have already built your methylation reference files', required=True)
    parser_align.add_argument('-s','--samDir', help='Full path to directory where your .sam files are located', required=True)
    parser_align.add_argument('-c','--consensus', help='consensus threshold to make consensus methylation call from all the reads with the same UMI (e.g. 0.9)', required=True)
    parser_align.add_argument('-p','--processes', help='number of processes', required=True)
    parser_align.set_defaults(func=run_subcommand)


    #######################################
    # parse args and call proper function #
    #######################################

    # args = parser.parse_args()

    if len(sys.argv) > 1:
        args = parser.parse_args()
        args.func(args)
    else:
        # if no subcommand is specified, print help menu
        parser.print_help()




if __name__ == "__main__":
     main()
