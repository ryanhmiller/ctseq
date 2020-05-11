from . import methylref
from . import addumis


def run(args):

    # make methylation reference
    methylref.run(args)

    # add umis to
