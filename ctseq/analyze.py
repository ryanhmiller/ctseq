from . import methylref
from . import addumis
from . import align
from . import callmolecules
from . import callmethylation


def run(args):

    # make methylation reference
    # methylref.run(args)

    # add umis
    addumis.run(args)

    # align
    align.run(args)

    # call molecules
    callmolecules.run(args)

    # call methylation
    callmethylation.run(args)
