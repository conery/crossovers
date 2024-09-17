#
# Visualizations of filtered blocks
#
# John Conery
# University of Oregon

import logging

import matplotlib.pyplot as plt

from .filters import SNPFilter

# List of plot types.  Imported by the top level app to define the names that
# can be entered on the command line, and used here to define the dispatch
# table.

plot_commands = ['count', 'length', 'location']

def set_params(filter, args):

    dispatch = {
        'chromosomes': filter.set_chromosome,
        'size':        filter.set_size,
        'length':      filter.set_length,
        'coverage':    filter.set_coverage,
        'match':       filter.set_matched,
    }

    params = vars(args)
    for f in dispatch.keys():
        if v := params.get(f):
            dispatch[f](v)


def count_histogram(df, args):
    fig, ax = plt.subplots()
    plt.hist(df.blk_size, bins=10, rwidth=0.8, align='left', range=(1,100), label=args.chromosomes)
    plt.title('Block Size')
    plt.xlabel('Number of SNPs')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    plt.show()

def length_histogram(df, args):
    fig, ax = plt.subplots()
    plt.hist(df.blk_len, bins=10, rwidth=0.8, label=args.chromosomes)
    plt.title('Block Length')
    plt.xlabel('Length (bp)')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    plt.show()

def location_histogram(df, args):   
    fig, ax = plt.subplots()
    plt.hist(df.blk_loc, bins=100, range=(0,1), label=args.chromosomes)
    plt.title('Block Location')
    plt.xlabel('Relative Position in the Chromosome')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    plt.show()


def visualize(args):

    dispatch = {
        'count': count_histogram,
        'length': length_histogram,
        'location':  location_histogram,
    }

    filter = SNPFilter()
    set_params(filter, args)

    logging.info('loading SNP data')
    filter.load_data(args.peaks)

    logging.info('filtering')
    summary = filter.summary()

    logging.info('plotting')
    dispatch[args.command](summary, args)

    logging.info('exit')
