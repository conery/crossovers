#
# Visualizations of filtered blocks
#
# John Conery
# University of Oregon

import logging

import matplotlib.pyplot as plt
import matplotlib

import pandas as pd

from .filters import SNPFilter

# List of plot types.  Imported by the top level app to define the names that
# can be entered on the command line, and used here to define the dispatch
# table.

plot_commands = ['count', 'length', 'location']

# Associate command line arguments with filter attributes

filter_params = {
    'chromosomes':  'chromosome',
    'size':         'size_range',
    'length':       'length_range',
    'coverage':     'coverage',
    'match':        'matched',
}

def set_params(filter, args):
    '''
    Scan the command line arguments to find filtering parameters,
    save values in the filter so they can be applied.

    Arguments:
      filter: the SNPFilter object that will do the filtering
      args: command line arguments
    '''
    for arg, attr in filter_params.items():
        if val := vars(args).get(arg):
            setattr(filter, attr, val)

def count_histogram(df, args):
    '''
    Use matplotlib to create and display at histogram of block sizes.

    Arguments:
      df:  the summary data frame created by applying filters
      args:  command line arguments
    '''
    fig, ax = plt.subplots()
    plt.hist(df.blk_size, bins=10, rwidth=0.8, align='left', range=(1,100), label=args.chromosomes)
    plt.title('Block Size')
    plt.xlabel('Number of SNPs')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    plt.show()

def length_histogram(df, args):
    '''
    Use matplotlib to create and display at histogram of block lengths.

    Arguments:
      df:  the summary data frame created by applying filters
      args:  command line arguments
    '''
    fig, ax = plt.subplots()
    plt.hist(df.blk_len, bins=10, rwidth=0.8, label=args.chromosomes)
    plt.title('Block Length')
    plt.xlabel('Length (bp)')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    plt.show()

def location_histogram(df, args):
    '''
    Use matplotlib to create and display at histogram of block locations.

    Arguments:
      df:  the summary data frame created by applying filters
      args:  command line arguments
    '''
    fig, ax = plt.subplots()
    plt.hist(df.blk_loc, bins=100, range=(0,1), label=args.chromosomes)
    plt.title('Block Location')
    plt.xlabel('Relative Position in the Chromosome')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    plt.show()


def visualize(args):
    '''
    Main function of the `vis` script.  Creates a SNPFilter object, saves
    filtering parameters found on the command line, loads and filters the
    SNP data, generates a histogram.

    Arguments:
      args:  command line arguments
    '''
    matplotlib.rcParams.update({'font.size': 12})

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
    _, summary = filter.apply()

    logging.info('plotting')
    dispatch[args.command](summary, args)

    if args.save:
        logging.info(f'writing summary frame to {args.save}')
        summary.to_csv(args.save)

    logging.info('exit')
