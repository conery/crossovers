#
# Visualizations of filtered blocks
#
# John Conery
# University of Oregon

import pandas as pd
import panel as pn
import numpy as np
import re
import logging

import matplotlib.pyplot as plt

from .gui import SNPFilter

# List of plot types.  Imported by the top level app to define the names that
# can be entered on the command line, and used here to define the dispatch
# table.

plot_commands = ['count', 'length', 'location']

def load_data(args):
    df = pd.read_csv(args.peaks)
    res = []
    for chr_id, chr in df.groupby('chrom_id'):
        if not re.search(args.chromosomes, chr_id):
            continue
        for blk_id, blk in chr.groupby('blk_id'):
            res.append(blk)
    return res

def set_params(filters, args):
    params = vars(args)
    for f in filters:
        if p := params.get(f):
            logging.info(f'set {f} to {p}')
            if isinstance(p, list):
                filters[f].value = tuple(p)
            else:
                filters[f].value = p

def count_histogram(blks, args):
    data = [len(b) for b in blks]
    fig, ax = plt.subplots()
    plt.hist(data, bins=10, rwidth=0.8, align='left', range=(1,100), label=args.chromosomes)
    plt.title('Block Size')
    plt.xlabel('Number of SNPs')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    plt.show()

def length_histogram(blks, args):
    data = [b.iloc[-1].position - b.iloc[0].position for b in blks if len(b) > 0]
    fig, ax = plt.subplots()
    plt.hist(data, bins=10, rwidth=0.8, label=args.chromosomes)
    plt.title('Block Length')
    plt.xlabel('Length (bp)')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    plt.show()

# Lengths of chromosomes, used by location histogram

chr_length = {
    1: 15114068,
    2: 15311845,
    3: 13819453,
    4: 17493838,
    5: 20953657,
    6: 17739129,
}

def relative_loc(i, n):
    '''
    Return the relative location of base i on chromosome n
    '''
    return i / chr_length[n]

def location_histogram(blks, args):
    chromosomes = set()
    data = []
    for b in blks:
        if len(b) == 0:
            continue
        chromosomes.add(b.iloc[0].chrom_id)
        chr_id = b.iloc[0].chromosome
        loc = b.position.apply(relative_loc, args=(chr_id,)).mean()
        data.append(loc)
    fig, ax = plt.subplots()
    plt.hist(data, bins=100, range=(0,1), label=args.chromosomes)
    plt.title('Block Location')
    plt.xlabel('Relative Position in the Chromosome')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    plt.show()

dispatch = {
    'count': count_histogram,
    'length': length_histogram,
    'location':  location_histogram,
}

def visualize(args):
    logging.info('making filters')
    filter = SNPFilter()
    set_params(filter.widget_map(), args)
    logging.info('loading data')
    blocks = [filter.apply(b) for b in load_data(args)]
    logging.info(f'read {len(blocks)} blocks')
    logging.info('creating plot')
    dispatch[args.command](blocks, args)

