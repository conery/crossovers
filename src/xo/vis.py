#
# Visualizations of filtered blocks
#
# John Conery
# University of Oregon

import concurrent.futures
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
    logging.info('loading data')
    df = pd.read_csv(args.peaks)
    logging.info(f'frame has {len(df)} rows')
    res = []
    groups = df.groupby('chrom_id')
    logging.info('chromosomes grouped')
    for chr_id, chr in groups:
        if not re.search(args.chromosomes, chr_id):
            continue
        for blk_id, blk in chr.groupby('blk_id'):
            res.append(blk)
    logging.info(f'{len(res)} blocks')
    return res

def set_params(filters, args):
    logging.info('making filters')
    params = vars(args)
    for f in filters:
        if p := params.get(f):
            logging.info(f'set {f} to {p}')
            if isinstance(p, list):
                filters[f].value = tuple(p)
            else:
                filters[f].value = p

def count_histogram(blks, args):
    logging.info(f'starting count_histogram')
    data = [len(b) for b in blks]
    logging.info(f'data array has {len(data)} points')
    fig, ax = plt.subplots()
    plt.hist(data, bins=10, rwidth=0.8, align='left', range=(1,100), label=args.chromosomes)
    plt.title('Block Size')
    plt.xlabel('Number of SNPs')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    logging.info('plot ready')
    plt.show()

def length_histogram(blks, args):
    logging.info(f'starting length_histogram')
    data = [b.iloc[-1].position - b.iloc[0].position for b in blks if len(b) > 0]
    logging.info(f'data array has {len(data)} points')
    fig, ax = plt.subplots()
    plt.hist(data, bins=10, rwidth=0.8, label=args.chromosomes)
    plt.title('Block Length')
    plt.xlabel('Length (bp)')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    logging.info('plot ready')
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
    # chromosomes = set()
    data = []
    for b in blks:
        if len(b) == 0:
            continue
        # chromosomes.add(b.iloc[0].chrom_id)
        chr_id = b.iloc[0].chromosome
        loc = b.position.apply(relative_loc, args=(chr_id,)).mean()
        data.append(loc)
    logging.info(f'data array has {len(data)} points')
    fig, ax = plt.subplots()
    plt.hist(data, bins=100, range=(0,1), label=args.chromosomes)
    plt.title('Block Location')
    plt.xlabel('Relative Position in the Chromosome')
    plt.ylabel('Number of Blocks')
    plt.legend(handlelength=0)
    logging.info('plot ready')
    plt.show()

# sequential version
def apply_filters(filter, blocks):
    return [filter.apply(b) for b in blocks] 

# # First parallel version -- send lots of small blocks, very high communication
# # overhead, much slower than sequential version.

# def apply_filters(filter, blks):
#     with concurrent.futures.ProcessPoolExecutor(8) as executor:
#         return executor.map(filter.apply, blks)

# # Second parallel version -- split list of blocks unto separate chunks, send
# # on chunk to each process.

# def task(filter, blks, start, end):
#     return [filter.apply(blks[i]) for i in range(start,end)]

# def apply_filters(filter, blks):

#     cap = len(blks)
#     chunksize = 10000
#     res = []

#     with concurrent.futures.ProcessPoolExecutor() as executor:
#         try:
#             futures = [executor.submit(task, filter, blks, b*chunksize, min(cap,(b+1)*chunksize)) for b in range(cap//chunksize+1)]
#             for future in concurrent.futures.as_completed(futures):
#                 res.extend(future.result())
#         except Exception as err:
#             print(err)

#     return res


dispatch = {
    'count': count_histogram,
    'length': length_histogram,
    'location':  location_histogram,
}

def visualize(args):
    filter = SNPFilter()
    set_params(filter.widget_map(), args)
    blocks = load_data(args)
    logging.info('start filtering')
    filtered = apply_filters(filter, blocks)
    logging.info('calling plot command')
    dispatch[args.command](filtered, args)
    logging.info('exit')
