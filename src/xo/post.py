#
# Postprocessing of filtered blocks
#
# John Conery
# University of Oregon

import logging

import numpy as np
from scipy.signal import find_peaks
from .filters import SNPFilter

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
            print(attr, val)
            setattr(filter, attr, val)

def save_block(df, loc, count):
    start = loc - count
    end = loc
    print(f'SNP {start}..{end-1}')
    nco = df.iloc[start:end]
    print(nco)    

def scan_block(block, args):
    '''
    Scan a block for runs that have satisfy heterozigosity requirements
    and have a minimum number of SNPs.
    '''

    def next_interval():
        nonlocal i, j
        j = i
        while j < len(p):
            # logging.debug(f'{i} {j} {p.iloc[j]}')
            if np.isnan(p.iloc[j]):
                if j - i >= args.min_snps:
                    break
                i = j+1
            j += 1
        logging.debug(f'interval {i}..{j} out of {len(p)}')
        return j - i >= args.min_snps
    
    res = []

    z = block.homozygosity
    if block.iloc[0].background == 'CB4856':
        rr = block.ref_reads
        p = z.where((z >= args.min_z) & (rr >= args.coverage))
    else:
        rr = block.ref_reads
        vr = block.var_reads
        p = z.where((abs(z - 0.5) <= args.delta_z) & (rr >= args.coverage) & (vr >= args.coverage))

    i = j = 0
    while next_interval():
        nco = block.iloc[i:j]
        res.append(nco)
        i = j+1

    if res:
        res.sort(key=lambda b: len(b))
        nco = res[-1]
    else:
        nco = None

    return nco

def postprocess(args):
    '''
    Main function of the `post` script.  

    Arguments:
      args:  command line arguments
    '''

    filter = SNPFilter()
    set_params(filter, args)

    logging.info('loading SNP data')
    filter.load_data(args.peaks)

    logging.info('filtering')
    blocks, summary = filter.apply()

    # summary is a data frame where each row is a block that passed all the
    # filters; use the indexes in that frame to fetch the blocks of SNPs where
    # the data resides
 
    for blk_id, _ in summary.iterrows():
        block = blocks.get_group(blk_id)
        nco = scan_block(block, args)
        if nco is not None:
            print(nco)
        
