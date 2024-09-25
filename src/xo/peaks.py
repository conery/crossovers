#
# Command line application to extract segments from the SNP table.  Uses
# SciPy's find_peaks function to look for regions of the table that change
# HMM state.
#
# Run the application using the top level xo script:
#
#  $ xo peaks OPTS
#

# John Conery
# University of Oregon

import pandas as pd
from scipy.signal import find_peaks

from rich.console import Console

def extract_blocks(chromosome, max_block_size):
    '''
    Use `find_peaks` from the SciPy signal processing library to look for
    sequences of SNPs.  Sequences that "stand out" are collected into a block,
    represented by a data frame with a new column appended to hold the block ID.

    Arguments:
      chromosome:  a data frame where each row describes a SNP
      max_block_size:  the maximum number of SNPs to include in a block

    Returns:
      a data frame containing all the SNPs in blocks.
    '''
    signal = ((chromosome.hmm_state1 == 'CB4856').cumsum() - (chromosome.hmm_state1 == 'N2').cumsum()).to_numpy()
    px, prop = find_peaks(signal, prominence=1)
    blocks = []
    for i in range(len(px)):
        if prop['left_bases'][i] == px[i] - prop['prominences'][i]:
            blk_start = prop['left_bases'][i] + 1
            blk_end = px[i]
        else:
            blk_start = px[i]+1
            blk_end = prop['right_bases'][i]
        if blk_end - blk_start > max_block_size:
            continue
        df = chromosome.iloc[blk_start:blk_end+1]
        blocks.append(df.assign(blk_id=i))
    return pd.concat(blocks) if blocks else None

def peak_finder(args):
    '''
    Top level function.  Reads the SNP data, groups it by chromosome, and
    calls `find_peaks` for each chromosome.  The results are collected in
    a data frame and written to a CSV file.

    Arguments:
      args:  command line arguments from `argparse`
    '''
    console = Console()
    with console.status(f'Processing SNPs', spinner='aesthetic') as status:
        console.log(f'Reading {args.snps}')
        snps = pd.read_pickle(args.snps, compression='gzip')
        console.log(f'read {len(snps)} SNPs')
        result = []
        for cname, sf in snps.groupby('chrom_id'):
            df = extract_blocks(sf, args.max_snps)
            if df is None:
                console.log(f'[red] no blocks in {cname}')
            else:
                console.log(f'{cname}: {len(sf)} SNPs {len(df)} in blocks')
                result.append(df)
        final = pd.concat(result)
        console.log(f'Writing to {args.output}')
        final.to_csv(args.output)
        console.log(f'Wrote {len(final)} records')

