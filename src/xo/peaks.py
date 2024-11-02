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

def majority_background(chr):
    '''
    Figure out whether a chromosome is mostly N2 or CB4856 to use as the
    default background (when a chromosome does not have a crossover).

    Arguments:
      chr: dataframe with SNPs in a chromosome

    Returns:
      'N2' or 'CB4856'
    '''
    grps = chr.groupby('hmm_state1')
    n2_group = grps.groups.get('N2')
    cb_group = grps.groups.get('CB4856')
    n2_count = len(n2_group) if 'N2' in grps else 0
    cb_count = len(cb_group) if 'CB4856' in grps else 0
    return 'N2' if n2_count > cb_count else 'CB4856'

# Chromosome lengths

chr_length = {
    1: 15114068,
    2: 15311845,
    3: 13819453,
    4: 17493838,
    5: 20953657,
    6: 17739129,
}

def add_background(cf, df):
    '''
    Use crossover locations to add a background location to a frame.

    Arguments:
      cf:  a data frame with crossover locations
      df:  the frame to update

    Returns:
      a copy of df with a new background column added
    '''
    # gen = df.iloc[0].hmm_state1
    gen = majority_background(df)
    # print('default', gen)
    col = pd.Series(gen, index=df.index)
    if cf is not None:
        left = 0
        right = chr_length[df.iloc[0].chromosome]
        for _, xo in cf.iterrows():
            mid = xo.start
            gen = 'CB4856' if xo.upstream_CB4856_purity > xo.downstream_CB4856_purity else 'N2'
            col[(df.position >= left) & (df.position <= mid)] = gen
            # print(f'{left}..{mid} = {gen}')
            left = mid
        gen = 'N2' if gen == 'CB4856' else 'CB4856'
        col[(df.position > mid) & (df.position < right)] = gen
        # print(f'{mid}..{right} = {gen}')
    df['background'] = col
    return df


def peak_finder(args):
    '''
    Top level function for the `peaks` command.
    Reads the SNP data, groups it by chromosome, and
    calls `extract_blocks` for each chromosome.  The results are collected in
    a data frame and written to a CSV file.

    Arguments:
      args:  command line arguments from `argparse`
    '''
    console = Console()
    with console.status(f'Processing SNPs', spinner='aesthetic') as status:
        console.log(f'Reading {args.snps}')
        snps = pd.read_pickle(args.snps, compression='gzip').groupby('chrom_id')
        console.log(f'read {len(snps)} SNP groups')
        xo = pd.read_pickle(args.crossovers, compression='gzip').groupby('chrom_id')
        console.log(f'read {len(xo)} crossover groups')
        result = []
        for cname, sf in snps:
            df = extract_blocks(sf, args.max_snps)
            if df is None:
                console.log(f'[red] no blocks in {cname}')
            else:
                console.log(f'{cname}: {len(sf)} SNPs {len(df)} in blocks')
                cf = xo.get_group(cname) if cname in xo.groups else None
                df = add_background(cf,df)
                result.append(df)
        final = pd.concat(result)
        console.log(f'Writing to {args.output}')
        final.to_csv(args.output)
        console.log(f'Wrote {len(final)} records')

