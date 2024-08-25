#! /usr/bin/env python3

# Command line application to extract segments from the SNP table.  Uses
# SciPy's find_peaks function to look for regions of the table that change
# HMM state.

# John Conery
# University of Oregon

import argparse
import pandas as pd
from scipy.signal import find_peaks

from rich.console import Console

def init_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--snps', metavar='F', default='BSP_TIGER.marker_dataframe.pickle.gzip', help='TIGER marker file')
    parser.add_argument('--output', metavar='F', default='peaks.pickle.gzip', help='output file')
    parser.add_argument('--sample', metavar='F', help='write a part of the input to file F and exit')
    parser.add_argument('--size', metavar='N', type=int, default=100000, help='number of records to write in sample')
    parser.add_argument('--compression', metavar='X', default='gzip', help='method used to compress input file')
    return parser.parse_args()

def extract_blocks(chromosome):
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
        df = chromosome.iloc[blk_start:blk_end+1]
        blocks.append(df.assign(blk_id=i))
        # print(chromosome.iloc[0].chrom_id, i, px[i], prop['prominences'][i], blk_start, blk_end)
    return pd.concat(blocks) if blocks else None

if __name__ == '__main__':
    args = init_cli()
    console = Console()
    with console.status(f'Processing SNPs', spinner='aesthetic') as status:
        console.log(f'Reading {args.snps}')
        snps = pd.read_pickle(args.snps, compression=args.compression)
        console.log(f'read {len(snps)} SNPs')
        if args.sample:
            status.update(status='Saving sample')
            snps[0:args.size].to_pickle(args.sample)
            console.log(f'Saved {args.size} records')
        else:
            result = []
            for cname, sf in snps.groupby('chrom_id'):
                df = extract_blocks(sf)
                if df is None:
                    console.log(f'[red] no blocks in {cname}')
                else:
                    console.log(f'{cname}: {len(sf)} SNPs {len(df)} in blocks')
                    result.append(df)
            final = pd.concat(result)
            console.log(f'Writing to {args.output}')
            final.to_pickle(args.output)
            console.log(f'Wrote {len(final)} records')

