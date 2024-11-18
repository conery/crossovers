#! /usr/bin/env python3

# Top level application for the crossovers project.  
#
# John Conery
# University of Oregon
# 
# To run the application type
#
#  $ xo CMND [OPTS]
#
# to run one of the commands -- peak finder, viewer, or filter.

import argparse
import logging
import os
import sys

import pandas as pd

from rich.console import Console
from rich.logging import RichHandler

from .gui import start_app
from .peaks import extract_blocks, add_background, peak_results
from .filters import SNPFilter

# Top level functions, called from main

def peak_finder(args):
    '''
    Top level function for the `peaks` command.  Reads the SNP data, groups 
    it by chromosome, and calls `extract_blocks` for each chromosome.  
    The results are collected in a data frame and written to a CSV file.

    Arguments:
      args:  command line arguments from `argparse`
    '''

    def _find_peaks():
        logging.info(f'Reading {args.snps}')
        snps = pd.read_pickle(args.snps, compression='gzip').groupby('chrom_id')
        logging.info(f'read {len(snps)} SNP groups')
        xo = pd.read_pickle(args.crossovers, compression='gzip').groupby('chrom_id')
        logging.info(f'read {len(xo)} crossover groups')
        result = []
        for cname, sf in snps:
            df = extract_blocks(sf, args.max_snps)
            if df is None:
                logging.info(f'[red]no blocks in {cname}')
            else:
                logging.info(f'{cname}: {len(sf)} SNPs {len(df)} in blocks')
                cf = xo.get_group(cname) if cname in xo.groups else None
                df = add_background(cf,df)
                result.append(df)
        final = pd.concat(result)
        logging.info(f'Writing to {args.output}')
        final.to_csv(args.output)
        logging.info(f'Wrote {len(final)} records')
        return len(snps), final

    logging.debug(f'peaks {vars(args)}')

    # If the log level is set to QUIET use a status spinner to show we're making progress

    if args.log == 'quiet':
        console = Console()
        with console.status(f'Processing SNPs', spinner='aesthetic') as status:
            res = _find_peaks()
    else:
        res = _find_peaks()

    peak_results(res)


# Chromosome lengths

chr_length = {
    1: 15114068,
    2: 15311845,
    3: 13819453,
    4: 17493838,
    5: 20953657,
    6: 17739129,
}

def filter_blocks(args):
    '''
    Top level function for the `filter` command.  Reads the peaks file
    (adding new columns for location and homozygosity), passes it to
    the filter, saves the result.

    Arguments:
      args:  command line arguments from `argparse`
    '''
    logging.debug(f'filter {vars(args)}')

    filter = SNPFilter(args)

    peaks = pd.read_csv(args.peaks).rename(columns={'Unnamed: 0': 'SNP'})
    peaks['chr_length'] = peaks.chromosome.map(lambda n: chr_length[n])
    peaks['location'] = peaks.position / peaks.chr_length
    peaks['homozygosity'] = peaks.ref_reads / (peaks.ref_reads + peaks.var_reads)

    res, _ = filter.apply(peaks)
    res.to_csv(args.output)


def postprocess(args):
    logging.debug(f'post {vars(args)}')

def init_cli():
    """
    Use argparse to create the command line API.

    Returns:
        a Namespace object with values of the command line arguments. 
    """
    snps_default = os.environ.get('XO_SNPS') or 'BSP_TIGER.marker_dataframe.pickle.gzip'
    intervals_default = os.environ.get('XO_INTERVALS') or 'BSP_TIGER.intervals_dataframe.pickle.gzip'
    crossovers_default = os.environ.get('XO_CROSSOVERS') or 'BSP_COs_final_set.pickle.gzip'
    peaks_default = os.environ.get('XO_PEAKS') or 'peaks.csv'
    filtered_default = os.environ.get('XO_BLOCKS') or 'filtered.csv'
    post_default = os.environ.get('XO_POST') or 'ncos.csv'

    # post_block_size = os.environ.get('XO_POST_BLOCK_SIZE') or (0,100)
    # post_block_length = os.environ.get('XO_POST_BLOCK_LENGTH') or (0,1000)
    # post_coverage = os.environ.get('XO_POST_COVERAGE') or 2
    # post_match = os.environ.get('XO_POST_MATCH') or False
    post_min_z = os.environ.get('XO_POST_MIN_Z') or 0.9
    post_delta_z = os.environ.get('XO_DELTA_HIGH_Z') or 0.1
    post_min_snps = os.environ.get('XO_POST_MIN_SNPS') or 2

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        title = 'subcommands',
        description = 'operation to perform',
        dest='cmnd'
    )

    parser.add_argument('--log', metavar='X', choices=['quiet','info','debug'])

    gui_parser = subparsers.add_parser('gui', help='explore blocks and NCOs')
    # gui_parser.add_argument('--intervals', metavar='F', default=intervals_default, help='SNP summaries')
    # gui_parser.add_argument('--peaks', metavar='F', default=peaks_default, help='blocks saved by peaks.py')
    # gui_parser.add_argument('--crossovers', metavar='F', default=crossovers_default, help='file with crossover locations')
    gui_parser.add_argument('--port', metavar='N', type=int, default=5006, help='local port for the Panel server')
    gui_parser.set_defaults(func=start_app)

    peak_parser = subparsers.add_parser('peaks', help='find blocks around peaks in the SNP data')
    peak_parser.add_argument('--snps', metavar='F', default=snps_default, help='input (TIGER marker) file')
    peak_parser.add_argument('--crossovers', metavar='F', default=crossovers_default, help='file with crossover locations')
    peak_parser.add_argument('--output', metavar='F', default=peaks_default, help='output file')
    peak_parser.add_argument('--max_snps', metavar='N', type=int, default=1000, help="max number of SNPs in a block")
    peak_parser.set_defaults(func=peak_finder)

    filter_parser = subparsers.add_parser('filter', help='apply filters to blocks')
    filter_parser.add_argument('--peaks', metavar='F', default=peaks_default, help='blocks saved by peaks.py')
    filter_parser.add_argument('--crossovers', metavar='F', default=crossovers_default, help='file with crossover locations')
    filter_parser.add_argument('--output', metavar='F', default=filtered_default, help='output file')
    filter_parser.add_argument('--chromosomes', metavar='P', default='BSP.*', help='chromosome name pattern')
    filter_parser.add_argument('--size', metavar='N', nargs=2, type=int, default=(0,100), help='block size range (#SNPs)')
    filter_parser.add_argument('--length', metavar='N', nargs=2, type=int, default=(0,10000), help='block length range (bp)')
    filter_parser.add_argument('--coverage', metavar='N', type=int, default=0, help='minimum coverage')
    filter_parser.add_argument('--match', action='store_true', help='require genome match')
    filter_parser.set_defaults(func=filter_blocks)

    post_parser = subparsers.add_parser('post', help='postprocessing of filtered blocks')
    post_parser.add_argument('--blocks', metavar='F', default=filtered_default, help='file with filtered blocks')
    post_parser.add_argument('--output', metavar='F', default=post_default, help='output file')
    post_parser.add_argument('--min_z', metavar='N', type=float, default=post_min_z, help='homozygosity for Type 2 blocks')
    post_parser.add_argument('--delta_z', metavar='N', type=float, default=post_delta_z, help='homozygosity for Type 1 blocks')
    post_parser.add_argument('--min_snps', metavar='N', type=int, default=post_min_snps, help='minimum number of SNPs of each type')
    post_parser.set_defaults(func=postprocess)

    if len(sys.argv) == 1:
        parser.print_help()
        exit()

    return parser.parse_args()

def setup_logging(args):
    """
    Configure the logging modile.  Uncomment one of the three levels
    for each command to define the logging level.
    """
    if args.cmnd == 'gui':
        level = logging.INFO
    else:
        match args.log:
            case 'info':
                level = logging.INFO
            case 'debug':
                level = logging.DEBUG
            case _:
                level = logging.WARNING
    logging.basicConfig(
        level=level,
        style='{',
        format='{relativeCreated:4.0f} msec: {message}',
        handlers = [RichHandler(markup=True, rich_tracebacks=True)],
    )

def main():
    """
    The argument parser associates a function with each script name, all
    we need to do here is call that function.
    """
    args = init_cli()
    setup_logging(args)
    try:
        args.func(args)
    except Exception as err:
        logging.exception(err)
