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
import os
import sys

import logging
from rich.logging import RichHandler

from .gui import start_app
from .peaks import peak_finder
# from .vis import visualize, plot_commands
from .post import postprocess

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
    match args.cmnd:
        case 'gui':
            level = logging.INFO
        case _:
            level = logging.INFO
            # level = logging.DEBUG
            # level = logging.WARNING
    logging.basicConfig(
        level=level,
        style='{',
        format='{relativeCreated:4.0f} msec: {message}',
        handlers = [RichHandler(markup=True)],
    )

def main():
    """
    The argument parser associates a function with each script name, all
    we need to do here is call that function.
    """
    args = init_cli()
    setup_logging(args)
    args.func(args)

