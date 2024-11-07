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

from .peaks import peak_finder
from .gui import start_app
from .vis import visualize, plot_commands
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
    save_default = os.environ.get('XO_SAVE') or 'summary.csv'

    post_block_size = os.environ.get('XO_POST_BLOCK_SIZE') or (0,100)
    post_block_length = os.environ.get('XO_POST_BLOCK_LENGTH') or (0,1000)
    post_coverage = os.environ.get('XO_POST_COVERAGE') or 2
    post_match = os.environ.get('XO_POST_MATCH') or True
    post_high_z = os.environ.get('XO_POST_HIGH_Z') or 0.9
    post_delta_z = os.environ.get('XO_DELTA_HIGH_Z') or 0.1

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        title = 'subcommands',
        description = 'operation to perform'
    )

    peak_parser = subparsers.add_parser('peaks', help='find peaks in the SNP data')
    peak_parser.add_argument('--snps', metavar='F', default=snps_default, help='input (TIGER marker) file')
    peak_parser.add_argument('--crossovers', metavar='F', default=crossovers_default, help='file with crossover locations')
    peak_parser.add_argument('--output', metavar='F', default=peaks_default, help='output file')
    peak_parser.add_argument('--max_snps', metavar='N', type=int, default=1000, help="max number of SNPs in a block")
    peak_parser.set_defaults(func=peak_finder)

    gui_parser = subparsers.add_parser('gui', help='explore blocks of SNPs')
    gui_parser.add_argument('--intervals', metavar='F', default=intervals_default, help='SNP summaries')
    gui_parser.add_argument('--peaks', metavar='F', default=peaks_default, help='blocks saved by peaks.py')
    gui_parser.add_argument('--crossovers', metavar='F', default=crossovers_default, help='file with crossover locations')
    gui_parser.add_argument('--port', metavar='N', type=int, default=5006, help='local port for the Panel server')
    gui_parser.set_defaults(func=start_app)

    vis_parser = subparsers.add_parser('vis', help='visualizations based on filtered blocks')
    vis_parser.add_argument('command', metavar='P', choices=plot_commands, help=f'type of plot to make {plot_commands}')
    vis_parser.add_argument('--peaks', metavar='F', default=peaks_default, help='blocks saved by peaks.py')
    vis_parser.add_argument('--chromosomes', metavar='P', default='BSP.*', help='chromosome name pattern')
    vis_parser.add_argument('--size', metavar='N', nargs=2, type=int, default=(0,100), help='block size range (#SNPs)')
    vis_parser.add_argument('--length', metavar='N', nargs=2, type=int, default=(0,10000), help='block length range (bp)')
    vis_parser.add_argument('--coverage', metavar='N', type=int, help='minimum coverage')
    vis_parser.add_argument('--match', action='store_true', help='require genome match')
    vis_parser.add_argument('--save', metavar='F', default=save_default, help='write summary dataframe to this file')
    vis_parser.set_defaults(func=visualize)

    post_parser = subparsers.add_parser('post', help='postprocessing of filtered blocks')
    post_parser.add_argument('--peaks', metavar='F', default=peaks_default, help='blocks saved by peaks.py')
    post_parser.add_argument('--chromosomes', metavar='P', default='BSP.*', help='chromosome name pattern')
    post_parser.add_argument('--size', metavar='N', nargs=2, type=int, default=post_block_size, help='block size range (#SNPs)')
    post_parser.add_argument('--length', metavar='N', nargs=2, type=int, default=post_block_length, help='block length range (bp)')
    post_parser.add_argument('--coverage', metavar='N', type=int, default=post_coverage, help='minimum coverage')
    post_parser.add_argument('--match', action='store_true', default=post_match, help='require genome match')
    post_parser.add_argument('--high_z', metavar='N', type=float, default=post_high_z, help='homozygosity for Type 2 blocks')
    post_parser.add_argument('--delta_z', metavar='N', type=float, default=post_delta_z, help='homozygosity for Type 1 blocks')
    post_parser.add_argument('--out', metavar='F', default=save_default, help='write processed data to this file')
    post_parser.set_defaults(func=postprocess)

    if len(sys.argv) == 1:
        parser.print_help()
        exit()

    return parser.parse_args()

def setup_logging(args):
    """
    Configure the logging modile.  Uncomment one of the first three
    lines to define the logging level.
    """
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

