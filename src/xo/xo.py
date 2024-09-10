#! /usr/bin/env python3

# Top level application for the crossovers project.  Type
#
#  $ xo CMND [OPTS]
#
# to run one of the commands -- peak finder, viewer, or filter.

import argparse
import os
import sys

from .peaks import peak_finder
from .gui import start_app
from .vis import visualize

def init_cli():
    snps_default = os.environ.get('XO_SNPS') or 'BSP_TIGER.marker_dataframe.pickle.gzip'
    intervals_default = os.environ.get('XO_INTERVALS') or 'BSP_TIGER.intervals_dataframe.pickle.gzip'
    peaks_default = os.environ.get('XO_PEAKS') or 'peaks.csv'

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        title = 'subcommands',
        description = 'operation to perform'
    )

    peak_parser = subparsers.add_parser('peaks', help='find peaks in the SNP data')
    peak_parser.add_argument('--snps', metavar='F', default=snps_default, help='input (IGER marker) file')
    peak_parser.add_argument('--output', metavar='F', default=peaks_default, help='output file')
    peak_parser.add_argument('--max_snps', metavar='N', type=int, default=1000, help="max number of SNPs in a block")
    peak_parser.set_defaults(func=peak_finder)

    gui_parser = subparsers.add_parser('gui', help='explore blocks of SNPs')
    gui_parser.add_argument('--intervals', metavar='F', default=intervals_default, help='SNP summaries')
    gui_parser.add_argument('--peaks', metavar='F', default=peaks_default, help='blocks saved by peaks.py')
    gui_parser.add_argument('--port', metavar='N', type=int, default=5006, help='local port for the Panel server')
    gui_parser.set_defaults(func=start_app)

    vis_parser = subparsers.add_parser('vis', help='visualizations based on filtered blocks')
    vis_parser.add_argument('--peaks', metavar='F', default=peaks_default, help='blocks saved by peaks.py')
    vis_parser.add_argument('--chromosomes', metavar='P', default='BSP*', help='names of chromosomes to use')
    vis_parser.add_argument('--size', metavar='N', nargs=2, help='block size range (#SNPs)')
    vis_parser.add_argument('--length', metavar='N', nargs=2, help='block length range (bp)')
    vis_parser.add_argument('--coverage', metavar='N', help='minimum coverage')
    vis_parser.add_argument('--match', action='store_true', help='require genome match')
    vis_parser.set_defaults(func=visualize)

    post_parser = subparsers.add_parser('post', help='postprocessing of filtered blocks')
    post_parser.add_argument('--peaks', metavar='F', default=peaks_default, help='blocks saved by peaks.py')
    post_parser.set_defaults(func=post)

    if len(sys.argv) == 1:
        parser.print_help()
        exit()

    return parser.parse_args()

def post(args):
    print('TBD')

def main():
    args = init_cli()
    args.func(args)

