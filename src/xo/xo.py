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
from .viewer import start_app

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

    view_parser = subparsers.add_parser('view', help='explore blocks of SNPs')
    view_parser.add_argument('--intervals', metavar='F', default=intervals_default, help='SNP summaries')
    view_parser.add_argument('--peaks', metavar='F', default=peaks_default, help='blocks saved by peaks.py')
    view_parser.add_argument('--port', metavar='N', type=int, default=5006, help='local port for the Panel server')
    view_parser.set_defaults(func=start_app)

    if len(sys.argv) == 1:
        parser.print_help()
        exit()

    return parser.parse_args()

def main():
    args = init_cli()
    args.func(args)

