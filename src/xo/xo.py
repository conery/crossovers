#! /usr/bin/env python3

# Top level application for the crossovers project.  Type
#
#  $ xo CMND [OPTS]
#
# to run one of the commands -- peak finder, viewer, or filter.

import argparse
from peaks import main
from viewer import start_app

def init_cli():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        title = 'subcommands',
        description = 'operation to perform'
    )

    peak_parser = subparsers.add_parser('peaks', help='find peaks in the SNP data')
    peak_parser.add_argument('--snps', metavar='F', default='BSP_TIGER.marker_dataframe.pickle.gzip', help='input (IGER marker) file')
    peak_parser.add_argument('--output', metavar='F', default='peaks.pickle.gzip', help='output file')
    peak_parser.set_defaults(func=main)

    view_parser = subparsers.add_parser('view', help='explore blocks of SNPs')
    view_parser.add_argument('--intervals', metavar='F', default='BSP_TIGER.intervals_dataframe.pickle.gzip', help='SNP summaries')
    view_parser.add_argument('--peaks', metavar='F', default='peaks.pickle.gzip', help='blocks saved by peaks.py')
    view_parser.add_argument('--port', metavar='N', type=int, default=5006, help='local port for the Panel server')
    view_parser.set_defaults(func=start_app)

    return parser.parse_args()

if __name__ == '__main__':
    args = init_cli()
    args.func(args)
