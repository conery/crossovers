#
# Postprocessing of filtered blocks
#
# John Conery
# University of Oregon

import logging

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
    # filters; use the indexs in that frame to fetch the blocks of SNPs where
    # the data resides
 
    for blk_id, _ in summary.iterrows():
        # print(blk_id)
        # print(blocks.get_group(blk_id))
        block = blocks.get_group(blk_id)
        print(block)

    logging.info('postprocess')

    logging.info('exit')
