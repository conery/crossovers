#
# Global settings shared across modules
#

import os

class Config:
    # names of files with SNP data from the BSP project
    snps_default = os.environ.get('XO_SNPS') or 'BSP_TIGER.marker_dataframe.pickle.gzip'
    intervals_default = os.environ.get('XO_INTERVALS') or 'BSP_TIGER.intervals_dataframe.pickle.gzip'
    crossovers_default = os.environ.get('XO_CROSSOVERS') or 'BSP_COs_final_set.pickle.gzip'

    # peak finder defaults
    peaks_output_default = os.environ.get('XO_PEAKS') or 'peaks.csv'
    peaks_max_snps = 1000

    # filter defaults
    filter_output_default = os.environ.get('XO_BLOCKS') or 'filtered.csv'
    filter_chromosome_pattern = 'BSP.*'
    filter_block_size = (0,100)
    filter_block_length = (0,10000)
    filter_coverage = 0

    # postprocessor defaults
    post_output_default = os.environ.get('XO_POST') or 'ncos.csv'
    post_min_z = 0.9
    post_delta_z = 0.1
    post_min_cover = 2
    post_block_size = 5

    # GUI settings
    sidebar_width = 350
    terminal_height = 200

# Chromosome lengths

chr_length = {
    1: 15114068,
    2: 15311845,
    3: 13819453,
    4: 17493838,
    5: 20953657,
    6: 17739129,
}
