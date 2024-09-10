#
# Visualizations of filtered blocks
#
# John Conery
# University of Oregon

import pandas as pd
import panel as pn
import numpy as np
import re

import matplotlib.pyplot as plt

from .gui import SNPFilter

def load_data(args):
    df = pd.read_csv(args.peaks)
    res = []
    for chr_id, chr in df.groupby('chrom_id'):
        if not re.search(args.chromosomes, chr_id):
            continue
        for blk_id, blk in chr.groupby('blk_id'):
            res.append(blk)
    return res

def visualize(args):
    blocks = load_data(args)
    filter = SNPFilter()
    print(filter.widgets())
    print(filter.widget_map())
