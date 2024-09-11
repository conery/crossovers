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

def set_params(filters, args):
    params = vars(args)
    for f in filters:
        if p := params.get(f):
            if isinstance(p, list):
                filters[f].value = tuple(p)
            else:
                filters[f].value = p

def visualize(args):
    filter = SNPFilter()
    set_params(filter.widget_map(), args)
    blocks = [filter.apply(b) for b in load_data(args)]
    for b in blocks:
        print(b)

