#! /usr/bin/env python3

#
# Panel application for viewing peaks in SNP data
#
# John Conery
# University of Oregon
# (conery@uoregon.edu)
#

import argparse
import pandas as pd
import panel as pn

import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle

def init_cli():
    """
    Use argparse to create the command line API.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--peaks', metavar='F', default='peaks.pickle.gzip', help='blocks saved by peaks.py')
    parser.add_argument('--port', metavar='N', type=int, default=5006, help='local port for the Panel server')

    return parser.parse_args()

class PeakViewerApp(pn.template.BootstrapTemplate):
    def __init__(self, **params):
        """
        Initialize the application.

        Arguments:
          params:  runtime options passed to the parent class constructor
        """
        super(PeakViewerApp, self).__init__(**params)

        self.back_button = pn.widgets.Button(name='◀︎')
        self.back_button.on_click(self.prev_chromosome)

        self.forward_button = pn.widgets.Button(name='▶︎')
        self.forward_button.on_click(self.next_chromosome)
        
        self.chromosome_id = pn.widgets.StaticText(name="", value="")

        chr_tab = pn.Column(
            pn.pane.HTML('<h3>Chromosome</h3>'),
            pn.Row(self.back_button, self.chromosome_id, self.forward_button),
            pn.pane.HTML('<p>Placeholder</p>'),
        )

        summ_tab = pn.Column(
            pn.pane.HTML('<h3>Summary</h3>'),
            pn.pane.HTML('<p>TBD</p>')
        )

        self.tabs = pn.Tabs(
            ('Chromosome', chr_tab),
            ('Summaary', summ_tab),
        )

        self.sidebar.append(pn.pane.HTML("<h3>Sidebar</h3>"))
        self.main.append(self.tabs)

        self.chr_index = 0

    def load_data(self, args):
        print('loading interval data')
        self.intervals = pd.read_pickle('BSP_TIGER.intervals_dataframe.pickle.gzip', compression='gzip')
        self.clist = sorted(self.intervals['chrom_id'].unique())
        print('loading peak data')
        self.peaks = pd.read_pickle(args.peaks)
        self.display_chromosome()

    def display_chromosome(self):
        chr_id = self.clist[self.chr_index]
        chrom = self.intervals[self.intervals.chrom_id == chr_id]
        rects = PatchCollection(self._make_patches(chrom), match_original=True)
        pks = self.peaks[self.peaks.chrom_id == chr_id]
        nrows, p = self._make_circles(pks)
        dots = PatchCollection(p, match_original=True)
        fig, ax = plt.subplots(figsize=(10,(nrows+1)/2))
        plt.box(False)
        plt.yticks([])
        ax.xaxis.set_ticks_position('top')
        ax.add_collection(rects)
        ax.add_collection(dots)
        plt.xlim(0,20000000)
        plt.ylim(-1000000*(nrows+1),2000000)
        plt.close(fig)
        self.chromosome_id.value = chr_id
        self.tabs[0].pop(-1)
        self.tabs[0].append(pn.pane.Matplotlib(fig, dpi=72, tight=True))

    def _make_patches(self, df):
        pcolor = {
            'CB4856': 'dodgerblue',
            'N2': 'indianred'
        }
        res = []
        for _, r in df.iterrows():
            c = pcolor.get(r.hmm_state) or 'lightgray'
            res.append(Rectangle((r.start,0), r.length, 500000, color=c))
        print(f'{len(res)} patches')
        return res

    def _make_circles(self, df):
        pcolor = {
            'CB4856': 'cornflowerblue',
            'N2': 'indianred',
            'uCB4856': 'lightsteelblue',
            'uN2': 'lightpink',
            'unknown': 'lightgray',
            'het': 'palegoldenrod',
        }
        res = []
        rownum = 0
        for grp_id, grp in df.groupby('blk_id'):
            w = grp.iloc[-1].position - grp.iloc[0].position
            m =  grp.iloc[0].position + w/2
            x0 = m - 1250000
            rownum += 1
            for _, snp in grp.iterrows():
                p = (snp.position - grp.iloc[0].position) / w
                x = x0 + p*2500000
                res.append(Circle((x,-1000000*rownum),100000,color=pcolor[snp.base_geno]))
        print(f'{rownum} rows')
        return rownum, res
    
    def next_chromosome(self,_):
        self.chr_index = (self.chr_index + 1) % len(self.clist)
        self.display_chromosome()

    def prev_chromosome(self,_):
        self.chr_index = (self.chr_index - 1) % len(self.clist)
        self.display_chromosome()

def make_app(args):
    """
    Instantiate the top level widget.

    Returns:
        a PeakViewerApp object
    """
    app = PeakViewerApp(
        title='NCO Viewer', 
        sidebar_width=450,
    )
    app.load_data(args)
    return app

def start_app(args):
    """
    Launch the Bokeh server.
    """
    pn.extension(design='native')
    app = make_app(args)
    pn.serve( 
        {'peaks': app},
        port = args.port,
        verbose = True,
        autoreload = True,
        websocket_origin= '*',
    )

if __name__ == '__main__':
    args = init_cli()
    start_app(args)

