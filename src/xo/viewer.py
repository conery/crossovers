
# Panel application for viewing peaks in SNP data
#
# John Conery
# University of Oregon

import pandas as pd
import panel as pn
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle

pn.extension('tabulator')

class BlockSizeFilter(pn.widgets.IntRangeSlider):
    def __init__(self):
        super(BlockSizeFilter,self).__init__(
            name = 'Block Size (#SNPs)',
            start = 0,
            end = 100,
        )
        self.value = (self.start, self.end)

    def filter(self, df):
        bsmin, bsmax = self.value
        res = set() if ((len(df) < bsmin) or (len(df) > bsmax)) else set(df.index)
        return res

class BlockLengthFilter(pn.widgets.IntRangeSlider):
    def __init__(self):
        super(BlockLengthFilter,self).__init__(
            name = 'Block Length (bp)',
            start = 0,
            end = 10000,
            step = 10,
        )
        self.value = (self.start, self.end)            

    def filter(self, df):
        w = df.iloc[-1].position - df.iloc[0].position
        blmin, blmax = self.value
        res = set() if ((w < blmin) or (w > blmax)) else set(df.index)
        return res

class CoverageFilter(pn.widgets.IntSlider):
    def __init__(self):
        super(CoverageFilter,self).__init__(
            name = 'Minimum Coverage',
            start = 0,
            end = 10,
        )
        self.value = self.start

    def filter(self, df):
        return set(df[(df.ref_reads + df.var_reads) >= self.value].index)

class SupportFilter(pn.widgets.Checkbox):
    def __init__(self):
        super(SupportFilter,self).__init__(
            name = 'Genome Match',
        )

    def filter(self, df):
        res = set(df[df.base_geno == df.hmm_state1].index) if self.value else set(df.index)
        return res

class SNPFilter(pn.Column):

    def __init__(self):
        super(SNPFilter, self).__init__()
        self.widgets = [cls() for cls in [BlockSizeFilter, BlockLengthFilter, CoverageFilter, SupportFilter]]
        self.append(pn.pane.HTML("<h3>Filters</h3>"))
        self.extend(self.widgets)

    def widgets(self):
        return self.widgets

    def apply(self, df):
        res = set(df.index)
        i = 0
        while i < len(self.widgets) and len(res) > 0:
            res &= self.widgets[i].filter(df)
            i += 1
        return res

class PeakViewerApp(pn.template.BootstrapTemplate):
    def __init__(self, **params):
        """
        Initialize the application.

        Arguments:
          params:  runtime options passed to the parent class constructor
        """
        super(PeakViewerApp, self).__init__(**params)

        self.filter = SNPFilter()
        for w in self.filter.widgets:
            w.param.watch(self.filter_cb, ['value'])

        button_style_sheet = ''':host(.solid) .bk-btn {
            --color: white;
        }
        '''

        self.back_button = pn.widgets.Button(name='◀︎', stylesheets=[button_style_sheet])
        self.forward_button = pn.widgets.Button(name='▶︎',styles={'background':'white'})
        for w in [self.back_button, self.forward_button]:
            w.param.watch(self.change_chromosome_cb, ['value'])
        
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
            ('Summary', summ_tab),
        )

        self.sidebar.append(self.filter)
        self.main.append(self.tabs)

        self.chr_index = 0

    def load_data(self, args):
        print('loading interval data')
        self.intervals = pd.read_pickle(args.intervals, compression='gzip')
        self.clist = sorted(self.intervals['chrom_id'].unique())
        print('loading peak data')
        self.peaks = pd.read_pickle(args.peaks)
        self.display_chromosome()

    def display_chromosome(self):
        chr_id = self.clist[self.chr_index]
        chrom = self.intervals[self.intervals.chrom_id == chr_id]
        peaks = self.peaks[self.peaks.chrom_id == chr_id]
        rects = PatchCollection(self._make_patches(chrom), match_original=True)
        self._make_dots(peaks)
        fig, ax = plt.subplots(figsize=(12,1))
        plt.box(False)
        plt.yticks([])
        plt.xticks(ticks=np.linspace(0,20000000,5), labels=[f'{int(n*20)}Mbp' for n in np.linspace(0,1,5)])
        ax.xaxis.set_ticks_position('top')
        ax.add_collection(rects)
        plt.xlim(0,20000000)
        plt.ylim(0,2000000)
        plt.close(fig)
        self.chromosome_id.value = chr_id
        graphic = pn.Column(pn.pane.Matplotlib(fig, dpi=72, tight=True))
        grid = self._make_grid()
        graphic.append(grid)
        self.tabs[0].pop(-1)
        self.tabs[0].append(graphic)

    def _make_patches(self, df):
        pcolor = {
            'CB4856': 'dodgerblue',
            'N2': 'indianred'
        }
        res = []
        for _, r in df.iterrows():
            c = pcolor.get(r.hmm_state) or 'lightgray'
            res.append(Rectangle((r.start,500000), r.length, 1000000, color=c))
            res.append(Circle((r.start,750000), 50000, color='black'))
        print(f'{len(res)} patches')
        return res

    def _make_dots(self, df):
        pcolor = {
            'CB4856': 'cornflowerblue',
            'N2': 'indianred',
            'uCB4856': 'lightsteelblue',
            'uN2': 'lightpink',
            'unknown': 'lightgray',
            'het': 'palegoldenrod',
        }
        self.dotfigs = {}
        self.blocks = {}
        for grp_id, grp in df.groupby('blk_id'):
            if snps := self.filter.apply(grp):
                blk = grp[grp.index.isin(snps)]
                fig, ax = plt.subplots(figsize=(10,0.8))
                plt.box(False)
                plt.xlim(0,10)
                plt.ylim(0,0.8)
                plt.yticks([])
                x0 = blk.iloc[0].position
                w = blk.iloc[-1].position - x0
                plt.xticks(ticks=np.linspace(0,10,5), labels=[f'{int(n*w)}bp' for n in np.linspace(0,1,5)])
                plt.suptitle(f'Block #{grp_id}\nStart: {(x0/1000000):.1f}Mbp\nSize: {len(grp)} SNPs\nLength: {w}bp', x=0, y=0.75, size='medium',ha='left')
                res = []
                for _, snp in blk.iterrows():
                    p = ((snp.position - x0) / w) if w > 0 else 0
                    x = p*10
                    res.append(Circle((x,0.2),0.1,color=pcolor[snp.base_geno]))
                dots = PatchCollection(res, match_original=True)
                ax.add_collection(dots)
                plt.close(fig)
                self.dotfigs[grp_id] = (x0, w, fig)
                self.blocks[grp_id] = blk[['position','base_geno','hmm_state1','reference','ref_reads','variant','var_reads']]

    def _make_grid(self):
        if len(self.dotfigs) == 0:
            return None
        styles = {'border': '1px solid black'}
        self.block_buttons = {}
        self.block_text = {}
        g = pn.Column()
        for i in self.dotfigs.keys():
            startloc, width, fig = self.dotfigs[i]
            self.block_buttons[i] = pn.widgets.Button(name='>', align='center', tags=[i])
            self.block_buttons[i].on_click(self.toggle_text)
            self.block_text[i] = pn.pane.DataFrame(self.blocks[i], visible=False)
            g.append(pn.Row(
                self.block_buttons[i],
                pn.pane.Matplotlib(fig, dpi=72, tight=True),
                styles={'background':'WhiteSmoke'},
            ))
            g.append(pn.Row(self.block_text[i]))
        return g
    
    def toggle_text(self, e):
        i = e.obj.tags[0]
        self.block_text[i].visible = not self.block_text[i].visible
        self.block_buttons[i].name = '∨' if self.block_text[i].visible else '>'

    def filter_cb(self, e):
        self.display_chromosome()        

    def change_chromosome_cb(self, e):
        delta = 1 if e.obj is self.forward_button else -1
        self.chr_index = (self.chr_index + delta) % len(self.clist)
        self.display_chromosome()        

def make_app(args):
    """
    Instantiate the top level widget.

    Returns:
        a PeakViewerApp object
    """
    app = PeakViewerApp(
        title='NCO Explorer', 
        sidebar_width=350,
    )
    app.load_data(args)
    return app

def start_app(args):
    """
    Launch the Bokeh server.
    """
    pn.extension(design='native')
    pn.config.throttled = True
    try:
        app = make_app(args)
        pn.serve( 
            {'peaks': app},
            port = args.port,
            verbose = True,
            autoreload = True,
            websocket_origin= '*',
        )
    except Exception as err:
        print(err)


