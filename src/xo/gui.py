
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
    '''
    A filter based on the count of the number of SNPs in a block.
    '''
    def __init__(self):
        super(BlockSizeFilter,self).__init__(
            name = 'Block Size (#SNPs)',
            start = 0,
            end = 100,
        )
        self.value = (self.start, self.end)
        self.tags = ['size', 1]

    def filter(self, df):
        '''
        If the block size (number of rows in the frame) is outside the range 
        return an empty block, otherwise return the block itself.
        '''
        bsmin, bsmax = self.value
        if ((len(df) < bsmin) or (len(df) > bsmax)):
            return pd.DataFrame(columns=df.columns)
        else:
            return df

class BlockLengthFilter(pn.widgets.IntRangeSlider):
    '''
    A filter based on the length of a block, defined as the number of bases
    between the first and last SNP.
    '''
    def __init__(self):
        super(BlockLengthFilter,self).__init__(
            name = 'Block Length (bp)',
            start = 0,
            end = 10000,
            step = 10,
        )
        self.value = (self.start, self.end)
        self.tags = ['length', 1]       

    def filter(self, df):
        '''
        If the block length is outside the range return an empty block, 
        otherwise return the block itself.
        '''
        w = df.iloc[-1].position - df.iloc[0].position
        blmin, blmax = self.value
        if ((w < blmin) or (w > blmax)):
            return pd.DataFrame(columns=df.columns)
        else:
            return df

class CoverageFilter(pn.widgets.IntSlider):
    '''
    Keep SNPs if the number of reads (in either the reference or variant column)
    is greater than a cutoff value. 
    '''
    def __init__(self):
        super(CoverageFilter,self).__init__(
            name = 'Minimum Coverage',
            start = 0,
            end = 10,
        )
        self.value = self.start
        self.tags = ['coverage', 0]

    def filter(self, df):
        return df[(df.ref_reads + df.var_reads) >= self.value]

class SupportFilter(pn.widgets.Checkbox):
    '''
    Keep SNPs if the base genome column matches the HMM state column.
    '''
    def __init__(self):
        super(SupportFilter,self).__init__(
            name = 'Genome Match',
        )
        self.tags = ['match', 0]

    def filter(self, df):
        if self.value:
            return df[df.base_geno == df.hmm_state1]
        else:
            return df

class SNPFilter(pn.Column):

    def __init__(self):
        '''
        Instantiate the widgets that will filter SNPs and put them in a Column
        in the order they will be displayed.  Then sort them by priority (the
        second value in each widget's tag list) so the list returned by calling
        widgets() is the order they are applied.  The widget map associates a
        filter name (the first item in the tag list) with the filter widget.
        '''
        super(SNPFilter, self).__init__()
        self._widgets = [cls() for cls in [BlockSizeFilter, BlockLengthFilter, CoverageFilter, SupportFilter]]
        self.append(pn.pane.HTML("<h3>Filters</h3>"))
        self.extend(self._widgets)
        self._widgets.sort(key = lambda w: w.tags[1])
        self._widget_map = { w.tags[0]: w for w in self._widgets }

    def widgets(self):
        return self._widgets
    
    def widget_map(self):
        return self._widget_map

    def apply(self, df):
        res = set(df.index)
        i = 0
        while i < len(self._widgets) and len(df) > 0:
            df = self._widgets[i].filter(df)
            i += 1
        return df

class PeakViewerApp(pn.template.BootstrapTemplate):
    def __init__(self, **params):
        """
        Initialize the application.

        Arguments:
          params:  runtime options passed to the parent class constructor
        """
        super(PeakViewerApp, self).__init__(**params)

        self.filter = SNPFilter()
        for w in self.filter.widgets():
            w.param.watch(self.filter_cb, ['value'])

        button_style_sheet = ''':host(.solid) .bk-btn {
            --color: white;
        }
        '''

        self.back_button = pn.widgets.Button(name='◀︎', stylesheets=[button_style_sheet])
        self.forward_button = pn.widgets.Button(name='▶︎',styles={'background':'white'})
        # self.chromosome_id = pn.widgets.StaticText(name="", value="")
        self.chromosome_id = pn.widgets.TextInput(name="", value="")

        for w in [self.back_button, self.forward_button]:
            w.param.watch(self.change_chromosome_cb, ['value'])

        self.chromosome_id.param.watch(self.chromosome_edited_cb, ['value'])      

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
        self.cmap = { name: i for i, name in enumerate(self.clist)}
        print('loading peak data')
        self.peaks = pd.read_csv(args.peaks)
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
            snps = self.filter.apply(grp)
            if len(snps) > 0:
                # blk = grp[grp.index.isin(snps)]
                fig, ax = plt.subplots(figsize=(10,0.8))
                plt.box(False)
                plt.xlim(0,10)
                plt.ylim(0,0.8)
                plt.yticks([])
                x0 = snps.iloc[0].position
                w = snps.iloc[-1].position - x0
                plt.xticks(ticks=np.linspace(0,10,5), labels=[f'{int(n*w)}bp' for n in np.linspace(0,1,5)])
                plt.suptitle(f'Block #{grp_id}\nStart: {(x0/1000000):.1f}Mbp\nSize: {len(snps)} SNPs\nLength: {w}bp', x=0, y=0.75, size='medium',ha='left')
                res = []
                for _, snp in snps.iterrows():
                    p = ((snp.position - x0) / w) if w > 0 else 0
                    x = p*10
                    res.append(Circle((x,0.2),0.1,color=pcolor[snp.base_geno]))
                dots = PatchCollection(res, match_original=True)
                ax.add_collection(dots)
                plt.close(fig)
                self.dotfigs[grp_id] = (x0, w, fig)
                self.blocks[grp_id] = snps[['position','base_geno','hmm_state1','reference','ref_reads','variant','var_reads']]

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

    def chromosome_edited_cb(self, e):
        if idx := self.cmap.get(e.obj.value):
            self.chr_index = idx
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


