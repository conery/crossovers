
# Panel application for viewing peaks in SNP data
#
# John Conery
# University of Oregon

import io
import pandas as pd
import panel as pn
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle

from xo.filters import SNPFilter

pn.extension('tabulator')

SIDEBAR_WIDTH = 350

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
        self.filter = None

    def filter_cb(self):
        self.filter.set_size(self.value)


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
        self.filter = None

    def filter_cb(self):
        self.filter.set_length(self.value)


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
        self.filter = None

    def filter_cb(self):
        self.filter.set_coverage(self.value)


class SupportFilter(pn.widgets.Checkbox):
    '''
    Keep SNPs if the base genome column matches the HMM state column.
    '''
    def __init__(self):
        super(SupportFilter,self).__init__(
            name = 'Genome Match',
        )
        self.tags = ['match', 0]
        self.filter = None

    def filter_cb(self):
        self.filter.set_matched(self.value)


class FilterBox(pn.Column):

    def __init__(self):
        '''
        Instantiate the widgets that will filter SNPs and put them in a Column
        in the order they will be displayed.  Then sort them by priority (the
        second value in each widget's tag list) so the list returned by calling
        widgets() is the order they are applied.  The widget map associates a
        filter name (the first item in the tag list) with the filter widget.
        '''
        super(FilterBox, self).__init__()
        self._widgets = [cls() for cls in [BlockSizeFilter, BlockLengthFilter, CoverageFilter, SupportFilter]]
        self.append(pn.pane.HTML("<h3>Filters</h3>"))
        self.extend(self._widgets)
        self._widgets.sort(key = lambda w: w.tags[1])
        # self._widget_map = { w.tags[0]: w for w in self._widgets }

    def widgets(self):
        return self._widgets
    
    # def widget_map(self):
    #     return self._widget_map
    
    def set_filter(self, f):
        for w in self._widgets:
            w.filter = f


class PeakViewerApp(pn.template.BootstrapTemplate):
    def __init__(self, **params):
        """
        Initialize the application.

        Arguments:
          params:  runtime options passed to the parent class constructor
        """
        super(PeakViewerApp, self).__init__(**params)

        self.filter = SNPFilter()
        self.filter_widgets = FilterBox()
        for w in self.filter_widgets.widgets():
            w.param.watch(self.filter_cb, ['value'])

        button_style_sheet = ''':host(.solid) .bk-btn {
            --color: white;
        }
        '''

        self.back_button = pn.widgets.Button(name='◀︎', stylesheets=[button_style_sheet])
        self.forward_button = pn.widgets.Button(name='▶︎', stylesheets=[button_style_sheet])
        self.chromosome_id = pn.widgets.TextInput(name="", value="")

        for w in [self.back_button, self.forward_button]:
            w.param.watch(self.change_chromosome_cb, ['value'])

        self.chromosome_pattern = pn.widgets.TextInput(name="Chromosomes", value="BSP.*")
        self.size_graph_button = pn.widgets.Button(name='Block Size', stylesheets=[button_style_sheet])
        self.length_graph_button = pn.widgets.Button(name='Block Length', stylesheets=[button_style_sheet])
        self.location_graph_button = pn.widgets.Button(name='Block Location', stylesheets=[button_style_sheet])

        self.download_button = pn.widgets.FileDownload(callback=self.download_cb, filename='summary.csv', align='center', icon='download')
        self.download_pane = pn.GridBox(self.download_button, height=200, width=SIDEBAR_WIDTH)
        self.download_button.visible = False

        for w in [self.size_graph_button, self.length_graph_button, self.location_graph_button]:
            w.param.watch(self.summary_plot_cb, ['value'])

        self.chromosome_id.param.watch(self.chromosome_edited_cb, ['value'])      

        chr_tab = pn.Column(
            pn.pane.HTML('<h3>Chromosome</h3>'),
            pn.Row(self.back_button, self.chromosome_id, self.forward_button),
            pn.pane.HTML('<p>Placeholder</p>'),
        )

        summ_tab = pn.Column(
            pn.pane.HTML('<h3>Summary</h3>'),
            self.chromosome_pattern,
            pn.Row(self.size_graph_button, self.length_graph_button, self.location_graph_button),
            pn.pane.HTML('<p>Click a button above to generate a plot summarizing all chromosomes.</p>')
        )

        self.tabs = pn.Tabs(
            ('Chromosome', chr_tab),
            ('Summary', summ_tab),
        )

        self.sidebar.append(
            pn.Column(
                self.filter_widgets,
                self.download_pane,
            )
        )
        self.main.append(self.tabs)

    def load_data(self, args):
        print('loading interval data')
        self.intervals = pd.read_pickle(args.intervals, compression='gzip').groupby('chrom_id')
        self.clist = list(self.intervals.groups.keys())
        self.cmap = { name: i for i, name in enumerate(self.clist)}
        print('loading peak data')
        self.filter.load_data(args.peaks)
        self.filter_widgets.set_filter(self.filter)

        # setting a value in the chromosome name widget triggers an update
        # to the graphic to display the first chromosome
        self.chr_index = 0
        self.chromosome_id.value = self.clist[0]

    def display_chromosome(self):
        chr_id = self.chromosome_id.value
        chrom = self.intervals.get_group(chr_id)
        rects = PatchCollection(self._make_patches(chrom), match_original=True)
        fig, ax = plt.subplots(figsize=(12,1))
        plt.box(False)
        plt.yticks([])
        plt.xticks(ticks=np.linspace(0,20000000,5), labels=[f'{int(n*20)}Mbp' for n in np.linspace(0,1,5)])
        ax.xaxis.set_ticks_position('top')
        ax.add_collection(rects)
        plt.xlim(0,20000000)
        plt.ylim(0,2000000)
        plt.close(fig)
        graphic = pn.Column(pn.pane.Matplotlib(fig, dpi=72, tight=True))
        if self.filter.has_chromosome_block(chr_id):
            self.blocks, self.summary = self.filter.apply(chr_id)
            # self._make_dots()
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
        return res

    # def _make_dots(self):
    def _make_grid(self):
        pcolor = {
            'CB4856': 'cornflowerblue',
            'N2': 'indianred',
            'uCB4856': 'lightsteelblue',
            'uN2': 'lightpink',
            'unknown': 'lightgray',
            'het': 'palegoldenrod',
        }
        self.block_buttons = {}
        self.block_text = {}
        g = pn.Column()
        for blk_id, blk_stats in self.summary.iterrows():
            block = self.blocks.get_group(blk_id)
            fig, ax = plt.subplots(figsize=(10,0.8))
            plt.box(False)
            plt.xlim(0,10)
            plt.ylim(0,0.8)
            plt.yticks([])
            x0 = block.iloc[0].position
            size = int(blk_stats.blk_size)
            length = int(blk_stats.blk_len)
            w = block.iloc[-1].position - x0
            plt.xticks(ticks=np.linspace(0,10,5), labels=[f'{int(n*w)}bp' for n in np.linspace(0,1,5)])
            plt.suptitle(f'Block #{blk_id}\nStart: {(x0/1000000):.1f}Mbp\nSize: {size} SNPs\nLength: {length}bp', x=0, y=0.75, size='medium',ha='left')
            res = []
            for _, snp in block.iterrows():
                p = ((snp.position - x0) / length) if length > 0 else 0
                x = p*10
                res.append(Circle((x,0.2),0.1,color=pcolor[snp.base_geno]))
            dots = PatchCollection(res, match_original=True)
            ax.add_collection(dots)
            plt.close(fig)
            self.block_buttons[blk_id] = pn.widgets.Button(name='>', align='center', tags=[blk_id])
            self.block_buttons[blk_id].on_click(self.toggle_text)
            df = block[['position','base_geno','hmm_state1','reference','ref_reads','variant','var_reads']]
            self.block_text[blk_id] = pn.pane.DataFrame(df, visible=False)
            g.append(pn.Row(
                self.block_buttons[blk_id],
                pn.pane.Matplotlib(fig, dpi=72, tight=True),
                styles={'background':'WhiteSmoke'},
            ))
            g.append(pn.Row(self.block_text[blk_id]))
        return g
    
    def toggle_text(self, e):
        i = e.obj.tags[0]
        self.block_text[i].visible = not self.block_text[i].visible
        self.block_buttons[i].name = '∨' if self.block_text[i].visible else '>'

    def filter_cb(self, e):
        e.obj.filter_cb()
        self.display_chromosome()        

    def change_chromosome_cb(self, e):
        delta = 1 if e.obj is self.forward_button else -1
        self.chr_index = (self.chr_index + delta) % len(self.clist)
        self.chromosome_id.value = self.clist[self.chr_index]

    def chromosome_edited_cb(self, e):
        idx = self.cmap.get(e.obj.value)
        if idx is not None:
            self.chr_index = idx
            self.chromosome_id.value = self.clist[idx]
            self.display_chromosome()

    histogram_params = {
        'Block Size': {
            'col':     'blk_size',
            'title':   'Block Size',
            'xlabel':  'Number of SNPs',
            'ylabel':  'Number of Blocks',
            'hist': {
                'bins': 10,
                'rwidth': 0.8,
                'align': 'left',
                'range': (1,100),
            },
        },
        'Block Length': {
            'col':     'blk_len',
            'title':   'Block Length',
            'xlabel':  'Length (bp)',
            'ylabel':  'Number of Blocks',
            'hist': {
                'bins': 10,
                'rwidth': 0.8,
            },
        },
        'Block Location': {
            'col':     'blk_loc',
            'title':   'Block Location',
            'xlabel':  'Relative Position in the Chromosome',
            'ylabel':  'Number of Blocks',
            'hist': {
                'bins': 100,
                'range': (0,1),
            },
        },
    }

    def summary_plot_cb(self, e):
        params = self.histogram_params[e.obj.name]
        self.tabs[1].loading = True
        self.filter.set_chromosome(self.chromosome_pattern.value)
        self.summary_df = self.filter.summary()
        fig, ax = plt.subplots(figsize=(7,5))
        plt.hist(self.summary_df[params['col']], label=self.chromosome_pattern.value, **params['hist'])
        plt.title(params['title'])
        plt.xlabel(params['xlabel'])
        plt.ylabel(params['ylabel'])
        plt.legend(handlelength=0)
        plt.close(fig)
        self.tabs[1].loading = False
        self.tabs[1].pop(-1)
        self.tabs[1].append(pn.pane.Matplotlib(fig, dpi=72, tight=True))
        self.download_button.visible = True

    def download_cb(self):
        sio = io.StringIO()
        self.summary_df.to_csv(sio)
        sio.seek(0)
        return sio


def make_app(args):
    """
    Instantiate the top level widget.

    Returns:
        a PeakViewerApp object
    """
    app = PeakViewerApp(
        title='NCO Explorer', 
        sidebar_width=SIDEBAR_WIDTH,
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


