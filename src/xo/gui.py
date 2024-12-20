
# Panel application for viewing peaks in SNP data
#
# John Conery
# University of Oregon
#
# Run the application using the top level xo command:
#
#  $ xo gui OPTS
#

import io
import logging
import numpy as np
import pandas as pd
import panel as pn

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle

from xo.filters import SNPFilter, NCOFilter
from xo.config import Config, chr_length

pn.extension('tabulator')

SIDEBAR_WIDTH = 350

class BlockSizeFilterWidget(pn.widgets.IntRangeSlider):
    """
    Use an integer range slider to provide settings for the
    filter based on the count of the number of SNPs in a block.
    """

    def __init__(self, f):
        '''
        Arguments:
          f:  the SNPFilter object that will do the filtering.
        '''
        c = Config()
        super(BlockSizeFilterWidget,self).__init__(
            name = 'Block Size (#SNPs)',
            start = c.filter_block_size[0],
            end = c.filter_block_size[1]
        )
        self.value = (self.start, self.end)
        self.tags = ['size', 1]
        self.filter = f

    def filter_cb(self):
        '''
        Callback activated when the slider changes.  Saves the new
        slider settings in the filter object.
        '''
        self.filter.size_range = self.value


class BlockLengthFilterWidget(pn.widgets.IntRangeSlider):
    """
    Use an integer range slider to provide settings for the
    filter based on length of a block, defined as the number
    of bases between the first and last SNP.
    """

    def __init__(self, f):
        '''
        Arguments:
          f:  the SNPFilter object that will do the filtering.
        '''
        c = Config()
        super(BlockLengthFilterWidget,self).__init__(
            name = 'Block Length (bp)',
            start = c.filter_block_length[0],
            end = c.filter_block_length[1],
            step = 100,
        )
        self.value = (self.start, self.end)
        self.tags = ['length', 1]     
        self.filter = f

    def filter_cb(self):
        '''
        Callback activated when the slider changes.  Saves the new
        slider settings in the filter object.
        '''
        self.filter.length_range = self.value


class CoverageFilterWidget(pn.widgets.IntSlider):
    """
    Use an integer slider to display the coverage cutoff value. 
    """

    def __init__(self, f):
        '''
        Arguments:
          f:  the filter object with the method that does the filtering
        '''
        super(CoverageFilterWidget,self).__init__(
            name = 'Minimum Coverage',
            start = 0,
            end = 10,
        )
        self.value = self.start
        self.tags = ['coverage', 0]
        self.filter = f

    def filter_cb(self):
        '''
        Callback activated when the slider changes.  Saves the new
        slider settings in the filter object.
        '''
        self.filter.coverage = self.value


class SupportFilterWidget(pn.widgets.Checkbox):
    """
    Use a checkbox to tell the filter to keep SNPs if the base genome
    column matches the HMM state column.
    """

    def __init__(self, f):
        '''
        Arguments:
          f:  the filter object with the method that does the filtering
        '''
        super(SupportFilterWidget,self).__init__(
            name = 'Genome Match',
        )
        self.tags = ['match', 0]
        self.filter = f

    def filter_cb(self):
        '''
        Callback activated when the slider changes.  Saves the new
        slider settings in the filter object.
        '''
        self.filter.matched = self.value


class TypeAMinFilterWidget(pn.widgets.TextInput):
    """
    Display a text entry box for the minimum homozygosity for Type A NCOs
    """

    def __init__(self, f):
        '''
        Arguments:
          f:  the NCOFilter object with the method that does the filtering
        '''
        c = Config()
        super(TypeAMinFilterWidget,self).__init__(
            name = 'Type A Minimum Homozygosity',
            width = 75,
        )
        # self.tags = ['match', 0]
        self.value = f'{c.post_min_z:0.2f}'
        self.filter = f

    def filter_cb(self):
        '''
        Callback activated when the text changes.  Saves the new
        value in the filter object.
        '''
        self.filter.min_z = float(self.value)


class TypeBDeltaFilterWidget(pn.widgets.TextInput):
    """
    Display a text entry box for the homozygosity range for Type B NCOs
    """

    def __init__(self, f):
        '''
        Arguments:
          f:  the NCOFilter object with the method that does the filtering
        '''
        c = Config()
        super(TypeBDeltaFilterWidget,self).__init__(
            name = 'Type B Homozygosity Range',
            width = 75,
        )
        # self.tags = ['match', 0]
        self.value = f'{c.post_delta_z:0.2f}'
        self.filter = f

    def filter_cb(self):
        '''
        Callback activated when the text changes.  Saves the new
        value in the filter object.
        '''
        self.filter.delta_z = float(self.value)

class NCOSizeFilterWidget(pn.widgets.TextInput):
    """
    Display a text entry box for the minimum number of SNPs in an NCO
    """

    def __init__(self, f):
        '''
        Arguments:
          f:  the NCOFilter object that will do the filtering.
        '''
        c = Config()
        super(NCOSizeFilterWidget,self).__init__(
            name = 'NCO Minimum Size (#SNPs)',
            width = 50,
        )
        self.value = str(c.post_block_size)
        # self.tags = ['size', 1]
        self.filter = f

    def filter_cb(self):
        '''
        Callback activated when the text changes.  Saves the new
        value in the filter object.
        '''
        self.filter.size = int(self.value)


class NCOCoverFilterWidget(pn.widgets.TextInput):
    """
    Display a text entry box for the minimum coverage in an NCO
    """

    def __init__(self, f):
        '''
        Arguments:
          f:  the NCOFilter object that will do the filtering.
        '''
        c = Config()
        super(NCOCoverFilterWidget,self).__init__(
            name = 'NCO Minimum Coverage',
            width=50,
        )
        self.value = str(c.post_min_cover)
        # self.tags = ['size', 1]
        self.filter = f

    def filter_cb(self):
        '''
        Callback activated when the text changes.  Saves the new
        value in the filter object.
        '''
        self.filter.min_cover = int(self.value)


class FilterBox(pn.Column):
    """
    A FilterBox is a column layout that has labels and an instance
    of each of the filter widgets.
    """

    def __init__(self, snp_filter):
        '''
        The `__init__` method instantiates the widgets that will
        filter SNPs and puts them in the Column
        in the order they will be displayed.

        Arguments:
          f:  the filter object with the method that does the filtering (passed to
              the constructors for each filter widget)
        '''
        super(FilterBox, self).__init__()

        self._widgets = []
        for cls in [BlockSizeFilterWidget, BlockLengthFilterWidget, CoverageFilterWidget, SupportFilterWidget]:
            w = cls(snp_filter)
            self._widgets.append(w)
            self.append(w)

    def widgets(self):
        '''
        Return a list of filter widgets.
        '''
        return self._widgets
        

class NCOFilterBox(pn.Column):
    """
    A FilterBox is a column layout that has labels and an instance
    of each of the filter widgets.
    """

    def __init__(self, nco_filter):
        '''
        The `__init__` method instantiates the widgets that will
        filter SNPs and puts them in the Column
        in the order they will be displayed.

        Arguments:
          f:  the filter object with the method that does the filtering (passed to
              the constructors for each filter widget)
        '''
        super(NCOFilterBox, self).__init__()

        self._widgets = []
       
        for cls in [TypeAMinFilterWidget, TypeBDeltaFilterWidget, NCOCoverFilterWidget, NCOSizeFilterWidget]:
            logging.debug(f'{cls}')
            w = cls(nco_filter)
            self._widgets.append(w)
            self.append(w)

    def widgets(self):
        '''
        Return a list of filter widgets.
        '''
        return self._widgets

class PeakViewerApp(pn.template.BootstrapTemplate):

    def __init__(self, **params):

        super(PeakViewerApp, self).__init__(**params)

        self.peaks = None               # groups of SNPs found by peak finder
        # self._filtered = None           # filtered blocks
        # self._ncos = None               # NCOs

        self.snp_filter = SNPFilter({})
        self.block_widgets = FilterBox(self.snp_filter)

        self.nco_filter = NCOFilter({})
        self.nco_widgets = NCOFilterBox(self.nco_filter)
        self.nco_widgets.visible = False

        self.nco_switch = pn.widgets.Switch(name='nco_switch')

        button_style_sheet = ''':host(.solid) .bk-btn {
            --color: white;
        }
        '''

        self.back_button = pn.widgets.Button(name='◀︎', stylesheets=[button_style_sheet])
        self.forward_button = pn.widgets.Button(name='▶︎', stylesheets=[button_style_sheet])
        self.chromosome_id = pn.widgets.TextInput(name="", value="")
        self.open_blocks = set()

        # self.chromosome_pattern = pn.widgets.TextInput(name="Chromosomes", value="BSP.*")
        # self.size_graph_button = pn.widgets.Button(name='Block Size', stylesheets=[button_style_sheet])
        # self.length_graph_button = pn.widgets.Button(name='Block Length', stylesheets=[button_style_sheet])
        # self.location_graph_button = pn.widgets.Button(name='Block Location', stylesheets=[button_style_sheet])

        # self.download_button = pn.widgets.FileDownload(callback=self.download_cb, filename='summary.csv', align='center', icon='download')
        # self.download_pane = pn.GridBox(self.download_button, height=200, width=SIDEBAR_WIDTH)
        # self.download_button.visible = False

        self.attach_callbacks()

        self.chromosome_panel = pn.Column(
            pn.pane.HTML('<h3>Chromosome Viewer</h3>'),
            pn.Row(self.back_button, self.chromosome_id, self.forward_button),
            pn.pane.HTML('<p>Placeholder</p>'),
            height=800,
        )

        self.sidebar.append(pn.pane.HTML("<h3>Block Parameters</h3>"))
        self.sidebar.append(self.block_widgets)
        self.sidebar.append(pn.layout.Divider())
        self.sidebar.append(pn.Row(
            pn.pane.HTML("<b>Find NCOs</b>"),
            self.nco_switch,
        ))
        self.sidebar.append(self.nco_widgets)

        self.main.append(self.chromosome_panel)

    def attach_callbacks(self):
        '''
        This method is called after all of the widgets have been created.
        It connects various widgets to functions that will be called when
        the widgets are activated.
        '''
        for w in self.block_widgets.widgets():
            w.param.watch(self.filter_cb, ['value'])

        for w in self.nco_widgets.widgets():
            w.param.watch(self.filter_cb, ['value'])

        for w in [self.back_button, self.forward_button]:
            w.param.watch(self.change_chromosome_cb, ['value'])

        # for w in [self.size_graph_button, self.length_graph_button, self.location_graph_button]:
        #     w.param.watch(self.summary_plot_cb, ['value'])

        self.chromosome_id.param.watch(self.chromosome_edited_cb, ['value'])      

        self.nco_switch.param.watch(self.find_ncos_cb, ['value'])

    def load_data(self, args):
        '''
        This method is called from the top level application after the GUI has
        been initialized. It reads the three data files needed by the application:
          * the interval data file has the names of all the chromosomes 
            (regardless of whether any SNPs were found); it's used to initialize the 
            list of chromosome names.  
          * the peaks data file (loaded by the filter object) has blocks of SNPs to
            visualize
          * the crossover data file has locations of main crossover points

        Chromosome names are saved in two instance vars:
          * clist is a list of all chromosome names; the forward and backward buttons
            move to the next or previous chromosome in this list
          * cmap is a dictionary that maps chromosome names to locations in clist; when
            the user types a new name in the chromosome name box the app uses the map
            to find the new list location

        Arguments:
          args:  command line arguments 
        '''
        c = Config()
        
        logging.info('Loading intervals')
        self.intervals = pd.read_pickle(c.intervals_default, compression='gzip').groupby('chrom_id')
        self.clist = list(self.intervals.groups.keys())
        self.cmap = { name: i for i, name in enumerate(self.clist)}
        logging.info(f'  read {len(self.intervals)} intervals')

        logging.info('Loading crossovers')
        self.crossovers = pd.read_pickle(c.crossovers_default, compression='gzip').groupby('chrom_id')
        logging.info(f'  read {len(self.crossovers)} crossovers')

        logging.info('Loading peaks')
        df = pd.read_csv(args.peaks)
        df['chr_length'] = df.chromosome.map(lambda n: chr_length[n])
        df['location'] = df.position / df.chr_length
        df['homozygosity'] = df.ref_reads / (df.ref_reads + df.var_reads)
        self.peaks = df.groupby('chrom_id')

        # self._chromosome_names = self.peaks.groups.keys()

        # p = Path(args.filtered)
        # if p.is_file():
        #     logging.info('Loading filtered')
        #     self._filtered = pd.read_csv(p)
        #     logging.info(f'  read {len(self._filtered)} filtered blocks from {args.filtered}')

        # p = Path(args.ncos)
        # if p.is_file():
        #     logging.info('Loading NCOs')
        #     self._ncos = pd.read_csv(p)
        #     logging.info(f'  read {len(self._ncos)} NCO blocks from {args.ncos}')

        # setting a value in the chromosome name widget triggers an update
        # to the graphic to display the first chromosome
        self.chr_index = 0
        self.chromosome_id.value = self.clist[self.chr_index]

    def display_chromosome(self):
        '''
        Update the chromosome display.  Called whenever the chromosome ID changes.  Shows
        a set of rectangular patches where the color is based on the region identified by 
        the HMM.  Below that is a grid with one row for each block of SNPs identifed by
        the peak finder.  Draw a vertical line at the crossover location (if there is one).
        '''
        chr_id = self.chromosome_id.value
        chrom = self.intervals.get_group(chr_id)
        logging.info(f'peak {chr_id}, {len(chrom)} SNPs')
    
        graphic = pn.Column()
        if chr_id in self.peaks.groups:
            df = self.peaks.get_group(chr_id)
            res, self.summary = self.snp_filter.apply(df)
            if self.nco_switch.value:
                res = self.nco_filter.apply(res)
                self.blocks = res.groupby('blk_id')
                self.nco_blocks = { n for n, grp in self.blocks if grp.nco.max() == 2 }
            else:
                self.blocks = res.groupby('blk_id')
                self.nco_blocks = set()
            grid = self._make_grid()
            graphic.append(grid)

        rects = PatchCollection(self._make_patches(chrom, chr_id), match_original=True)
        fig, ax = plt.subplots(figsize=(12,1))
        plt.box(False)
        plt.yticks([])
        plt.xticks(ticks=np.linspace(0,20000000,5), labels=[f'{int(n*20)}Mbp' for n in np.linspace(0,1,5)])
        ax.xaxis.set_ticks_position('top')
        ax.add_collection(rects)
        if chr_id in self.crossovers.groups:
            for _, xo in self.crossovers.get_group(chr_id).iterrows():
                logging.info(f'crossover {xo.start} {xo.is_CO}')
                if xo.is_CO:
                    plt.axvline(xo.start, color='green')
        plt.xlim(0,20000000)
        plt.ylim(0,2000000)
        plt.close(fig)
        graphic.insert(0,pn.pane.Matplotlib(fig, dpi=72, tight=True))

        self.chromosome_panel.pop(-1)
        self.chromosome_panel.append(graphic)

    def _make_patches(self, df, chr_id):
        '''
        Create a horizontal bar as a collection of rectangular patches, with one
        patch for each row in the data frame.  The rows have the starting coordinates
        lengths, and HMM states of chromosome regions, used to define the width and
        color of a patch.  Add black dots to indicate the locations of blocks.
        '''
        pcolor = {
            'CB4856': 'dodgerblue',
            'N2': 'indianred'
        }
        res = []
        for _, r in df.iterrows():
            c = pcolor.get(r.hmm_state) or 'lightgray'
            res.append(Rectangle((r.start,500000), r.length, 1000000, color=c))
        if chr_id in self.peaks.groups:
            for blk_index, _ in self.summary.iterrows():
                _, blk_id = blk_index
                block = self.blocks.get_group(blk_id)
                x0 = block.iloc[0].position
                res.append(Circle((x0,750000), 50000, color='black'))
                if blk_id in self.nco_blocks:
                    res.append(Circle((x0,1200000), 50000, color='yellow'))
        return res

    def _make_grid(self):
        '''
        Make a Column object that has a collection of figures, one for each block 
        in a chromosome (saved in an instance var).  The figures are saved in a grid.
        Below each figure is a text widget containing the data frame with the filtered 
        SNPs in the block, i.e. the grid for a chromosome with N blocks as 2*N rows.
        The frames are initially hidden.  The rows that have figures also have a toggle
        button; clicking this button will show or hide the frame.
        '''
        pcolor = {
            'CB4856': 'dodgerblue',
            'N2': 'indianred',
            'uCB4856': 'lightsteelblue',
            'uN2': 'lightpink',
            'unknown': 'lightgray',
            'het': 'palegoldenrod',
        }
        nco_sym = {
            0: ' ',
            1: '✓',
            2: '★',
        }
        self.block_buttons = {}
        self.block_text = {}
        g = pn.Column()
        for blk_index, blk_stats in self.summary.iterrows():
            _, blk_id = blk_index
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
            self.block_buttons[blk_id].on_click(self.toggle_text_cb)
            cols = ['position','base_geno','hmm_state1','reference','ref_reads','variant','var_reads','background','homozygosity']
            df = block[cols]
            if self.nco_switch.value:
                # df['nco'] = block.nco.map(lambda n: nco_sym[n])
                df = pd.concat([df, block.nco.map(lambda n: nco_sym[n])], axis=1)
            self.block_text[blk_id] = pn.pane.DataFrame(df, visible=False)
            if blk_id in self.open_blocks:
                self.block_text[blk_id].visible = True
                self.block_buttons[blk_id].name = '∨'
            g.append(pn.Row(
                self.block_buttons[blk_id],
                pn.pane.Matplotlib(fig, dpi=72, tight=True),
                styles={'background':'WhiteSmoke'},
            ))
            g.append(pn.Row(self.block_text[blk_id]))
        return g
    
    def toggle_text_cb(self, e):
        '''
        Callback function invoked when a toggle button in the chromosome display is clicked.
        Toggles the visibility of the frame and updates the button name based on the new
        visibility state.
        '''
        i = e.obj.tags[0]
        self.block_text[i].visible = not self.block_text[i].visible
        if self.block_text[i].visible:
            self.block_buttons[i].name = '∨'
            self.open_blocks.add(i)
        else:
            self.block_buttons[i].name = '>'
            self.open_blocks.remove(i)

    def filter_cb(self, e):
        '''
        Callback function invoked when any of the widgets used for filtering (sliders,
        checkbox, etc) is activated.
        '''
        e.obj.filter_cb()
        self.display_chromosome()       

    def find_ncos_cb(self, e):
        '''
        Callback function invoked when the 'Find NCOs' switch is toggled
        ''' 
        self.nco_widgets.visible = not self.nco_widgets.visible
        self.display_chromosome()

    def change_chromosome_cb(self, e):
        '''
        Callback function invoked when the left or right button next to the chromosome
        name is clicked.
        '''
        delta = 1 if e.obj is self.forward_button else -1
        self.chr_index = (self.chr_index + delta) % len(self.clist)
        self.chromosome_id.value = self.clist[self.chr_index]
        self.open_blocks = set()

    def chromosome_edited_cb(self, e):
        '''
        Callback function invoked whenever the chromosome ID is updated.  This will
        happen when the user edits the chromosome name or the name changes after
        a button click.
        '''
        idx = self.cmap.get(e.obj.value)
        if idx is not None:
            self.chr_index = idx
            self.chromosome_id.value = self.clist[idx]
            self.open_blocks = set()
            self.display_chromosome()

    # histogram_params = {
    #     'Block Size': {
    #         'col':     'blk_size',
    #         'title':   'Block Size',
    #         'xlabel':  'Number of SNPs',
    #         'ylabel':  'Number of Blocks',
    #         'hist': {
    #             'bins': 10,
    #             'rwidth': 0.8,
    #             'align': 'left',
    #             'range': (1,100),
    #         },
    #     },
    #     'Block Length': {
    #         'col':     'blk_len',
    #         'title':   'Block Length',
    #         'xlabel':  'Length (bp)',
    #         'ylabel':  'Number of Blocks',
    #         'hist': {
    #             'bins': 10,
    #             'rwidth': 0.8,
    #         },
    #     },
    #     'Block Location': {
    #         'col':     'blk_loc',
    #         'title':   'Block Location',
    #         'xlabel':  'Relative Position in the Chromosome',
    #         'ylabel':  'Number of Blocks',
    #         'hist': {
    #             'bins': 100,
    #             'range': (0,1),
    #         },
    #     },
    # }

    # def summary_plot_cb(self, e):
    #     '''
    #     Callback function invoked when the user clicks the name of one of the 
    #     histograms in the summary tab.  All histograms have the same basic parameters,
    #     those that are specific to a type of data are defined in the
    #     `histogram_params` dictionary.
    #     '''
    #     params = self.histogram_params[e.obj.name]
    #     self.tabs[1].loading = True
    #     # self.filter.set_chromosome(self.chromosome_pattern.value)
    #     self.filter.chromosome = self.chromosome_pattern.value
    #     self.summary_df = self.filter.summary()
    #     fig, ax = plt.subplots(figsize=(7,5))
    #     plt.hist(self.summary_df[params['col']], label=self.chromosome_pattern.value, **params['hist'])
    #     plt.title(params['title'])
    #     plt.xlabel(params['xlabel'])
    #     plt.ylabel(params['ylabel'])
    #     plt.legend(handlelength=0)
    #     plt.close(fig)
    #     self.tabs[1].loading = False
    #     self.tabs[1].pop(-1)
    #     self.tabs[1].append(pn.pane.Matplotlib(fig, dpi=72, tight=True))
    #     self.download_button.visible = True

    # def download_cb(self):
    #     '''
    #     Callback function invoked when the user clicks the download button (made 
    #     visible after plotting a histogram).
    #     '''
    #     sio = io.StringIO()
    #     self.summary_df.to_csv(sio)
    #     sio.seek(0)
    #     return sio


def make_app(args):
    """
    Instantiate the top level widget, load data from the data files.

    Arguments:
      args: command line arguments (including names of data files)

    Returns:
        a PeakViewerApp object
    """
    c = Config()

    app = PeakViewerApp(
        title='NCO Explorer', 
        sidebar_width=c.sidebar_width,
    )
    app.load_data(args)
    if args.log != 'quiet':
        logging.getLogger('bokeh').setLevel(logging.ERROR)
    return app

def start_app(args):
    """
    Main entry point, called from the top level `xo` script.

    Initialize the Panel library, instantiate the app, and pass
    the app to the server.
    """
    pn.extension(design='native')
    pn.config.throttled = True
    try:
        app = make_app(args)
        pn.serve( 
            app,
            port = args.port,
            verbose = True,
            autoreload = True,
            websocket_origin= '*',
        )
    except Exception as err:
        logging.exception(err)


