#
#
# Filter blocks of SNPs
#
# John Conery
# University of Oregon
#

import pandas as pd
import re

chr_length = {
    1: 15114068,
    2: 15311845,
    3: 13819453,
    4: 17493838,
    5: 20953657,
    6: 17739129,
}

class SNPFilter:

    def __init__(self):
        self._snps = None
        self._chromosome = 'BSP.*'
        self._size_range = (1,100)
        self._length_range = (1, 10000)
        self._coverage = 0
        self._matched = False
        # self._chr_names = []
        # self._chr_map = {}
        self._groups = None

    def load_data(self, fn):
        self._snps = pd.read_csv(fn).rename(columns={'Unnamed: 0': 'SNP'})
        self._snps['chr_length'] = self._snps.chromosome.map(lambda n: chr_length[n])
        self._snps['location'] = self._snps.position / self._snps.chr_length
        self._groups = self._snps.groupby('chrom_id')

    def has_chromosome_block(self, chr_id):
        return chr_id in self._groups.groups

    def set_chromosome(self, s):
        self._chromosome = s

    def set_size(self, limits):
        self._size_range = limits

    def set_length(self, limits):
        self._length_range = limits

    def set_matched(self, p):
        self._matched = p

    def set_coverage(self, n):
        self._coverage = n

    def apply(self, chr_id):
        '''
        Assumes caller has already checked to see if the chr_id is a block name
        '''
        df = self._groups.get_group(chr_id)
        if self._matched:
            df = df[df.base_geno == df.hmm_state1]
        if self._coverage:
            df = df[df.var_reads + df.ref_reads > self._coverage]
        return df.groupby('blk_id'), self.summary(df)
    
    def summary(self, df = None):
        if df is None:
            df = self._snps[self._snps.chrom_id.map(lambda s: bool(re.match(self._chromosome,s)))]
            groups = df.groupby(['chrom_id', 'blk_id'])
        else:
            groups = df.groupby('blk_id')

        data = pd.concat(
            [
                groups.size().rename('blk_size'), 
                (groups.max('position') - groups.min('position')).position.rename('blk_len'), 
                groups.mean('location').location.rename('blk_loc'),
            ],
            axis=1
        )
        min_size = data.blk_size >= self._size_range[0]
        max_size = data.blk_size <= self._size_range[1]
        min_len = data.blk_len >= self._length_range[0]
        max_len = data.blk_len <= self._length_range[1]
        return data[min_size & max_size & min_len & max_len]
    
