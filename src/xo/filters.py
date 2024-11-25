#
# This module defines two classes used by the command line and
# GUI of the crossover explorer application:
#
# A SNPFilter processes the raw blocks produced by the peak finder,
# filtering by block size, block length, etc.
#
# An NCOFilter does furter postprocessing of the filtered blocks
# based on homozygosity and coverage.
#
# John Conery
# University of Oregon
#

import logging
import numpy as np
import pandas as pd
import re

from rich.console import Console
from rich.table import Table

from .config import Config

class SNPFilter:
    """
    A SNPFilter applies filtering criteria to SNPs in blocks identified by the
    peak finder. Pass the constructor a dictionary that has settings filtering 
    critera. Call a method named `apply` to filter a data set.

    Attributes:
      chromosome:  when filtering return SNPs in chromosomes with names that match this pattern
      size_range:  a pair of integers with the minimum and maximum number of SNPs in a block
      length_range:  a pair of integers with the minimum and maxium block length (in base pairs)
      coverage:  a minimum number of reads required to use a SNP
      matched:  True if a SNP's base genome must be the same as the HMM state
    """

    def __init__(self, args):
        c = Config()
        self._chromosome = c.filter_chromosome_pattern
        self._min_size, self._max_size = c.filter_block_size
        self._min_length, self._max_length = c.filter_block_length
        self._coverage = c.filter_coverage
        self._matched = False
        self._result = None
        self._summary = None

        filter_params = {
            'chromosomes':  'chromosome',
            'size':         'size_range',
            'length':       'length_range',
            'coverage':     'coverage',
            'match':        'matched',
        }
        for arg, attr in filter_params.items():
            if val := args.get(arg):
                setattr(self, attr, val)

    def __repr__(self):
        res = f'chr {self._chromosome}'
        res += f' size {self._min_size}..{self._max_size}'
        res += f' len {self._min_length}..{self._max_length}'
        res += f' cover {self._coverage}'
        res += f' match {self._matched}'
        return res
  
    @property
    def chromosome(self):
        '''Filtered SNPs must have names that match this pattern.'''
        return self._chromosome
    
    @chromosome.setter
    def chromosome(self, s):
        self._chromosome = s

    @property
    def size_range(self):
        '''Minimum and maximum block size (number of SNPs)'''
        return (self._min_size, self._max_size)
    
    @size_range.setter
    def size_range(self, limits):
        self._min_size = limits[0]
        self._max_size = limits[1]

    @property
    def length_range(self):
        '''Minimum and maximum block length (bp)'''
        return (self._min_length, self._max_length)
    
    @length_range.setter
    def length_range(self, limits):
        self._min_length = limits[0]
        self._max_length = limits[1]

    @property
    def matched(self):
        '''If True SNPs must have matching HMM state'''
        return self._matched
    
    @matched.setter
    def matched(self, p):
        self._matched = p

    @property
    def coverage(self):
        '''Minimum number of reads for a SNP'''
        return self._coverage

    @coverage.setter    
    def coverage(self, n):
        self._coverage = n

    @property
    def result(self):
        '''Filtered data'''
        return self._result

    @property
    def summary(self):
        '''Summary of filtered data'''
        return self._summary

    def apply(self, df):
        '''
        Apply the filtering criteria to a set of SNPs.  The first step is to make a 
        subset of SNPs that pass the chromosome name, coverage, and match filters.  
        The resulting frame is then grouped by block.

        The second step makes a separate summary frame based on size and location of
        each block.  This frame is then filtered using the block size and block length
        filters.

        Arguments:
          df:  a data frame containing SNPs (including block IDs)

        Returns:
          res:  a frame containing the filtered data
          summary:  a data frame with size, length, and location of each block
        '''
        logging.info(f'SNPFilter.apply {len(df)} SNPs')
        df = df[df.chrom_id.map(lambda s: bool(re.match(self._chromosome,s)))]
        logging.info(f'  chromosome:  {len(df)} match {self._chromosome}')

        if self._matched:
            df = df[df.base_geno == df.hmm_state1]
            logging.info(f'  match: {len(df)} have genome match')

        if self._coverage:
            df = df[df.var_reads + df.ref_reads >= self._coverage]
            logging.info(f'  coverage:  {len(df)} have reads ≥ {self._coverage}')
               
        groups = df.groupby(['chrom_id','blk_id'])
        logging.info(f'Groupd into {len(groups)} blocks')

        sf = pd.concat(
            [
                groups.size().rename('blk_size'), 
                (groups.max('position') - groups.min('position')).position.rename('blk_len'), 
                groups.mean('location').location.rename('blk_loc'),
            ],
            axis=1
        )

        min_size = sf.blk_size >= self._min_size
        max_size = sf.blk_size <= self._max_size
        min_len = sf.blk_len >= self._min_length
        max_len = sf.blk_len <= self._max_length

        sf = sf[min_size & max_size & min_len & max_len]
        logging.info(f'Summary has {sf.blk_size.sum()} SNPs in {len(sf)} blocks')

        res = pd.concat(groups.get_group(n) for n in sf.index)

        return res, sf


class NCOFilter:
    """
    An NCOFilter applies additional filtering criteria to blocks of SNPs.  
    Pass the constructor a dictionary that has settings filtering critera.  
    Call a method named `apply` to filter a data set.

    Attributes:
      min_z:      homozygosity for type A blocks
      delta_z:    homozygosity for type B blocks
      min_cover:  minumum number of reads of each type
      size:     minimum block size
    """

    def __init__(self, args):
        c = Config()
        self._min_z = c.post_min_z
        self._delta_z = c.post_delta_z
        self._min_cover = c.post_min_cover
        self._size = c.post_block_size
        self._result = None

        filter_params = ['min_z','delta_z', 'min_cover', 'size']

        for attr in filter_params:
            if val := args.get(attr):
                logging.debug(f'{attr} {val}')
                setattr(self, attr, val)

    def __repr__(self):
        res = f' z {self._min_z}'
        res += f' ∆ {self._delta_z}'
        res += f' snps {self._min_cover}'
        res += f' len {self._size}'
        return res

    @property
    def min_z(self):
        '''Minimum homozygosity for Type A blocks'''
        return self._min_z
    
    @min_z.setter
    def min_z(self, x):
        self._min_z = x

    @property
    def delta_z(self):
        '''Homozygosity range for Type B blocks'''
        return self._delta_z
    
    @delta_z.setter
    def delta_z(self, x):
        self._delta_z = x

    @property
    def min_cover(self):
        '''Minimum number of SNPs of each type'''
        return self._min_cover
    
    @min_cover.setter
    def min_cover(self, n):
        self._min_cover = n

    @property
    def size(self):
        '''Minimum block size'''
        return self._size
    
    @size.setter
    def size(self, n):
        self._size = n

    @property
    def result(self):
        '''Filtered data'''
        return self._result

    def apply(self, df):
        '''
        Apply the NCO filtering criteria to a set of SNPs in filtered blocks SNPs.

        Arguments:
          df:  a data frame of filtered SNPs

        Returns:
          res:  a frame containing the filtered data
        '''
        logging.info(f'NCOFilter.appy {len(df)} SNPs')

        res = []
        for name, block in df.groupby(['chrom_id','blk_id']):
            logging.debug(f'{name} {len(block)} SNPs')
            block['nco'] = self._scan(block)
            res.append(block)

        self._result = pd.concat(res)
        return self._result


    def _scan(self, block):
        '''
        Scan a block for runs that have satisfy homozygosity requirements
        and have a minimum number of SNPs.  Add a column to a block to
        show which SNPs might be part of an NCO.
        '''

        def next_interval():
            nonlocal i, j
            j = i
            logging.debug(f'next_interval {i} {j}')
            while (j < len(p)):
                logging.debug(f'{i} {j} {p.iloc[j]}')
                if p.iloc[j]:
                    logging.debug(f'extend {j}')
                    j += 1
                elif j - i >= self.size:
                    logging.debug(f'long enough {i} {j}')
                    break
                else:
                    logging.debug(f'restart')
                    j += 1
                    i = j
            logging.debug(f'exit {i} {j} {self.size}')
        
        # block['nco'] = 0
        col = np.zeros(len(block), dtype=np.int8)

        z = block.homozygosity
        if block.iloc[0].background == 'CB4856':
            rr = block.ref_reads
            p = (z >= self.min_z) & (rr >= self.min_cover)
        else:
            rr = block.ref_reads
            vr = block.var_reads
            p = (abs(z - 0.5) <= self.delta_z) & (rr >= self.min_cover) & (vr >= self.min_cover)

        col[p] = 1
        
        res = []
        i = j = 0
        while i < len(block):
            next_interval()
            if (j - i) >= self.size:
                res.append((i,j))
            i = j + 1
            if len(res) > 3:
                break
        logging.debug(f'res {res}')

        if res:
            res.sort(key=lambda p: p[1]-p[0])
            i, j = res[-1]
            logging.debug(f'nco: {i}..{j}')
            col[i:j] = 2

        return col

    def print_summary(self):
        ncos = self._result.groupby(['chrom_id','blk_id'])
        counts = ncos.SNP.count()
        sf = counts.groupby(level='chrom_id').count()

        c = Console()
        g = Table.grid(padding=[0,2,0,1])
        g.add_column()
        g.add_column(justify='right')

        g.add_row('Chromosomes with NCOs', str(len(sf)))
        g.add_row('Mean number of NCOs per chromosome', f'{sf.mean():0.1f}')
        c.print(g)

        c.print()
        g = Table.grid(padding=[0,2,0,1])
        g.add_column()
        g.add_column(justify='right')
        for chr, n in sf.items():
            g.add_row(chr, str(n))
        c.print(g)
