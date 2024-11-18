#
# This module defines a SNPFilter class used by the command line and
# GUI of the crossover explorer application.
#
# John Conery
# University of Oregon
#

import logging
import pandas as pd
import re

class SNPFilter:
    """
    A SNPFilter applies filtering criteria to SNPs in blocks identified by the
    peak finder.  The constructor uses command line arguments to initialize 
    filtering critera.  Call a method named `apply` to filter a data set.

    Attributes:
      chromosome:  when filtering return SNPs in chromosomes with names that match this pattern
      size_range:  a pair of integers with the minimum and maximum number of SNPs in a block
      length_range:  a pair of integers with the minimum and maxium block length (in base pairs)
      coverage:  a minimum number of reads required to use a SNP
      matched:  True if a SNP's base genome must be the same as the HMM state
    """

    def __init__(self, args):
        self._chromosome = 'BSP.*'
        self._min_size = 1
        self._max_size = 100
        self._min_length = 1
        self._max_length = 10000
        self._coverage = 0
        self._matched = False

        filter_params = {
            'chromosomes':  'chromosome',
            'size':         'size_range',
            'length':       'length_range',
            'coverage':     'coverage',
            'match':        'matched',
        }
        for arg, attr in filter_params.items():
            if val := vars(args).get(arg):
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
        logging.info(f'Filtering {len(df)} SNPs')
        df = df[df.chrom_id.map(lambda s: bool(re.match(self._chromosome,s)))]
        logging.info(f'  chromosome:  {len(df)} match {self._chromosome}')

        if self._matched:
            df = df[df.base_geno == df.hmm_state1]
            logging.info(f'  match: {len(df)} have genome match')

        if self._coverage:
            df = df[df.var_reads + df.ref_reads > self._coverage]
            logging.info(f'  coverage:  {len(df)} have reads > {self._coverage}')
               
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
