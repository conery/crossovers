#
# This module defines a SNPFilter class used by the command line and
# GUI of the crossover explorer application.
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
    """
    A SNPFilter object is the main interface to the SNP data.  After creating the object
    call the `load_data` method to read a data frame with SNP calls produced by TIGER.
    Set filtering parameters (block length, _etc_), then call the `filter` method to 
    apply all the filters.

    Attributes:
      snps:  the input data frame
      chromosome:  when filtering return SNPs in chromosomes with names that match this pattern
      size_range:  a pair of integers with the minimum and maximum number of SNPs in a block
      length_range:  a pair of integers with the minimum and maxium block length (in base pairs)
      coverage:  a minimum number of reads required to use a SNP
      matched:  True if a SNP's base genome must be the same as the HMM state
    """

    def __init__(self):
        self._snps = None
        self._chromosome = 'BSP.*'
        self._size_range = (1,100)
        self._length_range = (1, 10000)
        self._coverage = 0
        self._matched = False
        self._groups = None

    def load_data(self, fn):
        '''
        Read the SNP data from a CSV file.  Add two new columns (chromosome length and 
        relative SNP location) used in summaries. SNPs are grouped by chromosome ID and
        the groups are saved in an instance variable.

        Arguments:
          fn: name of CSV file
        '''
        self._snps = pd.read_csv(fn).rename(columns={'Unnamed: 0': 'SNP'})
        self._snps['chr_length'] = self._snps.chromosome.map(lambda n: chr_length[n])
        self._snps['location'] = self._snps.position / self._snps.chr_length
        self._groups = self._snps.groupby('chrom_id')

    def has_chromosome_block(self, chr_id):
        '''
        Return True if the SNP data has blocks for a chromosome.

        Arguments:
          chr_id:  the name of the chromosome to look for.
        '''
        return chr_id in self._groups.groups
    
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
        return self._size_range
    
    @size_range.setter
    def size_range(self, limits):
        self._size_range = limits

    @property
    def length_range(self):
        '''"Minimum and maximum block length (bp)"'''
        return self._length_range
    
    @length_range.setter
    def length_range(self, limits):
        self._length_range = limits

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

    def apply(self, chr_id):
        '''
        Apply the match and coverage criteria to a single chromosome, return the 
        filtered SNPs grouped by blocks and a summary table.

        Arguments:
          chr_id:  the ID of the chromosome to use

        Returns:
          groups:  the blocks in the filtered chromosome
          summary:  a data frame with size, length, and location of each block
        '''
        df = self._groups.get_group(chr_id)
        if self._matched:
            df = df[df.base_geno == df.hmm_state1]
        if self._coverage:
            df = df[df.var_reads + df.ref_reads > self._coverage]
        return df.groupby('blk_id'), self.summary(df)
    
    def summary(self, df = None):
        '''
        Create a summary table (a data frame with columns for size, length, and location
        of blocks in a data set) with descriptions of the groups that satisfy the block
        size and length filters.  When called from the GUI, `df` is a frame with SNPs from
        the current chromosome, and the result has blocks from only this chromosome.  When
        called from the command line, `df` is None and the method uses the complete data set
        and returns a summary indexed by chromosome and block.
        '''
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
    
