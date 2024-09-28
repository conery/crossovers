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

# Chromosome lengths

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
        self._min_size = 1
        self._max_size = 100
        self._min_length = 1
        self._max_length = 10000
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

    def apply(self, chr_id = None):
        '''
        Apply the match and coverage criteria to a single chromosome, return the 
        filtered SNPs grouped by blocks and a summary table.

        If a chromosome ID is passed apply the filters to that chromosome only,
        otherwise use the complete set of all SNPs.

        The first step is to make a subset of SNPs that pass the coverage and match
        filters.  The resulting frame is then grouped by block.

        The second step makes a separate summary frame based on size and location of
        each block.  This frame is then filtered using the block size and block length
        filters.

        Arguments:
          chr_id:  the ID of the chromosome to use

        Returns:
          groups:  the blocks in the filtered chromosome
          summary:  a data frame with size, length, and location of each block
        '''
        if chr_id is None:
            df = self._snps[self._snps.chrom_id.map(lambda s: bool(re.match(self._chromosome,s)))]
        else:
            df = self._groups.get_group(chr_id)

        logging.info(f'Filtering {len(df)} SNPs')

        if self._matched:
            df = df[df.base_geno == df.hmm_state1]
            logging.info(f'{len(df)} match')

        if self._coverage:
            df = df[df.var_reads + df.ref_reads > self._coverage]
            logging.info(f'{len(df)} have coverage > {self._coverage}')
   
        if chr_id is None:
            groups = df.groupby(['chrom_id', 'blk_id'])
        else:
            groups = df.groupby('blk_id')
        logging.info(f'{len(groups)} groups')

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
        logging.info(f'summary has {sf.blk_size.sum()} SNPs in {len(sf)} blocks')

        return groups, sf
    