# Peak Finder

This script scans a chromosome to look for sequences of SNPs that are potentially associated with crossover or non-crossover events.

## Searching for Blocks of SNPs

The name "peak finder" comes from the method we use to look for these events.  We define a **signal function** that is based on counting the number of SNPs.  At any location in the chromosome, let *f* be the total number of CB4856 SNPs seen up to that point, and let *g* be the total number of N2 SNPs seen to that point.  Then our signal function is simply *f* − *g*.  As we move across the chromosome, the signal goes up whenever we encounter a CB4856 SNP and goes down whenever we see an N2 SNP.

As an example, suppose a chromosome has a single crossover point: before the location of the point all the SNPs are from the CB4856 parent, and after this location all the SNPs are from the N2 parent. The signal will increase monotonically up to location of the point and then decrease monotonically.  The signal has one peak, at the location where the crossover occurred.

### Non-Crossovers

If there are non-crossover events, we expect to see small regions where SNPs from one parent occur somewhere in the middle of a long sequence from the other parent.  These will appear as small "blips" in the signal:  a short drop in the region that is otherwise increasing, or a short rise in the region that is otherwise decreasing.

We use the `find_peaks` function from the NumPy signal processing library to look for these regions in our signal function.  We define a **block** to be the set of SNPs from the start to the end of one of the regions identified by the peak finder.

### Look for Valleys as Well as Peaks

Some chromosomes will have N2 SNPs before the main crossover point and CB4856 SNPs after that.  In this case the signal decreases before the crossover and increases after that point, so the crossover is at a location where there is a "valley" and not a peak.  The NumPy function returns locations of both peaks and valleys, and our code saves locations from each type of "blip".

### Maximum Block Size

If a chromosome has a single crossover, the peak finding algorithm will think the sequence from the start of the chromosome up to the crossover location is a block, and the sequence from the crossover point to the end of the chromosome is another block.  

To prevent these two blocks from being included in the output the script has a block size limit, which is the maximum number of SNPs that can be included in a block.

There are also "pathological" cases where there are small blips inside a longer region that also stands out from the background.  In these situations the peak finder reports the coordinates of the longer enclosing block along with the locations of the smaller blocks nested inside.  Choosing a smaller block size limit will prevent the enclosing block from being included in the output while still including the smaller blocks.

### No Crossovers

Some chromosomes in the data set have no crossover points.  In these situations the peak finder will not find a peak or a valley and the output file will have no blocks for this chromosome.

## Shell Command

To run the peak finder with the default options type

```bash
$ xo peaks
```

That will read the SNP data from TIGER and write blocks to an output file in CSV format.

- The default input file name is `BSP_TIGER.marker_dataframe.pickle.gzip`.  A different file name can be specified with the `--snps` option.
- The default output file name is `peaks.csv`.  A different name can be specified with `--output`.
- The default block size limit is 1000 SNPs.  A different value can be specified with the `--max_snps` option.

## Output File

The data written by this command is saved in a CSV file.

Each row in the file is a SNP that occurs in a block identified by the peak finder.  The lines are exact copies of the lines from the original input file, including the original SNP ID in the first column.

One new column has been appended to each line.  This column is a block ID number, which allows the GUI and visualization commands to group SNPs by blocks.





