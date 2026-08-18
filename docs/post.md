# Identify NCOs (`xo post`)

The `post` command reads blocks saved by the `filter` command, searches through each block to see if it contains a noncrossover, and writes an updated data set to the output file.

This command will look for contiguous groups of SNPs that satisfy our criteria for being parts of noncrossovers.
If a block has several such groups only the longest one is marked as being the noncrossover event in the block.

The new CSV file is the same as the input, but it contains one new column named `nco`:
- a 0 means the SNP is not part of a noncrossover
- a 1 means the SNP satisfies all the criteria for being part of a crossover but was not included (either because it's in a group that is too small or there is another longer group elsewhere in the block)
- a 2 means the SNP is part of the noncrossover found in this blocl

The `--help` option shows the various options that can be given on the command line:

```
$ xo post --help
usage: xo post [-h] [--blocks F] [--output F] [--min_z N] [--delta_z N] [--min_cover N] [--size N]

options:
  -h, --help     show this help message and exit
  --blocks F     file with filtered blocks
  --output F     output file
  --min_z N      homozygosity for Type A blocks
  --delta_z N    homozygosity for Type B blocks
  --min_cover N  minimum number of reads of each type
  --size N       minimum block size (bp)
```

The `--blocks` and `--output` options specify the names of the input and output files.
The defaults are `filtered.csv` for the input and `ncos.csv` for the output.

The remaining options all correspond to "widgets" in the "Find NCOs" controls in the GUI.

| Widget | Command Line | Default |
| --- | --- | --- |
| Type A Minimum Homozygosity | `--min_z` | `--min_z 0.9` |
| Type B Homozygosity Range | `--delta_z` | `--delta_z 0.1` |
| NCO Minimum Coverage | `--min_cover` | `--min_cover 2` |
| NCO Minimum Size | `--size` | `--size 5` |

If an option is not specified the same default values used in the GUI are used here by the script.
