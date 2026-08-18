# Filter Blocks (`xo filter`)

The `filter` command reads a CSV file containing descriptions of blocks identified by the peak finder, applies a series of filters, and writes a new set of blocks to an output file.

The `--help` option shows the various options that can be given on the command line:

```
$ xo filter --help
usage: xo filter [-h] [--peaks F] [--crossovers F] [--output F] [--chromosomes P] [--size N N] [--length N N] [--coverage N] [--match]

options:
  -h, --help       show this help message and exit
  --peaks F        blocks saved by peaks.py
  --output F       output file
  --crossovers F   file with crossover locations
  --chromosomes P  chromosome name pattern
  --size N N       block size range (#SNPs)
  --length N N     block length range (bp)
  --coverage N     minimum coverage
  --match
```

The `--peaks` and `--output` options specify the names of the input and output files.
The defaults are `peaks.csv` for the input and `filtered.csv` for the output.

The `--crossovers` option specifies the name of a file with crossover locations on each chromosome.
The default is `BSP_COs_final_set.pickle.gzip` (based on the output of the HMM that calls SNPs).

The `--chromosomes` option allows you to choose which chromosomes to include in the output.
By default all chromosomes are included, but this option can narrow the set to a specific group.
The value of this arguments is a pattern, and only chromosomes with names that match this pattern will be included.
For example, `BSP-OR.*` is a pattern that means "all chromosomes with names that start with "BSP-OR".
See the section on Chromsome Names below for more examples.

The remaining options all correspond to "widgets" in the control panel section of the GUI.

| Widget | Command Line | Default |
| --- | --- | --- |
| Slider labeled "Block Size" | `--size` | `--size 0 100` |
| Slider labeled "Block Length" | `--length` | `--length 0 10000` |
| Slider labeled "Minimum Coverage" | `--coverage` | `--coverage 0` |
| Checkbox labeled "Genome Match" | `--match` | |

If an option is not specified the same default values used in the GUI are used here by the script.

**Note:**  the `--size` and `--length` options both requre two values, since they correspond to sliders where both the minimum and maximum can be selected.

## Chromosome Names

The default setting for chromosome names will generate output files using all of the chromosomes.However, if you want to make separate the data into groups you can enter a name pattern using regular expression syntax as the argument to the `--chromosomes` option.

The default value `BSP.*`, meaning "any chromosome with a name that starts with BSP" (in other words, all chromosomes).  Some other examples of name patterns:

| pattern       | chromosomes used                                             |
| ------------- | ------------------------------------------------------------ |
| `BSP-OR.*`    | all oocytes                                                  |
| `BSP-SR.*`    | all spermatocytes                                            |
| `BSP-OR-10.*` | the 10 worms with names BSP-OR-001, BSP-OR-002, ... BSP-OR-009. |
| `BSP-SR-.*-1` | chromosome 1 for all spermatocytes                           |
