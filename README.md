# Crossovers

A command line application to view crossover and non-crossover events in chromosomes.

### Installation

This command will install a Python script named `xo`:

```
pip3 install git+https://github.com/conery/crossovers.git
```

To make sure it's installed:

```bash
$ xo --help
usage: xo [-h] {peaks,view} ...
...
```

### Usage

To run the script type `xo CMND` where `CMND` is one of the operations to perform.  In the current version there are only two operations:

- `peaks` will use the SciPy peak finding function to look for potentially interesting locations in chromosome
- `view` will display a GUI to allow users to browse the output of the `peaks` command

### Find Interesting Blocks of SNPs

To list all the options for this step:

```bash
$ xo peaks --help
usage: xo peaks [-h] [--snps F] [--output F] [--max_snps N]

options:
  -h, --help    show this help message and exit
  --snps F      input (IGER marker) file
  --output F    output file
  --max_snps N  max number of SNPs in a block
```

The `--snps` option is the name of the file with the SNP data.  The default is a (pickled) dataframe generated by TIGER:

```
BSP_TIGER.marker_dataframe.pickle.gzip
```

The `--output` option specifies the name of the output file, which will be a plain text CSV file (default:  `peaks.csv`).

The `--max_snps` option is a cutoff for the maximum block size to write to the output (default: 1000).

### View SNPs

The `view` command needs two data files.  One is the CSV file with blocks of SNPs from the `peaks` command, the other is a summary of the SNPs in another pickled dataframe.

- specify the path to the peak data with `--peaks` option; the default is `peaks.csv` (the default output name from the `peaks` command)
- specify the path the summarized data with `--intervals`; the default is `BSP_TIGER.intervals_dataframe.pickle.gzip`

### Example

The simplest workflow is to create a directory and copy (or link) the two pickled dataframes, using their default file names.

Run the `peaks` command with default options:

```bash
$ xo peaks
```

It will take a few minutes, but the script will show a status line to indicate it's still working.  When it's done there will be a file named `peaks.csv` in the current directory.

Run the `view` command:

```bash
$ xo view
loading interval data
loading peak data
Launching server at http://localhost:5006
```

That will read the two data files and start the GUI, which uses a web server connected to port 5606.  The server should start your default web browser and show a view of the first chromosome in the data set.

If you browser doesn't start automatically, just start the browser, open a new window, and enter the URL printed in the terminal window (`http://localhost:5006` in the example above).

#### Exiting the GUI

Close the web browser window, and type `^C` in the terminal window where you typed the `xo view` command.

