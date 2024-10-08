# Crossover Explorer

The **crossover explorer** is an application that filters and displays SNPs, with the goal of looking for evidence of crossover and non-crossover events in _C. elegans_ genomes.

The top-level application runs three different scripts that implement steps in the overall workflow:

- a "peak finder" script scans a set of SNPs generated by TIGER to look for chromosomal regions that might be the result of crossover or non-crossover events
- another script presents a graphical user interface that lets users browse the data generated by the peak finder
- a visualization script creates histograms that summarize the size and location of the regions

## Data Files

The application uses two sets of data:

* `BSP_TIGER.marker_dataframe.pickle.gzip` is a "pickled" and compressed Pandas data frame with the output from TIGER.  Each row in the frame describes a SNP, with columns for the
chromosome name, location, and the predicted parent genome (N2 or CB4856), and more.

* `BSP_TIGER.intervals_dataframe.pickle.gzip`, a summary of the SNP data, where each row defines a chromosome segment and its predicted parent.

We suggest creating a new folder to use for a project directory.  Move (or link) the two data files to this folder, then `cd` to the directory and run the crossover explorer scripts in that directory.

> **Note for Docker users**:  if you are running the GUI in a Docker container you will bind mount this same directory when starting the container.

## Shell Commands

To run one of the scripts type `xo` followed by one of the script names (`peaks`, `gui`, or `vis`).   Each of the scripts has its own help message, which you can see by adding `--help` to the command.

#### Examples

This command prints the help message for the top level `xo` application:

```bash
$ xo --help
usage: xo [-h] {peaks,gui,vis,post} ...

options:
  -h, --help            show this help message and exit

subcommands:
  operation to perform

  {peaks,gui,vis,post}
    peaks               find peaks in the SNP data
    gui                 explore blocks of SNPs
    vis                 visualizations based on filtered blocks
```

To see the help message for one of the scripts type `xo`, the script name, and then `--help`.  This shell command prints the help message for the `peaks` script:

```bash
$ xo peaks --help
usage: xo peaks [-h] [--snps F] [--output F] [--max_snps N]
...
```

## Abbreviating Options

The `xo` scripts, like most modern Unix command line applications, allow users to shorten option names, so that it is only necessary to type enough characters to distinguish one option from another.

The full script name (`peaks`, `gui`, or `vis`) must be entered completely, but after that any option names can be abbreviated.

#### Example

These two commands are equivalent:

```bash
$ xo peaks --output short_blocks.csv --max_snps 10
$ xo peaks --out short_blocks.csv --m 10
```

