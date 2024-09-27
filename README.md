# Crossover Explorer

The **crossover explorer** is a collection of Python applications that filter and display SNP data, with the goal of looking for evidence of crossover and non-crossover events in _C. elegans_ genomes.

There are three applications in the current version:

- a "peak finder" that searches chromosomes for blocks of potentially interesting SNPs
- a visualizer generates histograms that summarize sizes and locations of blocks
- a GUI allows users to explore blocks, filtering by block size and other attributes

The easiest way to run the GUI is through a Docker container.  Instructions for pulling the image and starting a container can be found in `README.Docker`. 

The GUI and the other two applications can be installed via Pip and run from the command line.  Install the apps using the instructions below, then read more about the commands and their options at [Crossover Explorer](https://conery.github.io/crossovers/).

## Installation

All of the applications are run via a top-level script named `xo`.  This command will install the script from GitHub:

```
pip3 install git+https://github.com/conery/crossovers.git
```

To make sure it's installed:

```bash
$ xo --help
usage: xo [-h] {peaks,gui,vis} ...
...
```

To run one of the applications type `xo` followed by the application names and any options that are valie for that application.  

All applications have a `--help` option.  This command prints the help message for the `gui` option:

```bash
$ xo gui --help
```

## Documentation

Detailed instructions for running the applications can be found at [Crossover Explorer](https://conery.github.io/crossovers/).  That site also has developer documentation for anyone who wants to modify or extend the applications.



