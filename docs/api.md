# Command Line API

The GUI and the filtering steps are all implemented as part of a command line application named `xo`.
To run one of the steps simply type `xo` followed by the name of the operation you want to perform (`peaks`, `gui`, `filter`, or `post`).
Each of the scripts has its own help message, which you can see by adding `--help` to the command.

## Help

This command prints the help message for the top level `xo` application:

```bash
$ xo --help
usage: xo [-h] [--log X] {gui,peaks,filter,post} ...

options:
  -h, --help            show this help message and exit
  --log X

subcommands:
  operation to perform

  {gui,peaks,filter,post}
    gui                 explore blocks and NCOs
    peaks               find blocks around peaks in the SNP data
    filter              apply filters to blocks
    post                postprocessing of filtered blocks
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

### Example

These two commands are equivalent:

```bash
$ xo peaks --output short_blocks.csv --max_snps 10
$ xo peaks --out short_blocks.csv --m 10
```
