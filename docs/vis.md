# Generate Plots from the Command Line

The `xo vis` command is an alternate method for generating the histograms that summarize blocks across all the chromosomes.

## Command Line Options

The command that runs the program has options to specify settings for each of the filters.  To see the complete list of options use the `--help` option:

```bash 
$ xo vis --help
usage: xo vis [-h] [--peaks F] [--chromosomes P] [--size N N] 
       [--length N N] [--coverage N] [--match] [--save F] P

positional arguments:
  P                type of plot to make ['count', 'length', 'location']

options:
  -h, --help       show this help message and exit
  --peaks F        blocks saved by peaks.py
  --chromosomes P  chromosome name pattern
  --size N N       block size range (#SNPs)
  --length N N     block length range (bp)
  --coverage N     minimum coverage
  --match          require genome match
  --save F         write summary dataframe to this file
```

According to this output, each time we run the program we have to specify a value for `P`, which stands for "plot type."  The possible values are `count`, `length`, and `location`, which correspond to the three kinds of histograms.

Five of the options correspond to the filter widgets and chromosome name widget displayed in the GUI:

- Use `--chromosomes` to specify a chromosome name pattern. If you plan to use this option read the section below on "name patterns in shell commands."  
- Use `--size` and `--length` to specify block size and length limits.  Note that both of these options require two integer arguments, corresponding to the settings of the double-ended slider in the GUI. 
- The `--coverage` option takes a single integer argument, which corresponds to the value of the coverage slider in the GUI.
- Including `--match` on the command line is the same as clicking the toggle button in the match widget.

The remaining two options are for specifying file names.

- The `--peaks` option can be used to specify an alternative to the default `peaks.csv` file with the output from the `peaks` command.
- Use `--save` to have the application write the size and location of each block to a CSV file.

### Default Filter Settings

If a filter option is not specified on the command line the default value will be used.  Each default is the same as the values shown in the GUI when it first starts.

| Option       | Description                             | Default           |
| ------------ | --------------------------------------- | ----------------- |
| `chromosome` | regular expression for chromosome names | `BSP.*`           |
| `size`       | block size (number of SNPs)             | 0 to 10           |
| `length`     | block length (base pairs)               | 0 1o 10000        |
| `coverage`   | minimum SNP coverage                    | 0                 |
| `match`      | require genome match                    | False (unchecked) |

The default is `BSP.*`, meaning "all chromosomes."

If these options are not specified the defaults are the same as the initial settings shown in the GUI:  0 to 10 for block size, and 0 to 10,000 for block length.

### Name Patterns in Shell Commands

The argument for the `--chromosome` option is a regular expression, and these expressions often contain asterisks.  For example, the pattern `BSP-.*-1` means "all names that start with `BSP-` and end with `-1` with any characters in between", in other words, "chromosome 1 from all worms."

We need to be careful when we use this pattern in a shell command.  When a shell sees an asterisk, it thinks you are typing the name of a file, and it tries to look for all files with names that match the pattern.

To prevent this behavior, we have to **enclose the pattern in single quotes**.  This is how we would create a histogram of block sizes using only chromosome 1 from each worm:

```bash
$ xo vis count --chromosomes 'BSP-.*-1'
```

> **Note:**  File name expansion -- called "globbing" in the Unix world -- is done by the shell.  The shell tries to look for files with matching names before it even starts the `xo` command.  If you leave out the quotes you're likely to see an error message like this:
>
> ```
> zsh: no matches found: BSP-.*-1
> ```
>
> Notice how this error is coming from the shell (`zsh` in this case, the default shell on macOS) and not from the `xo` program, which is never even started.

### Examples

For other examples of chromosome name patterns see the Chromosome Names section of the GUI documentation.

