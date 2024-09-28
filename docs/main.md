# Main Program

The top level application (`xo.py`) uses the `argparse` library to parse command line
arguments and then call the function that implemements the command specified
on the command line.

### `init_cli`

The default file names for the SNPs file, intervals file, and peaks file are
defined in the `init_cli` function.
First see if there is an environment variable for a file name
(there will be when the program is run in a Docker container), otherwise use
the predefined default.

The argument parser sets up a "subparser" for each command (peaks, gui, or vis).
The subparsers define the options valid for that command and specify which
function to call if the command is found on the command line.  These functions
are all imported from other modules.

**Source**

::: src.xo.xo.init_cli
    options:
      show_root_toc_entry: false

### `main`

::: src.xo.xo.main
    options:
      show_root_toc_entry: false
