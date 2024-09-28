# Chromosome Display

The function that displays a chromosome is a method named `display_chromosome`,
defined in the PeakViewerApp class.

The method is called whenver the user clicks the left or right button next to
the chromosome name or edits the chromosome name.

To break the code into manageable pieces there are two "helper methods":

* `_make_patches` builds a list of rectangles for the "ribbon" display
* `_make_grid` builds the column of blocks, interleaving graphical representations
of the SNPs with a text view of the dataframe

### `display_chromosome`

::: src.xo.gui.PeakViewerApp.display_chromosome
    options:
      show_root_toc_entry: false

### `_make_patches`

::: src.xo.gui.PeakViewerApp._make_patches
    options:
      show_root_toc_entry: false

### `_make_grid`

::: src.xo.gui.PeakViewerApp._make_grid
    options:
      show_root_toc_entry: false
