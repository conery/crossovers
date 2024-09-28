# `vis.py`

This script provides another way to generate a plot that summarizes blocks of SNPs.

The command line arguments are:

* the type of the histogram to make (`count`, `length`, `location`); these
correspond to the names on the buttons in the summary tab of the GUI

* filtering paramters that also correspond to widgets in the GUI; specify two
integers for size and length ranges, one integer for coverage


::: src.xo.vis
    options:
      show_root_toc_entry: false
      docstring_options:
        ignore_init_summary: true
      merge_init_into_class: true
      heading_level: 3
      filters: ""
      members_order: source
