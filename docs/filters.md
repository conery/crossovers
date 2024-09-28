# Filters

The `filter` module defines a class named SNPFilter.
An app creates an instance of this class to act as an interface to the full set of SNPs.
The general workflow is:

* call the method that loads SNP data from a CSV file
* assign values for filtering criteria (block size, _etc_)
* call the `apply` method to apply the filters

::: src.xo.filters
    options:
      show_root_toc_entry: false
      docstring_options:
        ignore_init_summary: true
      merge_init_into_class: true
      heading_level: 3
      filters: ""
      members_order: source
