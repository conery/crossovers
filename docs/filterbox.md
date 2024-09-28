# Filter Widgets

Five new classes are defined the beginning of `gui.py`.
Four of them are used to define the widgets that control filter parameters:

* BlockSizeFilterWidget is an integer range slider to specify minimum and maximum block sizes
* BlockLengthFilterWidget is another integer range slider, for block lengths
* CoverageFilterWidget is a single-value integer slider
* SupportFilterWidget is a checkbox that is either on or off

The initialization methods for each widget defines the widget's name (which will be
displayed above the widget when it is shown in the GUI) and sets its initial values.

Initialization methods are also passed a reference to a SNPFilter
object that has the code that will do the actual filtering before a chromosome is displayed.
When a filter widget is updated (_e.g._ when the user moves a slider) a "callback function"
is activated.  The function will record the new setting in the SNPFilter so it is used
the next time the filter is called.

The fifth class, called FilterBox, is a type of widget called a "layout".
It provides a way to collect all the filter widgets in a single location and provides a common interface to all the filters.

## FilterBox

::: src.xo.gui.FilterBox
    options:
      show_root_toc_entry: false
      docstring_options:
        ignore_init_summary: true
      merge_init_into_class: true
      heading_level: 3
      filters: ""
      members_order: source

## BlockSizeFilterWidget

::: src.xo.gui.BlockSizeFilterWidget
    options:
      show_root_toc_entry: false
      docstring_options:
        ignore_init_summary: true
      merge_init_into_class: true
      heading_level: 3
      filters: ""
      members_order: source

## BlockLengthFilterWidget

::: src.xo.gui.BlockLengthFilterWidget
    options:
      show_root_toc_entry: false
      docstring_options:
        ignore_init_summary: true
      merge_init_into_class: true
      heading_level: 3
      filters: ""
      members_order: source

## CoverageFilterWidget

::: src.xo.gui.CoverageFilterWidget
    options:
      show_root_toc_entry: false
      docstring_options:
        ignore_init_summary: true
      merge_init_into_class: true
      heading_level: 3
      filters: ""
      members_order: source

## SupportFilterWidget

::: src.xo.gui.SupportFilterWidget
    options:
      show_root_toc_entry: false
      docstring_options:
        ignore_init_summary: true
      merge_init_into_class: true
      heading_level: 3
      filters: ""
      members_order: source
