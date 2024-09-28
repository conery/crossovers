# PeakViewerApp

The PeakViewerApp class is the top level of the GUI.
When the class is instantiated, it creates all the components (buttons, sliders, labels, _etc_) and lays them out in the browser window.

**Layout**

The app is based on a Panel template named BootstrapTemplate.
Applications that use this template have a title at the top of the window,
a sidebar on the left side, and a main content area.

In our application, we put the filter widgets in the sidebar.
The main content is a "tab" widget with two tabs, one to show a single chromosome and the
other to show summary plots.
Each tab also needs its own layout:

* the chromosome tab has controls to navigate to different chromosomes
and a column that has graphics and tables for the blocks in the current chromosome
* the summary tab has a text widget for entering a chromosome name pattern and a set of buttons to make plots

The initialization code has a series of assignment statements that create objects and save them as instance variables, _e.g._
```
    self.back_button = pn.widgets.Button(name='◀︎', stylesheets=[button_style_sheet])
    self.forward_button = pn.widgets.Button(name='▶︎', stylesheets=[button_style_sheet])
```

After widgets are created they need to be placed in a layout.  This
statement lays out the structure of the chromosome tab:
```
    chr_tab = pn.Column(
        pn.pane.HTML('<h3>Chromosome</h3>'),
        pn.Row(self.back_button, self.chromosome_id, self.forward_button),
        ...
    )
```

Eventually all the layouts are saved in either the sidebar or the main section.

**Callbacks**

Another important role for the initialization function is to connect widgets to code that will be executed when the widget is activated, _i.e._ when the button on a scrollbar is
moved, or a button is clicked.

This is accomplished by telling the application to "watch" the value of a widget.
When the value changes we want to call a function to respond to the change.
This example shows how we tell the application to call a method named
`chromosome_edited_cb` whenver the text in the chromosome name widget is updated:
```
    self.chromosome_id.param.watch(self.chromosome_edited_cb, ['value'])  
```

**Initialization Method**

::: src.xo.gui.PeakViewerApp
    options:
      show_root_toc_entry: false
      docstring_options:
        ignore_init_summary: true
      merge_init_into_class: true
      heading_level: 3
      filters: ""
      members:
        - attach_callbacks
        - load_data
