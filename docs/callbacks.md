# Callback Functions

The functions that are called when a widget is activated all have the same form:
```
def cb(self, e):
    ...
```
The first argument is always `self` because the callbacks in this application are all
methods defined in the PeakViewerApp class.
That allows them to access the data and other parts of the GUI if needed.

The second argument is an **event** object.
It has several attributes, but the one we use most often is a reference to the widget
that caused the event.

Here is an example. 
When one of the buttons next to the chromosome name in the main display is clicked we want to update the display to show the next or previous chromsome in the data set.
The code that sets up the application attaches this callback to both buttons:
```
def change_chromosome_cb(self, e):
    delta = 1 if e.obj is self.forward_button else -1
    self.chr_index = (self.chr_index + delta) % len(self.clist)
    ...
```
The goal is to update the value of the variable that holds the current chromosome ID.
If the user clicked the forward button we want to add 1, otherwise we want to subtract 1.
The first statement looks at `e.obj`, which is a reference to the widget that triggered this event, and assigns `delta` to be either 1 or -1, depending on which button was clicked.
The second line updates the chromosome ID.

### `filter_cb`

::: src.xo.gui.PeakViewerApp.filter_cb
    options:
      show_root_toc_entry: false

### `change_chromosome_cb`

::: src.xo.gui.PeakViewerApp.change_chromosome_cb
    options:
      show_root_toc_entry: false

### `chromosome_edited_cb`

::: src.xo.gui.PeakViewerApp.chromosome_edited_cb
    options:
      show_root_toc_entry: false

### `toggle_text_cb`

::: src.xo.gui.PeakViewerApp.toggle_text_cb
    options:
      show_root_toc_entry: false

### `summary_plot_cb`

::: src.xo.gui.PeakViewerApp.summary_plot_cb
    options:
      show_root_toc_entry: false

### `download_cb`

::: src.xo.gui.PeakViewerApp.download_cb
    options:
      show_root_toc_entry: false

