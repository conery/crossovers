# Filters

All of the code for the crossover explorer GUI is in a single Python file
named `gui.py`.
The main sections of this file are:

* definitions of the "filter widgets" (sliders and checkboxes for setting filtering options)
* a class named `PeakViewerApp` that defines the Panel app and lays out the widgets;
the class also defines methods that display a chromosome and handle events by invoking "callback" functions
* code that creates and launches the viewer

