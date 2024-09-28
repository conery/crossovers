## Widget Classes

The GUI displayed in a user's browser shows dozens of graphical elements.
These "widgets" include sliders to specify block sizes and buttons that select
different chromosomes.

Our GUI uses a common technique that makes it easier to break up code and work on different widgets separately.
We define our own widget classes, using inheritance so our objects are a special case of an existing type of widget.

A good example is the block size widget.  We want to display a slider with
a label above it, and we want the slider to have values that range from 0 to 100.

The Python syntax for defining a new class that is derived from an existing class uses a `class` statement.
This is the statement that defines our BlockSizeFilterWidget class:
```
class BlockSizeFilterWidget(pn.widgets.IntRangeSlider):
    ...
```
`IntRangeSlider` is an existing class, defined in the Panel library.
The statement above means our new BlockSizeFilterWidget objects will be special types of IntRangeSliders.

The rest of the class defines how our slider should be initialized and what to do
when the slider is moved.


