# Source Directory

The source code is in a folder named `src`.  Inside that folder is another
folder named `xo`, which is a Python module (since it has an `__init__.py` file).

```bash
src
└── xo
    ├── __init__.py
    ├── filters.py
    ├── gui.py
    ├── peaks.py
    ├── vis.py
    └── xo.py
```

The top level program is in `xo.py`.

The code for the three commands -- `peaks`, `gui`, and `vis` -- are in files with
the same name as the command.

The file named `filters.py` has the definition of a SNPFilter class that is used
by the GUI and the `vis` script.

