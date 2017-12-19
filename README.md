Isothermal model for FCI simulations
====================================

[![License](https://img.shields.io/badge/license-GPL-blue.svg)](https://img.shields.io/badge/license-GPL-blue.svg)

Evolves plasma density, vorticity, parallel current and momentum.

This code was written using [BOUT++ version 4.1](https://github.com/boutproject/BOUT-dev/releases).

    Copyright 2017 Ben Dudson, University of York. Email: benjamin.dudson@york.ac.uk

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Tokamaks
--------

To generate a grid, use [Zoidberg](http://bout-dev.readthedocs.io/en/latest/user_docs/zoidberg.html)
with a 'g' EQDSK format file.

```python
import numpy as np
import zoidberg

field = zoidberg.field.GEQDSK("g014220.00200") # Read magnetic field

grid = zoidberg.grid.rectangular_grid(100, 10, 100,
       1.5-0.1, # Range in R (max - min)
       2*np.pi, # Toroidal angle
       3., # Range in Z
       xcentre=(1.5+0.1)/2, # Middle of grid in R
       yperiodic=True) # Periodic in toroidal angle

# Create the forward and backward maps
maps = zoidberg.make_maps(grid, field)

# Save to file
zoidberg.write_maps(grid, field, maps, gridfile="tokamak.fci.nc")

# Plot grid points and the points they map to in the forward direction
zoidberg.plot.plot_forward_map(grid, maps)
```

Then run the test case

```bash
$ make
$ ./isothermal-4field -d tokamak
```

