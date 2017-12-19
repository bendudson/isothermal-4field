Isothermal model for FCI simulations
====================================

Evolves plasma density, vorticity, parallel current and momentum.

This code was written for [BOUT++ version 4.1](https://github.com/boutproject/BOUT-dev/releases).

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

