import matplotlib
import matplotlib.pyplot as plt

from boutdata import collect
from boututils.datafile import DataFile

gridfile = "tokamak-200x10x200.fci.nc"
datapath = "filament-diffusion"
time_index = -1
psi_boundary = 0.02345850877  # Poloidal flux at the separatrix
variable = "n"
log = True
toroidal_extent = 360.  # Degrees
yoffset = 5 # Set this index to be zero degrees

#################

# Get the coordinates and poloidal flux from the input grid
with DataFile(gridfile) as grid:
    psi = grid["psi"][:,0,:]
    R = grid["R"][:,0,:]
    Z = grid["Z"][:,0,:]

data = collect(variable, path=datapath, tind=time_index)

ny = data.shape[2]

t_array = collect("t_array", path=datapath)
time = t_array[time_index]

# Color map
cmap = plt.cm.get_cmap("YlOrBr")

for y in range(ny):
    angle = toroidal_extent * float(y - yoffset)/ny
    
    savefile = variable+"_t%05d_y%02d.pdf" % (int(time),y)
    if log:
        savefile = "log_"+savefile
    
    fig, ax = plt.subplots()
    ax.contour(R, Z, psi, 50, colors='0.75', zorder=0)
    
    ax.contour(R, Z, psi, levels=[psi_boundary], colors='k', zorder=1)

    if log:
        cax = ax.contourf(R, Z, data[0, :, y, :], 50, cmap=cmap, zorder=2, alpha=0.5,
                          locator=matplotlib.ticker.LogLocator())
    else:
        cax = ax.contourf(R, Z, data[0, :, y, :], 50, cmap=cmap, zorder=2, alpha=0.5)
    
    for c in cax.collections:  # Remove white lines between levels
        c.set_edgecolor("face")
    cbar = fig.colorbar(cax)

    ax.set_title(variable + r' at $\phi = %.1f^o$, time = $%.1f/\Omega_{ci}$' % (angle, time))
    ax.set_aspect("equal")

    plt.savefig(savefile)
    #plt.show()
    plt.close()
    
    

