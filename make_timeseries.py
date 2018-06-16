#
# Make a series of images to turn into an animation

gridfile = "tokamak-ny10.fci.nc"
datapath = "filament"
yindex = 5
angle = 0.  # Toroidal angle (for title labelling)
variable = "n"
log = False

psi_boundary = 0.02345850877  # Poloidal flux at the separatrix


# Range of R and Z to plot
Rrange = [0, 2]
Zrange = [-2,2]

from boutdata import collect
from boututils.datafile import DataFile
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Read values from the grid file
with DataFile(gridfile) as g:
    R = g["R"][:,yindex,:]
    Z = g["Z"][:,yindex,:]
    try:
        psi = g["psi"][:,yindex,:]
    except KeyError:
        psi = None
    

# Turn the range of R and Z into index

R1d = R[:,0]
Z1d = Z[0,:]

xind = [np.argmin(abs(R1d - r)) for r in Rrange]
zind = [np.argmin(abs(Z1d - z)) for z in Zrange]

R = R[xind[0]:(xind[1]+1),zind[0]:(zind[1]+1)]
Z = Z[xind[0]:(xind[1]+1),zind[0]:(zind[1]+1)]
if psi is not None:
    psi = psi[xind[0]:xind[1],zind[0]:zind[1]]

# Read the data

data = collect(variable, path=datapath, yind=yindex, xind=xind, zind=zind)[:,:,0,:]
time = collect("t_array", path=datapath)

nt = len(time) # Number of time points

# Get the levels to be used for contourf
data_min = np.amin(data)
data_max = np.amax(data)

def ceil_limit(value):
    """
    Rounds up values to 1 sig fig

    0.0032 -> 0.004
    12 -> 20.0
    -0.14 -> -0.1
    """
    e = np.floor(np.log10(np.abs(value)))
    digit = np.ceil(value * 10**(-e))
    return digit * 10**e

def floor_limit(value):
    """
    Rounds down values to 1 sig fig

    0.0032 -> 0.003
    12 -> 10.0
    -0.14 -> -0.2
    """
    e = np.floor(np.log10(np.abs(value)))
    digit = np.floor(value * 10**(-e))
    return digit * 10**e

data_max = ceil_limit(data_max)
data_min = floor_limit(data_min)

print(data_max, data_min)

nlevels = 64
if log:
    log_min = np.log(data_min)
    log_max = np.log(data_max)
    levels = np.exp( np.arange(nlevels)/(nlevels-1) * (log_max - log_min) + log_min )
else:
    levels = np.arange(nlevels)/(nlevels-1) * (data_max - data_min) + data_min

# Color map
cmap = plt.cm.get_cmap("YlOrBr")

fileprefix = variable+"_y%02d" % (yindex,)
if log:
    fileprefix = "log_" + fileprefix

for t in range(nt):
    savefile = fileprefix + "_t%04d.png" % (t,)

    fig, ax = plt.subplots()
    if psi is not None:
        ax.contour(R, Z, psi, 50, colors='0.75', zorder=0)
        ax.contour(R, Z, psi, levels=[psi_boundary], colors='k', zorder=1)
    
    cax = ax.contourf(R, Z, data[t, :, :], levels=levels, cmap=cmap, zorder=2, alpha=0.5)

    for c in cax.collections:  # Remove white lines between levels
        c.set_edgecolor("face")
        
    cbar = fig.colorbar(cax)
    #cbar.solids.set_rasterized(True)
    #cbar.solids.set_edgecolor("face")

    ax.set_title(variable + r' at $\phi = %.1f^o$, time = $%.1f/\Omega_{ci}$' % (angle, time[t]))
    ax.set_xlabel("Major radius [m]")
    ax.set_ylabel("Height [m]")
    ax.set_aspect("equal")
    
    plt.title(variable + " t = %04d" % (t,))
    plt.savefig(savefile)
    #plt.show()
    plt.close()

print("Now run $ convert -loop 0 -delay 100 " + fileprefix + "_t*.png " + fileprefix + ".gif")

