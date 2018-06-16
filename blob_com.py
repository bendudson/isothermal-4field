#
# Centre-of-mass tracking of plasma blob.
# Calculates R, Z, psi as a function of time on a given toroidal slice (y index)

datasets = [("filament", "tokamak-ny10.fci.nc", "100eV")]
yindex = 5

from boutdata import collect
from boututils.datafile import DataFile
import numpy as np
import matplotlib.pyplot as plt

for datapath, gridfile, label in datasets:

    # Read values from the grid file
    with DataFile(gridfile) as g:
        R = g["R"][:,yindex,:]
        Z = g["Z"][:,yindex,:]
        #psi = g["psi"][:,yindex,:]
    
    
    # Read the density
    n = collect("n", yind=yindex, path=datapath)[:,:,0,:]
    
    time = collect("t_array", path=datapath)
    wci = collect("Omega_ci", path=datapath)
    time *= 1e6 / wci # microseconds
    
    nt = n.shape[0]
    
    # Arrays to contain the R, Z and psi coordinates of the blob
    Rcom = np.ndarray(nt)
    Zcom = np.ndarray(nt)
    #psicom = np.ndarray(nt)
    
    for t in range(nt):
        n_cell = n[t,:,:] * R # Note: Volume of cell proportional to R
        ntot = np.sum(n_cell) 
        Rcom[t] = np.sum(n_cell * R) / ntot
        Zcom[t] = np.sum(n_cell * Z) / ntot
        #psicom[t] = np.sum(n_cell * psi) / ntot

    Rcom *= 1e3 # mm
        
    plt.plot(time - time[0], Rcom - Rcom[0], label=label)

plt.xlabel(r"Time $\mu s$")
plt.ylabel("Major radius [mm]")
plt.legend()

plt.savefig("blob_com.png")
plt.savefig("blob_com.pdf")

plt.show()


