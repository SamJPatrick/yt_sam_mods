import yt
import os
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
#from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
import numpy as np


CLUMP_DIR = "Clumps"
INPUT_FILE = "DD0295_clump_info_pisn.h5"
OUTPUT_FILE = "clump_ratios_pisn_295.png"

def run_ratios(clump, vol_ratios, mass_ratios, densities):

    vol_ratio = clump[('clump', 'volume')] / \
        np.sum([child[('clump', 'volume')] for child in clump.children])
    mass_ratio = clump[('clump', 'cell_mass')] / \
        np.sum([child[('clump', 'cell_mass')] for child in clump.children])
    density = clump[('clump', 'min_number_density')]

    vol_ratios.append(vol_ratio)
    mass_ratios.append(mass_ratio)
    densities.append(density)

    for child in clump.children:
        if child.children != None:
            vol_ratios, mass_ratios, densities = \
                run_ratios(child, vol_ratios, mass_ratios, densities)
    return vol_ratios, mass_ratios, densities



if __name__ == "__main__":
            
    df = yt.load(os.path.join(CLUMP_DIR, INPUT_FILE))
    master_clump = df.tree
    vol_ratios = []
    mass_ratios = []
    densities = []
    vol_ratios, mass_ratios, densities = \
        run_ratios(master_clump, vol_ratios, mass_ratios, densities)

    plt.figure()
    plt.plot(densities, vol_ratios, color= 'red', label='volume')
    plt.plot(densities, mass_ratios, color= 'green', label='gas mass')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e0, 1e8)
    plt.xlim(1e-1, 1e10)
    plt.xlabel(r"$n_{\min} ~(cm^{-3})$")
    plt.ylabel(r"$\phi_{parent}/\Sigma \phi_{children}$")
    plt.legend(loc= 'upper left')
    plt.savefig(OUTPUT_FILE)

    

