import sys
import h5py
import glob
from unyt import Msun, pc
from unyt import unyt_array
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import get_title, get_dump_num



filenames = glob.glob("DD*_pressures.h5")
for name in filenames:
    title = get_title(name)
    prof = h5py.File(name, 'r')

    plt.figure()
    plt.title(title)
    
    radii = np.array(prof['radius/array_data'])
    ratio_1 = np.abs(np.array(prof['turbulent/array_data']) - np.array(prof['ram/array_data'])) / \
        np.abs(np.array(prof['turbulent/array_data']) - np.array(prof['thermal/array_data']))
    ratio_2 = np.array(prof['thermal/array_data']) / np.array(prof['hydrostatic/array_data'])
        
    plt.plot(radii, ratio_1, label= r"$\Delta_{ram} / \Delta_{thermal}$")
    plt.plot(radii, ratio_2, label= r"$p_{thermal} / p_{hydrostatic}$")
    
    plt.axhline(y= 1.0)
    plt.xlabel("Radius (pc)")
    plt.xscale('log')
    plt.xlim(1e-2, 5e3)
    plt.yscale('log')
    plt.ylim(1e-4, 1e3)
    plt.legend(loc= "upper right")
    plt.savefig(f"DD{get_dump_num(name)}_pressure_ratio_profiles.png")
