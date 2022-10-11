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



fields = ['hydrostatic', 'thermal', 'ram', 'turbulent']
filenames = glob.glob("DD*_pressures.h5")
for name in filenames:
    title = get_title(name)
    prof = h5py.File(name, 'r')

    plt.figure()
    plt.title(title)
    radii = unyt_array(list(prof['radius/array_data']), 'pc')
    for field in fields:
        p_profile = list(prof['/'.join([field, 'array_data'])])
        plt.plot(radii.in_units('pc'), p_profile, label= field)
    plt.xlabel("Radius (pc)")
    plt.ylabel("Pressure (Pa)")
    plt.xscale('log')
    plt.xlim(1e-2, 5e3)
    plt.yscale('log')
    plt.ylim(1e-19, 1e-11)
    plt.legend(loc= "upper right")
    plt.savefig(f"DD{get_dump_num(name)}_pressure_profiles.png")
