import sys
import h5py
import glob
import os
from unyt import Msun, pc
from unyt import unyt_array
import numpy as np

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import *


INDIR = "Profiles/Pressure_profiles"
OUTDIR = "Profiles/Pressure_graphs"


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass

fields = ['hydrostatic', 'thermal', 'ram', 'turbulent']
filenames = glob.glob(os.path.join(INDIR, "DD*_pressures.h5"))
for name in filenames:
    filename = name.split('/')[-1]
    prof = h5py.File(name, 'r')
    plt.figure()
    plt.title(get_title(filename, star_type))
    radii = unyt_array(list(prof['radius/array_data']), 'pc')
    for field in fields:
        p_profile = list(prof['/'.join([field, 'array_data'])])
        plt.plot(radii.in_units('pc'), p_profile, label= field)
    plt.xlabel("Radius (pc)")
    plt.ylabel("Pressure (Pa)")
    plt.xscale('log')
    plt.xlim(1e-2, 5e2)
    plt.yscale('log')
    plt.ylim(1e-19, 1e-11)
    plt.legend(loc= "upper right")
    plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(filename)}_pressure_profiles.png"))
    plt.close()
