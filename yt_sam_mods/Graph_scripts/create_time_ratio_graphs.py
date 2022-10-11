import sys
import h5py
import glob
from unyt import Msun, pc
from unyt import unyt_array
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import get_title, get_dump_num, get_time_z




be_info = h5py.File("be_time_mass_radius.h5", 'r')
be_times = unyt_array(list(be_info['time/array_data']), 'Myr')
be_norm_radii = unyt_array(list(be_info['radius_norm/array_data']), 'pc')
be_turb_radii = unyt_array(list(be_info['radius_turb/array_data']), 'pc')
be_norm_masses = unyt_array(list(be_info['mass_norm/array_data']), 'pc')
be_turb_masses = unyt_array(list(be_info['mass_turb/array_data']), 'pc')

filenames = glob.glob("DD*_bonnor_ebert.h5")
fields = ['tcool_ratio_norm', 'tcool_ratio_turb', 'tff_ratio_norm', 'tff_ratio_turb']
linestyles = ['-', '-.', '-', '-.']
colors = ['blue', 'blue', 'orange', 'orange']
assert len(colors) == len(linestyles) == len(colors)
for fn in filenames:
    
    i, time_diff = min(enumerate(np.abs(be_times - get_time_z(fn)[0])), key= lambda x: x[1])
    be_norm_radius = be_norm_radii[i]
    be_turb_radius = be_turb_radii[i]
    title = get_title(fn)
    profile = h5py.File(fn, 'r')

    plt.figure()
    plt.title(title)
    radii = unyt_array(list(profile['radius/array_data']), 'pc')
    for j, field in enumerate(fields):
        time_prof = list(profile['/'.join([field, 'array_data'])])
        plt.plot(radii.in_units('pc'), time_prof, linestyle= linestyles[j], color= colors[j], label= field)
    plt.xlabel("Radius (pc)")
    plt.ylabel("timescale")
    plt.xscale('log')
    plt.xlim(1e-2, 5e2)
    plt.yscale('log')
    plt.ylim(1e-1, 1e6)
    plt.axvline(x= be_norm_radius)
    plt.axvline(x= be_turb_radius, linestyle='-.')
    plt.axhline(y= 1.0)
    plt.legend(loc= "upper left")
    plt.savefig(f"DD{get_dump_num(fn)}_be_timescales.png")
    plt.close()
