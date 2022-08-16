import sys
import glob
import re
import yt
import h5py
from unyt import unyt_array
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt


def get_title(profile):

    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"

    dump = re.search(r'^(DD[0-9]{4})_enclosed_mass.h5', profile).group(1)
    fname = '/'.join([dump, dump]).encode()
    sim = yt.load(sim_path)
    time = sim.data['time'][np.where(sim.data['filename'] == fname)][0].in_units('Myr')
    z = sim.data['redshift'][np.where(sim.data['filename'] == fname)][0].value
    title = f"{dump}, z={z:.2f}, t = {time:.2f}"
    return (dump, title)


filenames = glob.glob("DD*_enclosed_mass.h5")
for name in filenames:
    dump, title = get_title(name)
    prof = h5py.File(name, 'r')

    baryonic_mass = unyt_array(list(prof['baryonic_mass/array_data']), 'Msun')
    dm_mass = unyt_array(list(prof['dm_mass/array_data']), 'Msun')
    total_mass = unyt_array(list(prof['total_mass/array_data']), 'Msun')
    radii = unyt_array(list(prof['radius/array_data']), 'pc')

    plt.figure()
    plt.title(title)
    plt.plot(radii, baryonic_mass, color='red', label='baryonic')
    plt.plot(radii, dm_mass, color='blue', label='dm')
    plt.plot(radii, total_mass, color='black', label='total')
    plt.xlabel("Radius (pc)")
    plt.ylabel(r'$Mass (M_{\odot})$')
    if (10 * (radii[1] - radii[0]) < (radii[-1] - radii[-2])):
        plt.xscale('log')
        plt.xlim(1e-2, 5e2)
    else :
        plt.xlim(0.0, None)
    plt.yscale('log')
    plt.ylim(1e-1, 1e7)
    plt.legend(loc='lower right')
    plt.savefig(dump + "_enclosed_mass.png")
    plt.close()
