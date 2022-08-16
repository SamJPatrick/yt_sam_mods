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

    metal3_mass = unyt_array(list(prof['metal3_mass/array_data']), 'Msun/yr')
    metallicity3 = unyt_array(list(prof['metallicity3/array_data']), 'Msun/yr')
    radii = unyt_array(list(prof['radius/array_data']), 'Msun/yr')

    plt.figure()
    plt.title(title)
    _, ax1 = plt.subplots()
    ax1.plot(radii, metal3_mass, color='red')
    if (10 * (radii[1] - radii[0]) < (radii[-1] - radii[-2])):
        ax1.set_xscale('log')
        ax1.set_xlim(1e-2, 5e2)
    else :
        ax1.set_xlim(0.0, None)
    ax1.set_xlabel("Radius (pc)")
    ax1.set_ylabel(r'$M_{metals} (M_{\odot})$', color='red')
    ax1.set_ylim(0.0, 45.0)
    ax2 = ax1.twinx()
    ax2.plot(radii, metallicity3, color='blue')
    ax2.set_ylabel(r'$Z (Z_{\odot})$', color='blue')
    ax2.set_ylim(1e-7, 1e-2)
    ax2.set_yscale('log')
    plt.savefig(dump + "_enclosed_metals.png")
    plt.close()
