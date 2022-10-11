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

from yt.extentions.sam_mods.graph_funcs import get_time_z




states = ['frag', 'support_turb', 'collapse', 'support_pressure']
filenames = glob.glob("DD*_state_info.h5")
filenames.sort(key= lambda x: get_time_z(x)[0])
times = unyt_array(np.zeros(len(filenames)), 'Myr')

gas_masses = [unyt_array(np.zeros(len(filenames)), 'Msun') for state in states]
gas_frac = [unyt_array(np.zeros(len(filenames)), '') for state in states]
mean_radii = [unyt_array(np.zeros(len(filenames)), 'pc') for state in states]
radius_frac = [unyt_array(np.zeros(len(filenames)), '') for state in states]
for i, name in enumerate(filenames):
    times[i] = get_time_z(name)[0]
    state_info = h5py.File(name, 'r')
    for j, state in enumerate(states):
        gas_masses[j][i] = state_info['gas_mass'][j]
        gas_frac[j][i] = state_info['gas_fraction'][j]
        mean_radii[j][i] = state_info['mean_radius'][j]
        radius_frac[j][i] = state_info['radius_fraction'][j]


plt.figure()
plt.title("Mass of different gas phases")
for i, state in enumerate(gas_masses):
    plt.plot(times.in_units('Myr'), state.in_units('Msun'), label= states[i])
plt.xlabel("Time (Myr)")
plt.ylabel(r"$Gas ~mass ~(M_{\odot})$")
plt.yscale('log')
plt.ylim(1e3, 2e5)
plt.legend(loc= "upper left")
plt.savefig("gas_masses_vs_times.png")
plt.close()

plt.figure()
plt.title("Mass weighted mean radius of different gas phases")
for i, state in enumerate(mean_radii):
    plt.plot(times.in_units('Myr'), state.in_units('pc'), label= states[i])
plt.xlabel("Time (Myr)")
plt.ylabel(r"Mean radius (pc)")
plt.ylim(0.0, 260.0)
plt.legend(loc= "upper left")
plt.savefig("mean_radii_vs_time.png")
plt.close()

plt.figure()
plt.title("Mass (as fraction of total) of different gas phases")
for i, state in enumerate(gas_frac):
    plt.plot(times.in_units('Myr'), state.in_units(''), label= states[i])
plt.xlabel("Time (Myr)")
plt.ylabel(r"Gas fraction")
plt.ylim(0.0, 1.0)
plt.legend(loc= "upper left")
plt.savefig("gas_fracion_vs_time.png")
plt.close()

plt.figure()
plt.title("Mass weighted mean radius (as fraction of virial) of different gas phases")
for i, state in enumerate(radius_frac):
    plt.plot(times.in_units('Myr'), state.in_units(''), label= states[i])
plt.xlabel("Time (Myr)")
plt.ylabel(r"Radius fraction")
plt.ylim(0.0, 1.0)
plt.legend(loc= "upper left")
plt.savefig("radius_fraction_vs_time.png")
plt.close()
