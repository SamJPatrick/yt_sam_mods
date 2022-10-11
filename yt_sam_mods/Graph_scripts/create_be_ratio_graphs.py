import sys
import h5py
import glob
from unyt import Msun, pc
from unyt import unyt_array
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import get_dump_num, get_title, get_time_z


FILENAME = "be_time_mass_radius.h5"


be_files = glob.glob("DD*_bonnor_ebert.h5")
be_files.sort(key= lambda x: int(get_dump_num(x)))
intcp_radius_norm = np.zeros(len(be_files))
intcp_radius_turb = np.zeros(len(be_files))
intcp_mass_norm = np.zeros(len(be_files))
intcp_mass_turb = np.zeros(len(be_files))
max_radius_norm = np.zeros(len(be_files))
max_radius_turb = np.zeros(len(be_files))
max_mass_norm = np.zeros(len(be_files))
max_mass_turb = np.zeros(len(be_files))
times = np.zeros(len(be_files))

for i, be_file in enumerate(be_files):
    
    be_prof = h5py.File(be_file, 'r')
    be_ratio_norm = list(be_prof['be_ratio_norm/array_data'])
    be_ratio_turb = list(be_prof['be_ratio_turb/array_data'])
    assert len(be_ratio_norm) == len(be_ratio_turb), "Error, ratio profiles of unequal length"
    nbins = len(be_ratio_norm)
    norm_max_index, norm_max_ratio = max([(j, ratio) for j, ratio in enumerate(be_ratio_norm)], key=lambda x: x[1])
    turb_max_index, turb_max_ratio = max([(j, ratio) for j, ratio in enumerate(be_ratio_turb)], key=lambda x: x[1])
    norm_intcp_index, norm_intcp_ratio = max([(j, np.abs(ratio - 1)) for j, ratio in enumerate(be_ratio_norm)], key=lambda x: x[1])
    turb_intcp_index, turb_inctp_ratio = max([(j, np.abs(ratio - 1)) for j, ratio in enumerate(be_ratio_turb)], key=lambda x: x[1])
    if (norm_intcp_index == nbins or norm_intcp_index == 0):
        norm_intcp_index = norm_max_index
    if (turb_intcp_index == nbins or turb_intcp_index == 0):
        turb_intcp_index = turb_max_index
    intcp_radius_norm[i] = be_prof['radius/array_data'][norm_intcp_index]
    intcp_radius_turb[i] = be_prof['radius/array_data'][turb_intcp_index]
    intcp_mass_norm[i] = be_prof['total_mass/array_data'][norm_intcp_index]
    intcp_mass_turb[i] = be_prof['total_mass/array_data'][turb_intcp_index]
    max_radius_norm[i] = be_prof['radius/array_data'][norm_max_index]
    max_radius_turb[i] = be_prof['radius/array_data'][turb_max_index]
    max_mass_norm[i] = be_prof['total_mass/array_data'][norm_max_index]
    max_mass_turb[i] = be_prof['total_mass/array_data'][turb_max_index]
    times[i] = get_time_z(be_file)[0]

    radii = list(be_prof['radius/array_data'])
    plt.figure()
    plt.title(get_title(be_file))
    plt.plot(radii, be_ratio_norm, label='pressure')
    plt.plot(radii, be_ratio_turb, label='turbulence')
    plt.axhline(y=1.0)
    plt.xlabel("Radius (pc)")
    plt.xscale('log')
    plt.xlim(1e-2, 5e3)
    plt.ylabel(r"$M_{tot} / M_{BE}$")
    plt.yscale('log')
    plt.ylim(1e-2, 1e2)
    plt.savefig(f"DD{get_dump_num(be_file)}_be_ratio_radius.png")
    plt.close()

    masses = list(be_prof['total_mass/array_data'])
    plt.figure()
    plt.title(get_title(be_file))
    plt.plot(masses, be_ratio_norm, label='pressure')
    plt.plot(masses, be_ratio_turb, label='turbulence')
    plt.axhline(y=1.0)
    plt.xlabel(r"$Enclosed ~baryon ~mass ~(M_{\odot})$")
    plt.xscale('log')
    plt.xlim(1e-2, 1e7)
    plt.ylabel(r"$M_{tot} / M_{BE}$")
    plt.yscale('log')
    plt.ylim(1e-2, 1e2)
    plt.savefig(f"DD{get_dump_num(be_file)}_be_ratio_mass.png")
    plt.close()
    

intcp_radius_norm = unyt_array(intcp_radius_norm)
intcp_radius_turb = unyt_array(intcp_radius_turb)
intcp_mass_norm = unyt_array(intcp_mass_norm)
intcp_mass_turb = unyt_array(intcp_mass_turb)
max_radius_norm = unyt_array(max_radius_norm)
max_radius_turb = unyt_array(max_radius_turb)
max_mass_norm = unyt_array(max_mass_norm)
max_mass_turb = unyt_array(max_mass_turb)
times = unyt_array(times, 'Myr')

intcp_radius_norm.write_hdf5(FILENAME, group_name='intcp_radius_norm')
intcp_radius_turb.write_hdf5(FILENAME, group_name='intcp_radius_turb')
max_radius_norm.write_hdf5(FILENAME, group_name='max_radius_norm')
max_radius_turb.write_hdf5(FILENAME, group_name='max_radius_turb')
intcp_mass_norm.write_hdf5(FILENAME, group_name='intcp_mass_norm')
intcp_mass_turb.write_hdf5(FILENAME, group_name='intcp_mass_turb')
max_mass_norm.write_hdf5(FILENAME, group_name='max_mass_norm')
max_mass_turb.write_hdf5(FILENAME, group_name='max_mass_turb')
times.write_hdf5(FILENAME, group_name='time')

'''
df = h5py.File(FILENAME, 'r')
intcp_radius_norm = df['intcp_radius_norm/array_data']
intcp_radius_turb = df['intcp_radius_turb/array_data']
max_radius_norm = df['max_radius_norm/array_data']
max_radius_turb = df['max_radius_turb/array_data']
intcp_mass_norm = df['intcp_mass_norm/array_data']
intcp_mass_turb = df['intcp_mass_turb/array_data']
max_mass_norm = df['max_mass_norm/array_data']
max_mass_turb = df['max_mass_turb/array_data']
times = df['time/array_data']
'''

plt.figure()
plt.title("Bonnor-Ebert radius")
plt.plot(times, intcp_radius_norm, label='normal', color='blue', linestyle='-')
plt.plot(times, intcp_radius_turb, label='turbulent', color='orange', linestyle='-')
plt.plot(times, max_radius_norm, label='normal', color='blue', linestyle='-.')
plt.plot(times, max_radius_turb, label='turbulent', color='orange', linestyle='-.')
plt.xlabel("Time (Myr)")
plt.ylabel("BE radius (pc)")
plt.yscale('log')
plt.ylim(5e-1, 5e2)
plt.legend(loc='upper left')
plt.savefig("be_radius_vs_time.png")
plt.close()

plt.figure()
plt.title("Bonnor-Ebert mass")
plt.plot(times, intcp_mass_norm, label='normal', color='blue', linestyle='-')
plt.plot(times, intcp_mass_turb, label='turbulent', color='orange', linestyle='-')
plt.plot(times, max_mass_norm, label='normal', color='blue', linestyle='-.')
plt.plot(times, max_mass_turb, label='turbulent', color='orange', linestyle='-.')
plt.xlabel("Time (Myr)")
plt.ylabel(r"$BE mass (M_{\odot})$")
plt.yscale('log')
plt.ylim(1e2, 5e6)
plt.legend(loc='upper left')
plt.savefig("be_mass_vs_time.png")
plt.close()
