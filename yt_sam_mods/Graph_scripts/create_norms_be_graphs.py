import os
import h5py
import glob
import yt
from unyt import Msun, pc
from unyt import unyt_array
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import get_title, get_dump_num, get_time_z
from yt.extensions.sam_mods.misc import transpose_unyt



BE_FOLDER = "Profiles/Bonnor_ebert_profiles"
BE_FILE = "be_time_mass_radius.h5"
NORM_FOLDER = "Profiles/Normal_profiles"
OUTFILE = "be_norms_vs_time.h5"

'''
be_info = h5py.File(os.path.join(os.getcwd(), BE_FOLDER, BE_FILE), 'r')
be_times = unyt_array(list(be_info['time/array_data']), 'Myr')
be_radii = unyt_array(list(be_info['max_radius_norm/array_data']), 'pc')
filenames = glob.glob(os.path.join(os.getcwd(), NORM_FOLDER, "DD*_profile_weight_field_cell_mass.h5"))
filenames.sort(reverse= True, key= lambda x: get_time_z(os.path.basename(x))[0])
assert len(be_times) == len(be_radii) == len(filenames), "Error, number of times and BE radii should be equal"

densities = [0.0] * len(be_times)
temperatures = [0.0] * len(be_times)
metallicities = [0.0] * len(be_times)
densities_final = [0.0] * len(be_times)
temperatures_final = [0.0] * len(be_times)
metallicities_final = [0.0] * len(be_times)

for fn in filenames:
    
    df = yt.load(fn)
    used = df.profile.used
    radii = df.profile.x[used]
    i, time_diff = min(enumerate(np.abs(be_times - get_time_z(os.path.basename(fn))[0])), key= lambda x: x[1])
    if (i == len(be_times) - 1):
        be_radius_final, rad_diff = min([(radius, np.abs(radius - be_radii[i])) for radius in radii], key= lambda x: x[1])
    rad_index, rad_diff = min(enumerate(np.abs(radii - be_radii[i])), key= lambda x: x[1])
    
    densities[i] = df.profile[('data', 'density')][used][rad_index]
    temperatures[i] = df.profile[('data', 'temperature')][used][rad_index]
    metallicities[i] = df.profile[('data', 'metallicity3')][used][rad_index]
    try :
        rad_index_final, rad_diff = min(enumerate(np.abs(be_radius_final - radii)), key= lambda x: x[1])
    except NameError:
        print("Error, could not obtain BE radius in final data dump")
    densities_final[i] = df.profile[('data', 'density')][used][rad_index_final]
    temperatures_final[i] = df.profile[('data', 'temperature')][used][rad_index_final]
    metallicities_final[i] = df.profile[('data', 'metallicity3')][used][rad_index_final]
    

densities = transpose_unyt(densities)
temperatures = transpose_unyt(temperatures)
metallicities = transpose_unyt(metallicities)
densities_final = transpose_unyt(densities_final)
temperatures_final = transpose_unyt(temperatures_final)
metallicities_final = transpose_unyt(metallicities_final)
be_times.write_hdf5(OUTFILE, group_name='time')
densities.write_hdf5(OUTFILE, group_name='density')
temperatures.write_hdf5(OUTFILE, group_name='temperature')
metallicities.write_hdf5(OUTFILE, group_name='metallicity')
densities_final.write_hdf5(OUTFILE, group_name='density_final')
temperatures_final.write_hdf5(OUTFILE, group_name='temperature_final')
metallicities_final.write_hdf5(OUTFILE, group_name='metallicity_final')
'''

df = h5py.File(OUTFILE, 'r')
be_times = df['time/array_data']
densities = df['density/array_data']
temperatures = df['temperature/array_data']
metallicities = df['metallicity/array_data']
densities_final = df['density_final/array_data']
temperatures_final = df['temperature_final/array_data']
metallicities_final = df['metallicity_final/array_data']

plt.figure()
plt.title("Density at BE radius vs time")
plt.plot(be_times, densities, color='blue', label='current')
plt.plot(be_times, densities_final, color='orange', label='final')
plt.xlabel("Time (Myr)")
plt.ylabel(r"$\rho ~(g ~cm^{-3})$")
plt.yscale('log')
plt.ylim(1e-25, 1e-19)
plt.legend(loc= 'upper left')
plt.savefig("be_density_vs_time.png")
plt.close()

plt.figure()
plt.title("Temperature at BE radius vs time")
plt.plot(be_times, temperatures, color='blue', label='current')
plt.plot(be_times, temperatures_final, color='orange', label='final')
plt.xlabel("Time (Myr)")
plt.ylabel(r"$T ~(K)$")
plt.yscale('log')
plt.ylim(1e2, 1e4)
plt.legend(loc= 'upper left')
plt.savefig("be_temperature_vs_time.png")
plt.close()

plt.figure()
plt.title("Metallicity at BE radius vs time")
plt.plot(be_times, metallicities, color='blue', label='current')
plt.plot(be_times, metallicities_final, color='orange', label='final')
plt.xlabel("Time (Myr)")
plt.ylabel(r"$Z ~(Z_{\odot})$")
plt.yscale('log')
plt.ylim(1e-10, 1e-2)
plt.legend(loc= 'upper left')
plt.savefig("be_metallicty_vs_time.png")
plt.close()
