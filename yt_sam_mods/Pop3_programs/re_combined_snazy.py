import h5py
import sys
import os

import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import yt
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from yt.extensions.sam_mods.graph_funcs import *
from unyt import unyt_quantity, unyt_array
from unyt import kb, mh, pc


NUM_RAYS = 20
Z0 = unyt_quantity(1e-7, 'Zsun') 

PRE_DIR = "Britton_sim_data"
POST_DIR = "Sobolev/Ray_profiles"
SIM_FILE = "simulation.h5"

color_dict = {'ccsn': 'blue', 'hn': 'orange', 'pisn': 'green'}
time_dict = {'ccsn': 119, 'hn': 165,'pisn': 175}
width_dict = {'ccsn': 17.0 * pc, 'hn': 12.0 * pc,'pisn': 33.0 * pc}

MU = 1.6
GAMMA = 5/3
CONV_CONST = kb / (MU * mh)

plt.figure()
for star_type in list(time_dict.keys()):
    
    ray_file = os.path.join(PRE_DIR, star_type.upper(), POST_DIR, f"DD{time_dict[star_type]:04d}_packed.h5")
    df = h5py.File(ray_file, 'r')
    distances = unyt_array(df['distances/array_data'], 'pc')
    arr_len = len(distances)

    mean_metals = unyt_array([np.mean(x) for x in zip(*[df[f'metallicity3_{n}/array_data'] for n in range (NUM_RAYS)])], 'Zsun')
    x_front = distances[np.argwhere(mean_metals < Z0)[0].item()]
      
    v_para = [transpose_unyt(x) for x in zip(*[unyt_array(df[f'velocity_para_{n}/array_data'], 'km/s') for n in range (NUM_RAYS)])]
    v_rad = [transpose_unyt(x) for x in zip(*[unyt_array(df[f'velocity_rad_{n}/array_data'], 'km/s') for n in range (NUM_RAYS)])]
    v_phi = [transpose_unyt(x) for x in zip(*[unyt_array(df[f'velocity_phi_{n}/array_data'], 'km/s') for n in range (NUM_RAYS)])]
    temps = [transpose_unyt(x) for x in zip(*[unyt_array(df[f'temperature_{n}/array_data'], 'K') for n in range (NUM_RAYS)])]
    del_phi = np.ptp(df['angle_list/array_data'])

    v_para_ave = transpose_unyt([np.nanmean(vels) for vels in v_para])
    v_rad_ave = transpose_unyt([np.nanmean(vels) / (distances[j] * del_phi) for j, vels in enumerate(v_rad)])
    v_phi_ave = transpose_unyt([np.nanmean(vels) / (distances[j] * del_phi) for j, vels in enumerate(v_phi)])
    t_ave = transpose_unyt([np.nanmean(temperature) for temperature in temps])

    v_para_std = transpose_unyt([np.nanstd(vels) / (distances[j] * del_phi) for j, vels in enumerate(v_para)])
    v_rad_std = transpose_unyt([np.nanstd(vels) / (distances[j] * del_phi) for j, vels in enumerate(v_rad)])
    v_phi_std = transpose_unyt([np.nanstd(vels) / (distances[j] * del_phi) for j, vels in enumerate(v_phi)])
    t_std = transpose_unyt([np.nanstd(temperature) for temperature in temps])

    v_para_diff = transpose_unyt([np.nanmean(transpose_unyt([v_para[j+1][n] - v_para[j][n] for n in range (NUM_RAYS)])) / (distances[j+1] - distances[j]) for j in range (arr_len - 1)])
    v_rad_diff = transpose_unyt([np.nanmean(transpose_unyt([v_rad[j+1][n] - v_rad[j][n] for n in range (NUM_RAYS)])) / (distances[j+1] - distances[j]) for j in range (arr_len - 1)])
    v_phi_diff = transpose_unyt([np.nanmean(transpose_unyt([v_phi[j+1][n] - v_phi[j][n] for n in range (NUM_RAYS)])) / (distances[j+1] - distances[j]) for j in range (arr_len - 1)])
    t_diff = transpose_unyt([np.nanmean(transpose_unyt([temps[j+1][n] - temps[j][n] for n in range (NUM_RAYS)])) / (distances[j+1] - distances[j]) for j in range (arr_len - 1)])

    cooling_times = transpose_unyt([np.nanmean(transpose_unyt(times)) for times in zip(*[unyt_array(df[f'cooling_time_{n}/array_data'], 's') for n in range (NUM_RAYS)])])
    densities = transpose_unyt([np.nanmean(transpose_unyt(dens)) for dens in zip(*[unyt_array(df[f'density_{n}/array_data'], 'g/cm**3') for n in range (NUM_RAYS)])])
    numer = (1/(GAMMA-1)) * ((t_ave[:-1] / cooling_times[:-1]) + v_rad_ave[:-1] * t_std[:-1] + v_para_ave[:-1] * t_diff) + t_ave[:-1] * (v_rad_std[:-1] + v_rad_ave[:-1] + v_para_diff)
    denom = v_rad_std[:-1]**2 + v_para_diff**2 + 0.5 * v_rad_diff**2 + 0.5 * (v_phi_std[:-1] - v_phi_ave[:-1])**2 + 0.5 * (v_rad_diff + v_para_std[:-1])**2 + 0.5 * v_phi_diff**2 
    nu = (0.5 * CONV_CONST * (numer / denom)).to('m**2/s')

    time = get_time_z(os.path.basename(ray_file), star_type, sim_file= os.path.join(PRE_DIR, star_type.upper(), SIM_FILE))[0].to('Myr')
    index_span = np.argwhere((distances - width_dict[star_type]) > 0)[0].item()
    index_start = int(index_span / 2)
    index_end = arr_len - int(index_span / 2) - 1
    distances_len = distances[index_start:index_end]
    integrals = [0.0 for distance in distances_len]
    for k, index in enumerate(range(index_start, index_end)):
        for l in range (index - int(index_span / 2), index + int(index_span / 2)):
            integrals[k] += (distances[l+1] - distances[l]) / nu[l] 
    plt.plot(distances_len, integrals, color= color_dict[star_type], label=f"{star_type} @ {time:.2f}" + r"$~(\omega=$" + f"{width_dict[star_type]})")
    plt.axvline(x= x_front, linestyle= '--', color= color_dict[star_type])

plt.ylabel(r"$ Re/v ~(s ~m^{-1})$")
plt.xlabel("Distance (pc)")
plt.xlim(0.0, 400.0)
plt.yscale('log')
plt.ylim(1e-19, 1e-10)
plt.grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
plt.tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)
plt.xticks(np.arange(0, 401, 50))
plt.legend(loc='upper left')
plt.savefig("model_profiles.png")
plt.close()
