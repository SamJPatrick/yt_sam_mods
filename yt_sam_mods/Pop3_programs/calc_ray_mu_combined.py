import glob
import csv
import h5py
import sys
import os

import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from unyt import unyt_quantity, unyt_array
from unyt import kb, mh
import yt
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from yt.extensions.sam_mods.graph_funcs import *




NUM_RAYS = 20
FILTER_WINDOW = 7

INDIR_RAYS = "Sobolev/Ray_profiles"
OUTDIR_NU = "Sobolev/"
DISTANCE_FILE = "ray_distances.txt"

COLORS = ['blue', 'orange', 'magenta', 'cyan', 'brown', 'red']
#DATASET_NUMS = [120, 130, 140, 150, 160]
#DATASET_NUMS = [165, 175, 190]
DATASET_NUMS = [153, 157, 161, 173, 185]

MU = 1.6
GAMMA = 5/3
CONV_CONST = kb / (MU * mh)
Z0 = unyt_quantity(1e-7, 'Zsun')
NU_LIMS = unyt_array([1e29, 1e38], 'm**2/s')


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


with open(DISTANCE_FILE, newline='\n') as myfile:
    reader = list(csv.reader(myfile, delimiter= '\t'))
dumps, distances = zip(*[entry for entry in reader])
distances = unyt_array([float(distance) for distance in distances], 'pc')
dumps = np.array(dumps)
x_halo = [distances[np.where(f"DD{num:04d}" == dumps)[0][0]] for num in DATASET_NUMS]


ray_files = sorted([os.path.join(INDIR_RAYS, f"DD{num:04d}_packed.h5") for num in DATASET_NUMS])
x_front = [0.0 for i in range (len(DATASET_NUMS))]
plt.figure()
plt.title(star_type.upper())
for i, ray_file in enumerate(ray_files):
    
    df = h5py.File(ray_file, 'r')
    dump_name = os.path.basename(ray_file)
    distances = unyt_array(df['distances/array_data'], 'pc')
    arr_len = len(distances)
    time = get_time_z(dump_name, star_type)[0].to('Myr')

    mean_metals = unyt_array([np.mean(x) for x in zip(*[df[f'metallicity3_{n}/array_data'] for n in range (NUM_RAYS)])], 'Zsun')
    indices = np.where(mean_metals < Z0)[0]
    if (len(indices) == 0):
        x_front[i] = x_front[i-1]
        continue
    for index in indices:
        if (i==0 or distances[index] >= x_front[i-1]):
            x_front[i] = distances[index]
            break
        if (x_front[i] == 0.0):
            x_front[i] = x_front[i-1]
    
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
    plt.plot(distances[:-1], signal.medfilt(nu, kernel_size= FILTER_WINDOW), color= COLORS[i], label=f"{time:.2f}")

plt.scatter(x_halo, [NU_LIMS[0] for i in range (len(x_halo))], color=COLORS[:len(x_halo)], marker= 'X')
plt.scatter(x_front, [NU_LIMS[0] for i in range (len(x_front))], color=COLORS[:len(x_front)], marker= 'o')    
plt.legend(loc='upper right')
plt.ylabel(r"$\nu ~(m^{2} ~s^{-1})$")
plt.xlabel("Radius (pc)")
plt.xlim(0.0, 400.0)
plt.yscale('log')
plt.ylim(*NU_LIMS)
plt.savefig(os.path.join(OUTDIR_NU, "nu_combined_st.png"))
plt.close()
