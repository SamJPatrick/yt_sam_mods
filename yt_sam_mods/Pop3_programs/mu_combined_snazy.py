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
from grid_figure import GridFigure
from unyt import unyt_quantity, unyt_array
from unyt import kb, mh




NUM_RAYS = 20
FILTER_WINDOW = 7
NU_LIMS = unyt_array([1e29, 1e38], 'm**2/s')

PRE_DIR = "Britton_sim_data"
POST_DIR = "Sobolev/Ray_profiles"
SIM_FILE = "simulation.h5"

COLORS = ['blue', 'orange', 'magenta', 'cyan', 'brown', 'red']
time_dict = {'ccsn': [106, 110, 119, 137],
             'hn': [151, 156, 165, 183],
             'pisn': [110, 150, 175, 205]}
MU = 1.6
GAMMA = 5/3
CONV_CONST = kb / (MU * mh)


my_fig = GridFigure(3, 1, figsize=(11, 9),
                    left_buffer=0.09, right_buffer=0.08,
                    bottom_buffer=0.1, top_buffer=0.1,
                    vertical_buffer=0, horizontal_buffer=0.12)


for i, star_type in enumerate(list(time_dict.keys())):
    ray_files = np.sort([os.path.join(PRE_DIR, star_type.upper(), POST_DIR, f"DD{num:04d}_packed.h5")
                         for num in time_dict[star_type]])
    for j, ray_file in enumerate(ray_files):

        df = h5py.File(ray_file, 'r')
        distances = unyt_array(df['distances/array_data'], 'pc')
        arr_len = len(distances)
    
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
        my_fig[i].plot(distances[:-1], signal.medfilt(nu, kernel_size= FILTER_WINDOW), color= COLORS[j], label=f"{time:.2f}")

    my_fig[i].set_ylim(*NU_LIMS)
    my_fig[i].set_yscale("log")
    my_fig[i].yaxis.set_label_text(r"$\nu_{eff} ~$" + f"{star_type}" + r"$~$(m$^{2}$ s$^{-1}$)")
    my_fig[i].set_xlim(0.0, 400.0)
    my_fig[i].grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
    my_fig[i].tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)
    my_fig[i].legend(loc='upper right')

my_fig[0].xaxis.set_ticklabels(['' for i in range (len(np.arange(0, 401, 20)))])
my_fig[1].xaxis.set_ticklabels(['' for i in range (len(np.arange(0, 401, 20)))])
my_fig[2].xaxis.set_ticks(np.arange(0, 401, 20), minor=True)
my_fig[2].xaxis.set_label_text(r"Distance (pc)")
plt.savefig("model_profiles.png")
plt.close()
