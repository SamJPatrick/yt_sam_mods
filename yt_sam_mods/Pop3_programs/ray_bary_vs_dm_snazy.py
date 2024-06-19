import os
import yt
import sys
import glob
import h5py
import numpy as np
import scipy.signal as sig

from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from grid_figure import GridFigure

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc



NUM_RAYS = 20
FILT_WINDOW = 11
DENS_CUTOFF = 1e-27

POST_DIR = "Britton_sim_data"
RAY_DIR = "Sobolev/Ray_profiles"
dataset_dict = {'pisn': [98, 100, 110, 120],
                #'hn': [137, 139, 141, 145, 151],
                'ccsn': [99, 100, 101, 103, 106]}
COLORS = ['green', 'blue', 'orange', 'magenta', 'cyan', 'brown']

my_fig = GridFigure(2, 1, figsize=(8, 6),
                    left_buffer=0.10, right_buffer=0.02,
                    bottom_buffer=0.08, top_buffer=0.02,
                    vertical_buffer=0.00)


for i, star_type in enumerate(list(dataset_dict.keys())):
    datasets = np.sort([os.path.join(POST_DIR, star_type.upper(), RAY_DIR, f"DD{num:04d}_packed.h5") \
                        for num in dataset_dict[star_type]])
    for j, dataset in enumerate(datasets):
    
        dump_name = os.path.basename(dataset)
        time = get_time_z(dump_name, star_type, sim_file= os.path.join(POST_DIR, star_type.upper(), "simulation.h5"))[0].to('Myr')
        ds = h5py.File(dataset, 'r')
        distances = unyt_array(ds['distances/array_data'], 'pc')
        bary_densities = transpose_unyt([np.mean(transpose_unyt(x)) for x in \
                                         zip(*[unyt_array(ds[f'density_{n}/array_data'], 'g/cm**3') for n in range (NUM_RAYS)])])

        dm_arrs = list(zip(*[unyt_array(ds[f'dark_matter_density_{n}/array_data'], 'g/cm**3') for n in range (NUM_RAYS)]))
        dm_densities = [0.0 for j in range (len(dm_arrs))]
        for k, arr in enumerate(dm_arrs):
            arr_sane = [k for k in arr if k > unyt_quantity(DENS_CUTOFF, 'g/cm**3')]
            if (len(arr_sane) == 0):
                dm_mean = unyt_quantity(0.0, 'g/cm**3')
            else :
                dm_mean = np.mean(transpose_unyt(arr_sane))
            dm_densities[k] = dm_mean
        dm_densities = transpose_unyt(dm_densities)

        my_fig[i].plot(distances, sig.medfilt(dm_densities, kernel_size= FILT_WINDOW), color= COLORS[j], linestyle='--')
        my_fig[i].plot(distances, bary_densities, color= COLORS[j], linestyle='-', label=f'{time:.2f}')

    my_fig[i].yaxis.set_label_text(r"$\rho~$" + f"{star_type}" + r"$~$(g cm$^{-3}$)")
    my_fig[i].set_yscale('log')
    my_fig[i].set_ylim(1e-27, 1e-18)
    my_fig[i].set_xlim(0.0, 400.0)
    my_fig[i].xaxis.set_ticks(np.arange(0, 401, 20), minor= True)
    my_fig[i].grid(visible=True, axis="both", zorder=1, linestyle=":", color="black", alpha=0.6)
    my_fig[i].tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)
    my_fig[i].legend(loc='upper right')

my_fig[0].xaxis.set_ticklabels(['' for n in range(len(np.arange(0, 401, 20)))])
my_fig[1].xaxis.set_label_text("Radius (pc)")
plt.savefig("model_profiles.png")
