import os
import sys
import csv
import glob
import h5py
import numpy as np

from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from grid_figure import GridFigure

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc, kb, mh



NUM_RAYS = 20
X_LIM = unyt_quantity(400.0, 'pc')

PRE_DIR = "Britton_sim_data"
RAY_DIR = "Sobolev/Ray_profiles"
#dataset_dict = {'ccsn': [106, 110, 119],
#                'hn': [151, 156, 165],
#                'pisn': [110, 150, 175]}
#dataset_dict = {'ccsn': [106, 113, 119],
#                'hn': [151, 159, 165],
#                'pisn': [110, 140, 160]}
dataset_dict = {'ccsn': [106, 110, 113, 119],
                'hn': [151, 156, 159, 165],
                'pisn': [110, 140, 150, 160]}



COLORS = ['green', 'blue', 'orange', 'magenta', 'cyan', 'brown']
FIELDS = ['velocity_para', 'temperature', 'metallicity3']

my_fig = GridFigure(3, 3, figsize=(11, 9),
                    left_buffer=0.08, right_buffer=0.02,
                    bottom_buffer=0.06, top_buffer=0.08,
                    horizontal_buffer=0.02, vertical_buffer=0.02)
assert len(my_fig) == len(dataset_dict) * len(FIELDS), \
    "Error, number of figure pannels doesn't match number of stars * fields"

'''
#FIELDS = ['density', 'El_fraction', 'temperature', 'velocity_para']
#DATASET_NUMS = [100, 110]
#DATASET_NUMS = [120, 130, 140, 150, 160]
#DATASET_NUMS = [139, 141, 145, 151]
#DATASET_NUMS = [153, 157, 161, 173, 185]
#DATASET_NUMS = [100, 101, 103, 106]
#DATASET_NUMS = [107, 111, 115, 127, 139]
#COLORS = ['blue', 'orange', 'magenta', 'cyan', 'brown', 'red']
'''

for i, star_type in enumerate(list(dataset_dict.keys())):

    data_rays = np.sort([os.path.join(PRE_DIR, star_type.upper(), RAY_DIR, f"DD{num:04d}_packed.h5") for num in dataset_dict[star_type]])
    times = transpose_unyt([get_time_z(os.path.basename(filename), star_type,
                                       sim_file= os.path.join(PRE_DIR, star_type.upper(), "simulation.h5"))[0].to('Myr')for filename in data_rays])

    for k, field in enumerate(FIELDS):
        field_dict = get_field_dict(field)
        gnum = 3 * k + i 
        for j in range (len(data_rays)):
            
            try :
                ds_rays = h5py.File(data_rays[j], 'r')
            except Exception:
                print(f"Error, could not load dataset {data_rays[j]}, continuing")
                continue
            distances_ray = unyt_array(ds_rays["distances/array_data"], 'pc')
            arr_ray = transpose_unyt([np.nanmean(transpose_unyt(x)) for x in \
                                      zip(*[unyt_array(ds_rays[f'{field}_{n}/array_data'], field_dict['units']) for n in range (NUM_RAYS)])])
            my_fig[gnum].plot(distances_ray, arr_ray, color= COLORS[j], label= f'{times[j]:.2f}')
        
        if (field_dict['log'] == True):
            my_fig[gnum].set_yscale('log')
        if ('velocity' in field):
            my_fig[gnum].axhline(y= 0.0, color='black')
        my_fig[gnum].tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)
        my_fig[gnum].grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
        my_fig[gnum].set_xlim(0.0, X_LIM)
        my_fig[gnum].set_ylim(field_dict['limits'])


for i, my_axes in enumerate(my_fig.left_axes):
    field = FIELDS[i]
    field_dict = get_field_dict(field)
    if ('velocity' in field):
        my_axes.yaxis.set_label_text(r"Velocity (km s$^{-1}$)")
    elif ('metallicity' in field):
        my_axes.yaxis.set_label_text(r"Metallicity (Z$_{\odot}$)")
    else :
        label = ' '.join(np.char.capitalize(field.split('_')))
        my_axes.yaxis.set_label_text(f"{label} ({field_dict['units']})")
for my_axes in my_fig.right_axes:
    my_axes.tick_params(labelleft=False)
for my_axes in my_fig.center_axes:
    my_axes.tick_params(labelleft=False)
for my_axes in my_fig.bottom_axes:
    my_axes.xaxis.set_ticks(np.arange(0, 400, 50), minor=True)
    my_axes.xaxis.set_label_text("Radius (pc)")
for i, my_axes in enumerate(my_fig.top_axes):
    my_axes.xaxis.set_label_position("top")
    my_axes.xaxis.set_label_text(list(dataset_dict.keys())[i].upper(), fontsize=14)
    my_axes.tick_params(labelbottom=False)
    my_axes.legend(loc='upper right')
for my_axes in my_fig.middle_axes:
    my_axes.tick_params(labelbottom=False)

plt.savefig("model_profiles.png")
