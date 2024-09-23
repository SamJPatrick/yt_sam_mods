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
from matplotlib.collections import LineCollection
from matplotlib.colorbar import ColorbarBase


import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc, kb, mh



NUM_RAYS = 20
X_LIM = unyt_quantity(400.0, 'pc')

PRE_DIR = "Britton_sim_data"
DISTANCE_FILE = "halo_distances.txt"
SIM_FILE = "simulation.h5"
RAY_DIR = "Sobolev/Ray_profiles"

dataset_dict = {'ccsn': [106, 110, 113, 119],
                'hn': [151, 156, 159, 165],
                'pisn': [110, 140, 150, 160]}


CMAP = "viridis"
cmap = plt.get_cmap(CMAP)
#cmap = matplotlib.colormaps[CMAP]
COLORS = cmap(np.linspace(0, 1, len(dataset_dict['pisn'])))
#COLORS = ['magenta', 'cyan', 'blue', 'green', 'orange', 'brown']
FIELDS = ['velocity_para', 'temperature', 'density', 'metallicity3']

my_fig = GridFigure(4, 3, figsize=(11, 11),
                    left_buffer=0.08, right_buffer=0.08,
                    bottom_buffer=0.06, top_buffer=0.06,
                    horizontal_buffer=0.06, vertical_buffer=0.02)
assert len(my_fig) == len(dataset_dict) * len(FIELDS), \
    "Error, number of figure pannels doesn't match number of stars * fields"


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




for i, star_type in enumerate(list(dataset_dict.keys())):

    with open(os.path.join(PRE_DIR, star_type.upper(), DISTANCE_FILE), newline='\n') as myfile:
        reader = list(csv.reader(myfile, delimiter= '\t'))
    dumps, distances = zip(*[entry for entry in reader])
    distances = unyt_array([float(distance) for distance in distances], 'pc')
    times = transpose_unyt([get_time_z(dump, star_type, sim_file= os.path.join(PRE_DIR, star_type.upper(), SIM_FILE))[0] for dump in dumps])
    index = np.argmin(np.abs(times))
    distances = distances[index:]
    times = times[index:]
    #x_halo = [distances[np.where(f"DD{num:04d}" == dumps)[0][0]] for num in dataset_dict[star_type]]

    for k, field in enumerate(FIELDS):
        
        gnum = 3 * k + i 
        my_y = get_field_dict(field)['limits'][0] * np.ones(distances.size)
        points = np.array([distances, my_y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        #my_norm = plt.Normalize(times[0], times[-1])
        if (star_type == 'pisn'):
            my_norm = plt.Normalize(0.0, 5.0)
        else :
            my_norm = plt.Normalize(0.0, 7.0)
        lc = LineCollection(segments, cmap= CMAP, norm=my_norm)
        lc.set_array(times)
        lc.set_linewidth(20)
        #my_line = my_axes.add_collection(lc)
        my_fig[gnum].add_collection(lc)

        if (k == 0):
            new_ax = my_fig.add_cax(my_fig[gnum], "right", buffer=0.01, width=0.02)
            cbar = ColorbarBase(new_ax, cmap= plt.get_cmap(CMAP), norm=my_norm, orientation='vertical')
            #cbar.set_label("Time (Myr)")
            cbar.set_label("")


for i, my_axes in enumerate(my_fig.left_axes):
    field = FIELDS[i]
    field_dict = get_field_dict(field)
    if ('velocity' in field):
        my_axes.yaxis.set_label_text(r"Velocity (km s$^{-1}$)")
    elif ('metallicity' in field):
        my_axes.yaxis.set_label_text(r"Metallicity (Z$_{\odot}$)")
    elif ('density' in field):
        my_axes.yaxis.set_label_text(r"Density (g cm$^{-3}$)")
    else :
        label = ' '.join(np.char.capitalize(field.split('_')))
        my_axes.yaxis.set_label_text(f"{label} ({field_dict['units']})")
        
for my_axes in my_fig.right_axes:
    my_axes.tick_params(labelleft=False)
for my_axes in my_fig.center_axes:
    my_axes.tick_params(labelleft=False)
for my_axes in my_fig.middle_axes:
    my_axes.tick_params(labelbottom=False)
    
for my_axes in my_fig.bottom_axes:
    my_axes.xaxis.set_ticks(np.arange(0, 400, 50), minor=True)
    my_axes.xaxis.set_label_text("Radius (pc)")

for i, my_axes in enumerate(my_fig.top_axes):
    my_axes.xaxis.set_label_position("top")
    my_axes.xaxis.set_label_text(star_type.upper(), fontsize=14)
    my_axes.tick_params(labelbottom=False)
    my_axes.legend(loc='upper right')


'''
for i, my_axes in enumerate(my_fig.top_axes):

    star_type = list(dataset_dict.keys())[i]
    
    my_axes.xaxis.set_label_position("top")
    my_axes.xaxis.set_label_text(star_type.upper(), fontsize=14)
    my_axes.tick_params(labelbottom=False)
    my_axes.legend(loc='upper right')

    with open(os.path.join(PRE_DIR, star_type.upper(), DISTANCE_FILE), newline='\n') as myfile:
        reader = list(csv.reader(myfile, delimiter= '\t'))
    dumps, distances = zip(*[entry for entry in reader])
    distances = unyt_array([float(distance) for distance in distances], 'pc')
    times = transpose_unyt([get_time_z(dump, star_type, sim_file= os.path.join(PRE_DIR, star_type.upper(), SIM_FILE))[0] for dump in dumps])
    index = np.argmin(np.abs(times))
    distances = distances[index:]
    times = times[index:]
    #x_halo = [distances[np.where(f"DD{num:04d}" == dumps)[0][0]] for num in dataset_dict[star_type]]

    my_y = get_field_dict(FIELDS[0])['limits'][0] * np.ones(distances.size)
    points = np.array([distances, my_y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #my_norm = plt.Normalize(times[0], times[-1])
    if (star_type == 'pisn'):
        my_norm = plt.Normalize(0.0, 5.0)
    else :
        my_norm = plt.Normalize(0.0, 7.0)
    lc = LineCollection(segments, cmap= CMAP, norm=my_norm)
    lc.set_array(times)
    lc.set_linewidth(20)
    #my_line = my_axes.add_collection(lc)
    my_axes.add_collection(lc)
    
    new_ax = my_fig.add_cax(my_axes, "right", buffer=0.01, width=0.02)
    cbar = ColorbarBase(new_ax, cmap= plt.get_cmap(CMAP), norm=my_norm, orientation='vertical')
    #cbar.set_label("Time (Myr)")
    cbar.set_label("")

for my_axes in my_fig:
    my_axes.add_collection(lc)
'''    

plt.savefig("model_profiles.png")
