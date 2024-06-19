import os
import sys
import h5py
import glob

from unyt import Msun, pc, Myr, pc, kb, mh, yr
from unyt import unyt_array
import numpy as np

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from grid_figure import GridFigure


X_LIM = unyt_quantity(400.0, 'pc')
N_H2 = unyt_quantity(1e4, 'cm**(-3)')

PROF_DIR = "."
COLORS = ['green', 'blue', 'orange', 'magenta', 'cyan', 'brown']
FIELDS_ENC = {'mass_enc': 'Msun', 'energy_enc': 'erg'}
FIELDS_TIME = ['free-fall', 'sound-crossing'] 

prof_files = glob.glob("*_clump_density_profiles.h5")
num_stars = len(prof_files)
my_fig = GridFigure(2, num_stars, figsize=(11, 9),
                    left_buffer=0.08, right_buffer=0.08,
                    bottom_buffer=0.06, top_buffer=0.02,
                    horizontal_buffer=0.08, vertical_buffer=0.02)


for i, prof_file in enumerate(prof_files):
    
    df = h5py.File(prof_files[i], 'r')
    number_densities = unyt_array(list(df['density/array_data']), 'cm**(-3)')
    star_type = prof_file.split('_')[0]
    
    for j, field in enumerate(FIELDS_TIME):
        data = unyt_array(list(df[f'{field}/array_data']), 'Myr')
        my_fig[i*2].plot(number_densities, data, color= COLORS[j], label=FIELDS_TIME[j])
    my_fig[i*2].set_yscale('log')
    my_fig[i*2].yaxis.set_label_text(f"Timescale {star_type} (Myr)")
    my_fig[i*2].set_ylim(1e-4, 1e3)
    my_fig[i*2].yaxis.set_ticks(np.logspace(-4, 3, 8), minor=True)
    my_fig[i*2].axvline(x= N_H2, color='black', linestyle='--')
    my_fig[i*2].legend(loc='upper right')

    my_axes = [my_fig[i*2+1], my_fig[i*2+1].twinx()]
    for j, field in enumerate(list(FIELDS_ENC.keys())):
        data = unyt_array(list(df[f'{field}/array_data']), FIELDS_ENC[field])
        my_axes[j].plot(number_densities, data, color= COLORS[j], label= list(FIELDS_ENC.keys())[j])
    my_axes[0].set_yscale('log')
    my_axes[0].yaxis.set_label_text(r"Enclosed mass %s (M$_{\odot}$)" % star_type)
    my_axes[0].set_ylim(1e-4, 1e6)
    my_axes[0].yaxis.set_ticks(np.logspace(-4, 6, 11), minor=True)
    my_axes[0].axvline(x= N_H2, color='black', linestyle='--')
    lines, labels = my_axes[0].get_legend_handles_labels()
    my_axes[1].set_yscale('log')
    my_axes[1].yaxis.set_label_text(f"Enclosed thermal energy {star_type} (erg)")
    my_axes[1].set_ylim(1e40, 1e50)
    my_axes[1].yaxis.set_ticks(np.logspace(40, 50, 11), minor=True)
    my_axes[1].axvline(x= N_H2, color='black', linestyle='--')
    lines1, labels1 = my_axes[1].get_legend_handles_labels()
    my_axes[1].legend(lines + lines1, labels + labels1, loc='upper right')


for my_axes in my_fig:
    my_axes.set_xscale('log')
    my_axes.set_xlim(1e-5, 1e11)
    my_axes.xaxis.set_ticks(np.logspace(-5, 11, 17), minor=True)
    my_axes.tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)
    my_axes.grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
for my_axes in my_fig.left_axes:
    my_axes.tick_params(axis="y", left=True, direction="inout", which="both")
    my_axes.tick_params(axis="y", right=True, direction="in", which="both")
for my_axes in my_fig.right_axes:
    my_axes.tick_params(axis="y", right=True, direction="inout", which="both", labelright=True)
    my_axes.tick_params(axis="y", left=True, direction="inout", which="both", labelleft=True)
for my_axes in my_fig.bottom_axes:
    my_axes.xaxis.set_label_text(r"Number density (cm$^{-3}$)")
for my_axes in my_fig.top_axes:
    my_axes.tick_params(axis="x", top=True, direction="inout", which="both", labeltop=False)
    my_axes.tick_params(axis="x", bottom=True, direction="inout", which="both", labelbottom=False)
    #my_axes.set_xlabel(xlabel, labelpad=8)
try :
    for my_axes in my_fig.middle_axes:
        my_axes.tick_params(axis="x", top=True, bottom=True, labeltop=False, labelbottom=False)
except RuntimeError:
    pass

plt.savefig("model_profiles.png")
