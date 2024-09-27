import sys
import glob
import os
import h5py
import yt
import numpy as np

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from grid_figure import GridFigure

from unyt import Msun, pc
from unyt import unyt_array
from yt.extensions.sam_mods.graph_funcs import *


PRE_DIR = "Britton_sim_data"
INDIR = "Profiles/Pressure_profiles"
OUTDIR = "Profiles/Pressure_graphs_test"
FIELDS = ['hydrostatic', 'thermal', 'ram', 'turbulent']
FIELD_DICT = {'hydrostatic': ('blue', '-'),
              'thermal': ('red', '--'),
              'ram': ('green', '-'),
              'turbulent': ('orange', '--')}
DATA_DICT = {'ccsn': [151, 163, 168],
             'hn': [197, 217, 224],
             'pisn': [278, 284, 295]}


my_fig = GridFigure(3, 3, figsize=(13, 9),
                    left_buffer=0.08, right_buffer=0.02,
                    bottom_buffer=0.06, top_buffer=0.06,
                    vertical_buffer=0.04, horizontal_buffer=0.06)


for i, (star_type, dir_nums) in enumerate(list(DATA_DICT.items())):
    for j, dir_num in enumerate(dir_nums):

        gnum = 3 * j + i
        filename = os.path.join(PRE_DIR, star_type.upper(), INDIR, f'DD{dir_num:04d}_pressures.h5')
        sim_file = os.path.join(PRE_DIR, star_type.upper(), "simulation.h5")
        prof_df = h5py.File(filename, 'r')
        time = get_time_z(os.path.basename(filename), star_type, sim_file= sim_file)[0].to('Myr')
        df_be = yt.load(os.path.join(PRE_DIR, star_type.upper(), "star_None_mass_normal_profiles.h5"))
        ind_time = np.argmin(np.abs(df_be.data[('data', 'time')].to('Myr') - get_time_offset(star_type) - time))
        ind_max = np.argmax(df_be.data[('data', 'bonnor_ebert_ratio')][ind_time])
        radius_be = df_be.data[('data', 'radius')][ind_time][ind_max].to('pc')
        radii = unyt_array(list(prof_df['radius/array_data']), 'pc')

        for field, props in FIELD_DICT.items():
            p_profile = list(prof_df['/'.join([field, 'array_data'])])
            my_fig[gnum].plot(radii.in_units('pc'), p_profile, color= props[0], linestyle= props[1], label= field)
        my_fig[gnum].axvline(x= radius_be)
        my_fig[gnum].annotate(f"time={time:.2f}", xy= (1e1, 5e-18), backgroundcolor= 'white', fontsize= 'medium')

for my_axes in my_fig:
    my_axes.set_xscale('log')
    my_axes.set_xlim(1e-1, 5e2)
    my_axes.set_yscale('log')
    my_axes.set_ylim(1e-18, 1e-10)
    my_axes.grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
for my_axes in my_fig.left_axes:
    my_axes.set_ylabel("Pressure (Pa)")
for i, my_axes in enumerate(my_fig.top_axes):
    star_type = list(DATA_DICT.keys())[i]
    my_axes.xaxis.set_label_position("top")
    my_axes.xaxis.set_label_text(star_type.upper(), fontsize=14)
for my_axes in my_fig.bottom_axes:
    my_axes.set_xlabel("Radius (pc)")
    
my_fig[0].legend(loc= "upper right")
plt.savefig("model_profiles.png")
