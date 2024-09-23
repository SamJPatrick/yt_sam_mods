import os
import yt
import sys
import numpy as np
from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from grid_figure import GridFigure

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc


PRE_DIR = "Britton_sim_data"
MASS_FILE = "star_None_mass.h5"
SIMULATION_FILE = "simulation.h5"
OUTDIR = "."

FIELDS = ['accretion_rate', 'accretion_rate_z']
COLOR_DICT = {'ccsn': 'blue', 'hn': 'orange', 'pisn': 'green'}





my_fig = GridFigure(3, 1, figsize=(13, 9),
                    left_buffer=0.08, right_buffer=0.02,
                    bottom_buffer=0.06, top_buffer=0.06,
                    vertical_buffer=0.04, horizontal_buffer=0.00)


for i, my_axes in enumerate(my_fig):
    
    star_type = list(COLOR_DICT.keys())[i]
    df_norm = yt.load(os.path.join(PRE_DIR, star_type.upper(), "star_None_mass_normal_profiles.h5"))
    df_vel = yt.load(os.path.join(PRE_DIR, star_type.upper(), "star_None_mass_velocity_profiles_test.h5"))
    ind_max_be = np.argmax(df_norm.data[('data', 'bonnor_ebert_ratio')][-1])
    ind_zero_be = ind_max_be #+ np.argmin(np.abs(df_norm.data[('data', 'bonnor_ebert_ratio')][-1][ind_max_be:] - 1.0))
    times, t_indicies = get_xaxis_data(star_type, df_norm)
    acc_rate = 4 * np.pi * transpose_unyt([-df_vel.data[('data', 'velocity_spherical_radius')][t_index][ind_max_be] * \
                                           df_norm.data[('data', 'density')][t_index][ind_max_be] * \
                                           df_norm.data[('data', 'radius')][t_index][ind_max_be]**2 for t_index in t_indicies]).to('Msun/yr')
    my_axes.plot(times, acc_rate, color= 'blue', label= "inward")
    my_axes.plot(times, -acc_rate, color= 'orange', label= "outward")
    my_axes.yaxis.set_label_text(f"Accretion rate {star_type.upper()} " + r"(M$_{\odot}$ yr$^{-1}$)")

    if (star_type == 'pisn'):
        my_axes.set_xlim(0.0, 130.0)
        my_axes.xaxis.set_ticks(np.arange(0, 131, 10))
    else :
        my_axes.set_xlim(0.0, 30.0)
        my_axes.xaxis.set_ticks(np.arange(0, 31, 2))
    my_axes.set_yscale('log')
    my_axes.set_ylim(1e-5, 1e0)
    my_axes.grid(visible=True, axis="both", zorder=1, linestyle=":", color="black", alpha=0.6)
    my_axes.tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)
    
my_fig[0].legend(loc= 'upper right')
my_fig[-1].xaxis.set_label_text("Time (Myr)")
plt.savefig("model_profiles.png")
