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
MASS_FILE = "star_None_mass_velocity_profiles.h5"
STAR_DIRS = ['CCSN', 'HN', 'PISN']

my_fig = GridFigure(1, 3, figsize=(13, 5),
                    left_buffer=0.10, right_buffer=0.06,
                    bottom_buffer=0.12, top_buffer=0.06,
                    vertical_buffer=0.00, horizontal_buffer=0.08)


for i, my_axes in enumerate(my_fig):
    
    my_axes.set_xlabel("Radius (pc)", labelpad=8)
    my_axes.set_ylabel(f"Radial velocity {STAR_DIRS[i]}" + r" (km s$^{-1}$)", labelpad=8)
    my_axes.set_xlim(1e-3, 5e2)
    my_axes.set_xscale('log')
    #my_axes.xaxis.set_ticks(np.logspace(5e-3, 5e2, 10))
    my_axes.grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
    my_axes.tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)        
    if (i == 2):
        my_axes.set_ylim(-0.5, 12.0)
    else :
        my_axes.set_ylim(-1.0, 3.0)

    df_vel = yt.load(os.path.join(PRE_DIR, STAR_DIRS[i], MASS_FILE))
    used = df_vel.data[('data', 'used')][-1].astype(bool)
    velocities = - df_vel.data[('data', 'velocity_spherical_radius')][-1][used].to('km/s')
    radii = df_vel.data[('data', 'radius')][-1][used].to('pc')
    my_axes.plot(radii, velocities)
    my_axes.axhline(y= 0.0, color= 'red')

plt.savefig("model_profiles.png")
