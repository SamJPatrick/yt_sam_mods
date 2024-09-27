import os
import yt
import sys
import numpy as np

from yt.extensions.sam_mods.misc import *
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.graph_funcs import *

from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc
from grid_figure import GridFigure

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colorbar import ColorbarBase

import pdb



RADIUS_FILE = "star_None_radius_normal_profiles.h5"
PRE_DIR = "Britton_sim_data"

X_LIMS = unyt_array([1e-1, 2e2], 'pc')
Y_LIMS = unyt_array([1e0, 1e4], 'Myr')
DELTA_Y = 1.5

DATA_DIRS = ['DD0254', 'DD0258', 'DD0273', 'DD0278', 'DD0284', 'DD0295']
fields = ["turbulent_sound_crossing_time", "sound_crossing_time", "total_dynamical_time", "cooling_time"]
labels = ["turbulent sound-crossing", "sound-crossing", "free-fall", "cooling"]
colors = ["orange", "red", "green", "blue"]
assert len(fields) == len(labels) == len(colors), "Error, field length do not match"



try :
    star_type = sys.argv[1]
except IndexError:
    star_type = "pisn"
    pass
time_offset = get_time_offset(star_type)


# blue = 0 = frag;  green = 1 = free-fall; yellow = 2 = turb-support; red = 3 = therm-support
def get_state(profile_dict, n):
    if (profile_dict['cooling_time'][n] < profile_dict['total_dynamical_time'][n]):
        bkgcolor = 'blue'
        state = 0
    elif (profile_dict['total_dynamical_time'][n] < profile_dict['turbulent_sound_crossing_time'][n]
          and profile_dict['total_dynamical_time'][n] < profile_dict['sound_crossing_time'][n]):
        bkgcolor = 'green'
        state = 1
    elif (profile_dict['total_dynamical_time'][n] < profile_dict['sound_crossing_time'][n]):
        bkgcolor = 'yellow'
        state = 2
    else :
        bkgcolor = 'red'
        state = 3
    return bkgcolor



df_sim = yt.load(os.path.join(PRE_DIR, star_type.upper(), "simulation.h5"))
df_radius = yt.load(os.path.join(PRE_DIR, star_type.upper(), RADIUS_FILE))
idx_max = np.argmax(df_radius.data[('data', 'bonnor_ebert_ratio')][-1])
my_fig = GridFigure(2, 3, figsize=(15, 9),
                    top_buffer= 0.02, bottom_buffer= 0.06,
                    left_buffer= 0.06, right_buffer= 0.02,
                    horizontal_buffer= 0.08, vertical_buffer= 0.04)

for i, my_axes in enumerate(my_fig):
    
    time = get_time_z(DATA_DIRS[i], star_type,
                      sim_file= os.path.join(PRE_DIR, star_type.upper(), "simulation.h5"))[0]
    index = np.argmin(np.abs(df_radius.data[('data', 'time')] - get_time_offset(star_type) - time))
    radii = df_radius.data[('data', 'radius')][1:].to('pc')
    radius_be = radii[idx_max].to('pc')
    
    profile_dict = {}
    for field, label, color in zip(fields, labels, colors):        
        profile_dict[field] = df_radius.data[('data', field)][index].to('Myr')
        my_axes.plot(radii, df_radius.data[('data', field)][index].to('Myr'), color=color,
                     alpha=0.7, linewidth=1.5, label=label)

    states = [get_state(profile_dict, i) for i in range (len(radii))]
    for j in range (len(radii) - 1):
        my_axes.fill_between(radii, Y_LIMS[0] * DELTA_Y,
                             where=((radii >= radii[j]) & (radii <= radii[j+1])), color= get_state(profile_dict, j))
    
    my_axes.set_xscale('log')
    my_axes.set_yscale('log')
    my_axes.set_xlim(*X_LIMS)
    my_axes.set_ylim(*Y_LIMS)                
    my_axes.yaxis.set_label_text(f"Timescales (Myr) @ t= {time:.1f}")
    my_axes.axvline(x= radius_be, color='black')
    #my_axes.xaxis.set_ticks(np.logspace(5e-3, 5e2, 10))
    my_axes.grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
    my_axes.tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)
    if (i > 2):
        my_axes.xaxis.set_label_text(f"Radius (pc)")

my_fig[0].legend(loc='upper left')
plt.savefig("model_profiles.png")
plt.close()
