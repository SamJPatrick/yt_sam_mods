import os
import yt
import sys
import numpy as np
from  yt.extensions.sam_mods.graph_funcs import get_time_offset

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc


MASS_FILE = "star_None_mass.h5"
SIMULATION_FILE = "simulation.h5"


try :
    star_mode = sys.argv[1]
except IndexError:
    star_mode = ""
    pass
time_offset = get_time_offset(star_mode)

fields = ["turbulent_sound_crossing_time", "sound_crossing_time", "total_dynamical_time", "cooling_time"]
labels = ["turbulent sound-crossing", "sound-crossing", "free-fall", "cooling"]
colors = ["orange", "red", "green", "blue"]
assert len(fields) == len(labels) == len(colors), "Error, field lengths do not match"

df_mass = yt.load(MASS_FILE)
index = np.argmax(df_mass.data[('data', 'bonnor_ebert_ratio')][-1])
times = df_mass.data[('data', 'time')].to('Myr')
timescale_dict = {field: unyt_array([0.0 for i in range (len(times))], 'Myr') for field in fields}
for field in fields:
    for i, time in enumerate(times):
        timescale_dict[field][i] = df_mass.data[('data', field)][i][index].to('Myr')

plt.figure()
plt.title("Timescales at BE peak vs time")
for field, color, label in zip(fields, colors, labels):
    plt.plot(times - time_offset, timescale_dict[field], color=color, alpha=0.7, linewidth=1.5, label=label)
plt.ylabel("Timescales (Myr)")
plt.xlabel("Time (Myr)")
#plt.xlim(0.0, 40.0)
plt.xlim(0.0, 130.0)
plt.yscale('log')
plt.ylim(1e0, 1e6)
plt.legend(loc='upper right')
plt.savefig("timescales_at_peak.png")
plt.close()
