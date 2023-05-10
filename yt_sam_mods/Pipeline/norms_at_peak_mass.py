import os
import yt
import sys
import numpy as np
from  yt.extensions.sam_mods.graph_funcs import *

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc



MASS_FILE = "star_None_mass.h5"
SIMULATION_FILE = "simulation.h5"
OUTDIR = "Profiles/Norms_at_peak_mass"
FIELDS = ['density', 'temperature', 'velocity_magnitude', 'sound_speed',
          'accretion_rate', 'accretion_rate_z', 'metallicity3', 'bonnor_ebert_ratio']


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass
time_offset = get_time_offset(star_type)

df_mass = yt.load(MASS_FILE)
index = np.argmax(df_mass.data[('data', 'bonnor_ebert_ratio')][-1])
times = df_mass.data[('data', 'time')].to('Myr')
norms_dict = {field: unyt_array([0.0 for i in range (len(times))], get_field_dict(field)['units']) for field in FIELDS}
for field in FIELDS:
    for i, time in enumerate(times):
        norms_dict[field][i] = np.abs(df_mass.data[('data', field)][i][index])
first_index = np.argwhere(times >= time_offset)[0].item()

for field, arr in norms_dict.items():
    field_dict = get_field_dict(field)
    label = ' '.join(np.char.capitalize(field.split('_')))
    plt.figure()
    plt.title(f"{label} at BE Peak")
    plt.plot(times[first_index:] - time_offset, norms_dict[field][first_index:])
    plt.ylabel(f"{label} ({field_dict['units']})")
    plt.xlabel("Time (Myr)")
    #plt.xlim(0.0, 40.0)
    plt.xlim(0.0, 130.0)
    if (field_dict['log'] == True):
        plt.yscale('log')
    plt.ylim(field_dict['limits'])
    plt.savefig(os.path.join(OUTDIR, f"{field}_at_peak.png"))
    plt.close()
