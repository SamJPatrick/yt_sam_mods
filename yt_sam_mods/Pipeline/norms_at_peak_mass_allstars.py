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
OUTDIR = "Norms_at_peak_allstars"
FIELDS = ['density', 'temperature', 'metallicity3', 'bonnor_ebert_ratio']
STAR_FILES = ['CCSN_data', 'HN_data', 'PISN_data']
STAR_TYPES = ['ccsn', 'hn', 'pisn']


for field in FIELDS:
    plt.figure()
    label = ' '.join(np.char.capitalize(field.split('_')))
    field_dict = get_field_dict(field)
    plt.title(f"{label}")
    for i, star_type in enumerate(STAR_TYPES):
        time_offset = get_time_offset(star_type)
        df_mass = yt.load(os.path.join(STAR_FILES[i], MASS_FILE))
        index = np.argmax(df_mass.data[('data', 'bonnor_ebert_ratio')][-1])
        times = df_mass.data[('data', 'time')].to('Myr')
        first_index = np.argwhere(times >= time_offset)[0].item()
        data = [df_mass.data[('data', field)][i][index] for i in range (first_index, len(times))]
        plt.plot(times[first_index:] - time_offset, data, label= star_type)
    plt.ylabel(f"{label} ({field_dict['units']})")
    plt.xlabel("Time (Myr)")
    #plt.xlim(0.0, 40.0)
    plt.xlim(0.0, 130.0)
    if (field_dict['log'] == True):
        plt.yscale('log')
    if (field == 'temperature'):
        plt.ylim(1e1, 1e5)
    else :
        plt.ylim(field_dict['limits'])
    if (field == 'bonnor_ebert_ratio'):
        plt.axhline(y= 0.0)
    plt.legend(loc = 'upper right')
    plt.savefig(os.path.join(OUTDIR, f"{field}_at_peak.png"))
    plt.close()
