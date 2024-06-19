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



RADIUS_FILE = "star_None_radius.h5"
SIMULATION_FILE = "simulation.h5"
OUTDIR = "Profiles/Normal_graphs"
#FIELDS = ['density', 'temperature', 'velocity_magnitude', 'sound_speed']
FIELDS = ['metallicity3']

try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass
time_offset = get_time_offset(star_type)


df_sim = yt.load(SIMULATION_FILE)
df_radius = yt.load(RADIUS_FILE)


times = df_radius.data[('data', 'time')].to('Myr')
radii = df_radius.data[('data', 'radius')][1:].to('pc')
index_be = np.argmax(df_radius.data[('data', 'bonnor_ebert_ratio')][-1])
radius_be = radii[index_be]
for field in FIELDS:
    field_dict = get_field_dict(field)
    for i, time in enumerate(times):
        index_sim = np.argwhere(df_sim.data[('data', 'time')].to('Myr') == time).item()
        dsfn = df_sim.data[('data', 'filename')].astype(str)[index_sim].split('/')[-1]
        plt.figure()
        plt.title(get_title(dsfn, star_type))
        plt.plot(radii, df_radius.data[('data', field)][i])
        label = ' '.join(np.char.capitalize(field.split('_')))
        plt.ylabel(f"{label} ({field_dict['units']})")
        plt.xlabel("Radius (pc)")
        plt.xlim(1e-2, 5e2)
        plt.xscale('log')
        if (field_dict['log'] == True):
            plt.yscale('log')
        plt.ylim(field_dict['limits'])
        plt.axvline(x= radius_be)
        plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dsfn)}_{field}.png"))
        plt.close()
