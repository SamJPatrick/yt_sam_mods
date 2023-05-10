import os
import yt
import sys
import glob
import numpy as np
from  yt.extensions.sam_mods.graph_funcs import *

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc



RADIUS_FILE = "star_None_radius.h5"
PROFILE_PATH = "Profiles/Normal_profiles"
SIMULATION_FILE = "simulation.h5"
OUTDIR = "Profiles/Test_graphs"
FIELDS = ['density', 'temperature', 'velocity_magnitude', 'sound_speed']


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass
time_offset = get_time_offset(star_type)


files = sorted(glob.glob(os.path.join(PROFILE_PATH, "DD*_profile_weight_field_cell_mass.h5")))
for df_path in files:
    for field in FIELDS:
        field_dict = get_field_dict(field)
        df = yt.load(df_path)
        df_name = os.path.basename(df_path)
        time = get_time_z(df_name, star_type)[0]
        plt.figure()
        plt.title(get_title(df_name, star_type))
        plt.plot(df.data[('data', 'radius')], df.data[('data', f'{field}')])
        label = ' '.join(np.char.capitalize(field.split('_')))
        plt.ylabel(f"{label} ({field_dict['units']})")
        plt.xlabel("Radius (pc)")
        plt.xlim(1e-2, 5e2)
        plt.xscale('log')
        if (field_dict['log'] == True):
            plt.yscale('log')
            plt.ylim(field_dict['limits'])
        plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(df_name)}_{field}.png"))
        plt.close()
