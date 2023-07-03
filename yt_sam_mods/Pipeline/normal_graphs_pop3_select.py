import os
import yt
import sys
import glob
import numpy as np
from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc



INDIR = "Profiles_pop3/Normal_profiles"
OUTDIR = "Profiles_pop3"
FIELDS = ['density', 'temperature']
DATASET_NUMS = [98, 100, 110]   # <----- MODIFY SELECTED DATA DUMPS HERE!!


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


datasets = [os.path.join(INDIR, f"DD{num:04d}_profile_weight_field_cell_mass.h5") for num in DATASET_NUMS]
for field in FIELDS:
    plt.figure()
    plt.title(star_type.upper())        
    for j, dataset in enumerate(datasets):
        dump_name = os.path.basename(dataset)
        ds = yt.load(dataset)
        field_dict = get_field_dict(field)
        used = ds.profile.used
        radii = ds.profile.x[used].to('pc')
        plt.plot(radii, ds.data[('data', field)][used], \
                 label= f"{get_time_z(os.path.basename(dataset), star_type)[0]:.2f}")               
    label = ' '.join(np.char.capitalize(field.split('_')))
    plt.ylabel(f"{label} ({field_dict['units']})")
    plt.xlabel("Radius (pc)")
    plt.xlim(1e-1, 1e4)
    plt.xscale('log')
    plt.legend()
    if (field_dict['log'] == True):
        plt.yscale('log')
    plt.ylim(field_dict['limits'])
    plt.savefig(os.path.join(OUTDIR, f"combined_times_{field}.png"))
    plt.close()
