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



INDIR = "Profiles_pop3/Velocity_profiles_lin"
#INDIR = "Profiles_pop3/Velocity_profiles_log"
OUTDIR = "Profiles_pop3/Velocity_graphs_lin"
#OUTDIR = "Profiles_pop3/Velocity_graphs_log"

FIELDS = ['tangential_velocity_magnitude', 'velocity_spherical_radius']


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


datasets = sorted(glob.glob(os.path.join(INDIR, "DD*_1D_profile_radius_cell_mass.h5")))
for i, dataset in enumerate(datasets):
    dump_name = os.path.basename(dataset)
    ds = yt.load(dataset)
    used = ds.profile.used
    radii = ds.profile.x[used].to('pc')

    for field in FIELDS:
        field_dict = get_field_dict(field)
        plt.figure()
        plt.title(get_title(dump_name, star_type))
        plt.plot(radii, ds.data[('data', field)][used].to('km/s'))
        label = ' '.join(np.char.capitalize(field.split('_')))
        plt.ylabel(f"{label} (km/s)")
        plt.xlabel("Radius (pc)")
        plt.xlim(0, 400)
        #plt.xlim(1e-1, 1e4)
        #plt.xscale('log')
        if (field_dict['log'] == True):
            plt.yscale('log')
        if ("velocity" in field):
            plt.axhline(y= 0.0)
        plt.ylim(field_dict['limits'])
        plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dump_name)}_{field}.png"))
        plt.close()
