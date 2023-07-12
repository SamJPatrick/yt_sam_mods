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



INDIR = "Profiles_pop3/Normal_profiles_lin"
#INDIR = "Profiles_pop3/Normal_profiles_log"
OUTDIR = "Profiles_pop3/Normal_graphs_lin"
#OUTDIR = "Profiles_pop3/Normal_graphs_log"

FIELDS = ['density', 'dark_matter_density', 'metallicity3', 'H2_p0_fraction', 'El_fraction', \
          'temperature', 'pressure', 'entropy', 'sound_speed', 'velocity_magnitude']


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


datasets = sorted(glob.glob(os.path.join(INDIR, "DD*_profile_weight_field_cell_mass.h5")))
for i, dataset in enumerate(datasets):
    dump_name = os.path.basename(dataset)
    ds = yt.load(dataset)
    used = ds.profile.used
    radii = ds.profile.x[used].to('pc')

    for field in FIELDS:
        field_dict = get_field_dict(field)
        plt.figure()
        plt.title(get_title(dump_name, star_type))
        plt.plot(radii, ds.data[('data', field)][used].to(field_dict['units']))
        label = ' '.join(np.char.capitalize(field.split('_')))
        plt.ylabel(f"{label} ({field_dict['units']})")
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
