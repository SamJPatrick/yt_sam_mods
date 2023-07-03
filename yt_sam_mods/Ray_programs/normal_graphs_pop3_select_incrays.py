import os
import yt
import sys
import glob
import h5py
import numpy as np
from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc



NUM_RAYS = 10
INDIR_PROFS = "Profiles_pop3/Normal_profiles"
INDIR_VELS = "Profiles_pop3/Velocity_profiles"
INDIR_RAYS = "Sobolev_old/Ray_profiles/Packed_rays"
OUTDIR = "Profiles_pop3"
FIELDS = ['density', 'temperature', 'H2_p0_fraction', 'velocity_spherical_radius']
DATASET_NUMS = [120, 130, 140, 150, 160, 165]   # <----- MODIFY SELECTED DATA DUMPS HERE!!
#DATASET_NUMS = [98, 100, 110]
COLORS = ['b', 'g', 'r', 'c', 'm', 'y']


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


data_profs = sorted([os.path.join(INDIR_PROFS, f"DD{num:04d}_profile_weight_field_cell_mass.h5") for num in DATASET_NUMS])
data_vels = sorted([os.path.join(INDIR_VELS, f"DD{num:04d}_1D_profile_radius_cell_mass.h5") for num in DATASET_NUMS])
data_rays = sorted([os.path.join(INDIR_RAYS, f"DD{num:04d}_packed.h5") for num in DATASET_NUMS])
assert len(data_profs) == len(data_vels), "Error, unequal number of datasets found"
assert len(data_vels) == len(data_rays), "Error, unequal number of datasets found"


for field in FIELDS:
    plt.figure()
    plt.title(star_type.upper())
    field_dict = get_field_dict(field)
    for i in range (len(data_profs)):
        if ("velocity" in field):
            ds_prof = yt.load(data_vels[i])
        else :
            ds_prof = yt.load(data_profs[i])
        used = ds_prof.profile.used
        radii = ds_prof.profile.x[used].to('pc')
        ds_ray = h5py.File(data_rays[i], 'r')
        distances = unyt_array(ds_ray['distances'][:], 'pc')
        arr_prof = ds_prof.data[('data', field)][used]
        if ("velocity" in field):
            arr_ray = transpose_unyt([np.mean(transpose_unyt(x)) for x in \
                                      zip(*[unyt_array(ds_ray[f'velocity_{n}'][:], field_dict['units']) for n in range (NUM_RAYS)])])
        else :
            arr_ray = transpose_unyt([np.mean(transpose_unyt(x)) for x in \
                                      zip(*[unyt_array(ds_ray[f'{field}_{n}'][:], field_dict['units']) for n in range (NUM_RAYS)])])
        plt.plot(radii, arr_prof, color= COLORS[i], linestyle= ':')
        plt.plot(distances, arr_ray, label= f"{get_time_z(os.path.basename(data_rays[i]), star_type)[0]:.2f}", color= COLORS[i], linestyle= '-')
    label = ' '.join(np.char.capitalize(field.split('_')))
    plt.ylabel(f"{label} ({field_dict['units']})")
    plt.xlabel("Radius (pc)")
    #plt.xlim(1e-1, 1e4)
    #plt.xscale('log')
    plt.xlim(0.0, 400.0)
    plt.legend(loc='upper right', fontsize='x-small')
    if ("velocity" in field):
        plt.axhline(y= 0.0)
    if (field_dict['log'] == True):
        plt.yscale('log')
    plt.ylim(field_dict['limits'])
    plt.savefig(os.path.join(OUTDIR, f"combined_times_{field}_sedov.png"))
    plt.close()
