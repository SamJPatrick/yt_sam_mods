import os
import sys
import csv
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
from unyt import Myr, yr, pc, kb, mh


GAMMA = 5/3
MU = 1.6

INDIR_RAYS = "Sobolev/Ray_profiles"
INDIR_VELS = "Profiles_pop3/Velocity_profiles_lin"
INDIR_PROF = "Profiles_pop3/Normal_profiles_lin"
OUTDIR = "Sobolev/Ray_graphs_linefill"
DISTANCE_FILE = "ray_distances.txt"

NUM_RAYS = 20
X_LIM = 400.0
FIELDS = ['density', 'El_fraction', 'temperature', 'velocity_para']
DATASET_NUMS = [100, 110]
#DATASET_NUMS = [120, 130, 140, 150, 160]
COLORS = ['blue', 'orange', 'magenta', 'cyan', 'brown', 'red']


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


with open(DISTANCE_FILE, newline='\n') as myfile:
    reader = list(csv.reader(myfile, delimiter= '\t'))
dumps, distances = zip(*[entry for entry in reader])
distances = unyt_array([float(distance) for distance in distances], 'pc')
dumps = np.array(dumps)

#RHO_BKGRD = unyt_quantity(1e-25, 'g/cm**3')
#E_FACTOR = (get_sn_energy(star_type) / RHO_BKGRD)**(1/5)
TIME_OFFSET = unyt_quantity(0.15, 'Myr')

data_rays = sorted([os.path.join(INDIR_RAYS, f"DD{num:04d}_packed.h5") for num in DATASET_NUMS])
data_prof = sorted([os.path.join(INDIR_PROF, f"DD{num:04d}_profile_weight_field_cell_mass.h5") for num in DATASET_NUMS])
data_vels = sorted([os.path.join(INDIR_VELS, f"DD{num:04d}_1D_profile_radius_cell_mass.h5") for num in DATASET_NUMS])
assert len(data_prof) == len(data_vels), "Error, unequal number of datasets found"
assert len(data_vels) == len(data_rays), "Error, unequal number of datasets found"

for field in FIELDS:
    field_dict = get_field_dict(field)
    plt.figure()
    plt.title(field.replace('_', ' ').capitalize())
    for i in range (len(DATASET_NUMS)):
        try :
            ds_rays = h5py.File(data_rays[i], 'r')
            ds_norm = yt.load(data_prof[i])
            ds_vels = yt.load(data_vels[i])
        except Exception:
            print(f"Error, could not load dataset {DATASET_NUMS[i]}, continuing")
            continue

        dump_name = os.path.basename(data_rays[i])
        time = get_time_z(dump_name, star_type)[0].to('Myr')
        star_lifetime = get_lifetime_offset(star_type)
        if (time > TIME_OFFSET):
            if (time > star_lifetime):
                grad_field = "temperature"
            else :
                grad_field = "El_fraction"
            means = [np.mean(x) for x in zip(*[ds_rays[f'{grad_field}_{n}/array_data'] for n in range (NUM_RAYS)])]
            grads = [((means[i+1] - means[i]) / means[i]) for i in range (len(means) - 1)]
            indices = np.argsort(grads) 
            #dist_theo = (E_FACTOR * (get_time_z(dump_name, star_type)[0] - get_lifetime_offset(star_type))**(2/5)).to('pc')

        rays_max = transpose_unyt([np.nanmax(transpose_unyt(x)) for x in \
                              zip(*[unyt_array(ds_rays[f'{field}_{n}/array_data'], field_dict['units']) for n in range (NUM_RAYS)])])
        rays_min = transpose_unyt([np.nanmin(transpose_unyt(x)) for x in \
                              zip(*[unyt_array(ds_rays[f'{field}_{n}/array_data'], field_dict['units']) for n in range (NUM_RAYS)])])
        arr_ray = transpose_unyt([np.nanmean(transpose_unyt(x)) for x in \
                              zip(*[unyt_array(ds_rays[f'{field}_{n}/array_data'], field_dict['units']) for n in range (NUM_RAYS)])])
        if (field == 'velocity_rad'):
            used = ds_vels.profile.used
            radii = ds_vels.profile.x[used]
            arr_ray = arr_ray.to('km/s')
            arr_prof = ds_vels.profile[('data', 'tangential_velocity_magnitude')][used].to('km/s')
        elif (field == 'velocity_para'):
            used = ds_vels.profile.used
            radii = ds_vels.profile.x[used]
            arr_ray = arr_ray.to('km/s')
            arr_prof = ds_vels.profile[('data', 'velocity_spherical_radius')][used].to('km/s')
        else :
            used = ds_norm.profile.used
            radii = ds_norm.profile.x[used]
            arr_prof = ds_norm.profile[('data', field)][used].to(field_dict['units'])
        plt.plot(ds_rays[f"distances/array_data"], arr_ray, color= COLORS[i], label= f'{time:.2f}')
        plt.plot(radii, arr_prof, color= COLORS[i], linestyle=':')
        if ("velocity" in field):
            plt.fill_between(ds_rays[f"distances/array_data"], rays_min.to('km/s'), rays_max.to('km/s'), alpha=0.2)
        else :
            plt.fill_between(ds_rays[f"distances/array_data"], rays_min, rays_max, alpha=0.2)

        if (time > TIME_OFFSET):
            if (i == 4):
                plt.axvline(x= ds_rays[f"distances/array_data"][indices[4]+1], color=COLORS[i], linestyle='--')
            else :
                plt.axvline(x= ds_rays[f"distances/array_data"][indices[0]+1], color=COLORS[i], linestyle='--')
        #plt.axvline(x= distances[np.argwhere(f"DD{get_dump_num(dump_name)}" == dumps).item()], color=COLORS[i], linestyle='-')
        x_halo = distances[np.argwhere(f"DD{get_dump_num(dump_name)}" == dumps).item()] / unyt_quantity(X_LIM, 'pc')
        plt.annotate('', xy= (x_halo, 0.0), xycoords='axes fraction', xytext= (x_halo + 0.05, -0.1), \
                     textcoords= 'axes fraction', arrowprops={'facecolor': COLORS[i], 'shrink': 0.05})
        
    if (field_dict['log'] == True):
        plt.yscale('log')
    if ('velocity' in field):
        plt.axhline(y= 0.0)
    label = ' '.join(np.char.capitalize(field.split('_')))
    plt.ylabel(f"{label} ({field_dict['units']})")
    plt.xlabel("Radius (pc)")
    plt.xlim(0.0, X_LIM)
    plt.ylim(field_dict['limits'])
    plt.legend(loc='upper right')
    plt.savefig(os.path.join(OUTDIR, f"{field}_linefill_hii.png"))
    plt.close()
