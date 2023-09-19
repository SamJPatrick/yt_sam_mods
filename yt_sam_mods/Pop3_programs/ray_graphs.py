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

INDIR = "Sobolev/Ray_profiles"
OUTDIR = "Sobolev/Ray_graphs"
#DISTANCE_FILE = "ray_distances.txt"

NUM_RAYS = 20
FIELDS = ['density', 'dark_matter_density', 'H2_p0_fraction', 'El_fraction', 'cooling_time', \
          'temperature', 'pressure', 'metallicity3', 'entropy', 'velocity_para', 'velocity_rad', 'velocity_phi']



try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


#with open(DISTANCE_FILE, newline='\n') as myfile:
#    reader = list(csv.reader(myfile, delimiter= '\t'))
#dumps, distances = zip(*[entry for entry in reader])
#distances = unyt_array([float(distance) for distance in distances], 'pc')

RHO_BKGRD = unyt_quantity(1e-25, 'g/cm**3')
E_FACTOR = (get_sn_energy(star_type) / RHO_BKGRD)**(1/5)
TIME_OFFSET = unyt_quantity(0.15, 'Myr')

datasets = sorted(glob.glob(os.path.join(INDIR, "DD*_packed.h5")))
for i, dataset in enumerate(datasets):
    
    dump_name = os.path.basename(dataset)
    try :
        ds = h5py.File(dataset, 'r')
    except Exception:
        continue
    time = get_time_z(dump_name, star_type)[0].to('Myr')
    star_lifetime = get_lifetime_offset(star_type)
    if (time > TIME_OFFSET):
        if (time > star_lifetime):
            grad_field = "temperature"
        else :
            grad_field = "El_fraction"
        means = [np.mean(x) for x in zip(*[ds[f'{grad_field}_{n}/array_data'] for n in range (NUM_RAYS)])]
        grads = [((means[i+1] - means[i]) / means[i]) for i in range (len(means) - 1)]
        index_max = np.argmin(grads) + 1
        dist_theo = (E_FACTOR * (get_time_z(dump_name, star_type)[0] - get_lifetime_offset(star_type))**(2/5)).to('pc')
    
    for field in FIELDS:    
        field_dict = get_field_dict(field)
        plt.figure()
        plt.title(get_title(dump_name, star_type))
        for n in range (NUM_RAYS):
            plt.plot(ds[f"distances/array_data"], ds[f"{field}_{n}/array_data"])
        label = ' '.join(np.char.capitalize(field.split('_')))
        plt.ylabel(f"{label} ({field_dict['units']})")
        plt.xlabel("Radius (pc)")
        plt.xlim(0.0, 400.0)
        if (field_dict['log'] == True):
            plt.yscale('log')
        if ('velocity' in field):
            plt.axhline(y= 0.0)
        if (time > TIME_OFFSET):
            if (time > star_lifetime):
                plt.axvline(x= dist_theo, color='red', linestyle='--')
                plt.axvline(x= ds[f"distances/array_data"][index_max], color='red')
            else :
                plt.axvline(x= ds[f"distances/array_data"][index_max], color='green')
        #plt.axvline(x= distances[i], color='blue')
        plt.ylim(field_dict['limits'])
        plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dump_name)}_{field}.png"))
        plt.close()
