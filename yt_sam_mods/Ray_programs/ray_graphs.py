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
DISTANCE_FILE = "ray_distances.txt"
NUM_RAYS = 10
FIELDS = ['density', 'dark_matter_density', 'H2_p0_fraction', 'El_fraction', 'cooling_time', \
          'temperature', 'pressure', 'metallicity3', 'entropy', 'velocity_norm', 'velocity_para']
CUSTOM_LIMS = (1e-29, 1e-21)


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


with open(DISTANCE_FILE, newline='\n') as myfile:
    reader = list(csv.reader(myfile, delimiter= '\t'))
dumps, distances = zip(*[entry for entry in reader])
distances = unyt_array([float(distance) for distance in distances], 'pc')

RHO_BKGRD = unyt_quantity(1e-25, 'g/cm**3')
E_FACTOR = (get_sn_energy(star_type) / RHO_BKGRD)**(1/5)

datasets = glob.glob(os.path.join(INDIR, "DD*_packed.h5"))
datasets = sorted(datasets)
for i, dataset in enumerate(datasets):
    dump_name = os.path.basename(dataset)
    ds = h5py.File(dataset, 'r')
    if (get_time_z(dump_name, star_type)[0] > get_lifetime_offset(star_type)):
        grad_field = "temperature"
    else :
        grad_field = "El_fraction"
    means = transpose_unyt([np.mean(transpose_unyt(x)) for x in \
                            zip(*[unyt_array(ds[f'{grad_field}_{n}/array_data'], 'g/cm**3') for n in range (NUM_RAYS)])])
    grads = [((means[i+1] - means[i]) / means[i]) for i in range (len(means) - 1)]
    index_max = np.argmin(transpose_unyt(grads)) + 1
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
        if (get_time_z(dump_name, star_type)[0] > get_lifetime_offset(star_type)):
            plt.axvline(x= dist_theo, color='red', linestyle='--')
            plt.axvline(x= ds[f"distances/array_data"][index_max], color='red')
        else :
            plt.axvline(x= ds[f"distances/array_data"][index_max], color='green')
        plt.axvline(x= distances[i], color='blue')
        plt.ylim(field_dict['limits'])
        plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dump_name)}_{field}.png"))
        plt.close()


'''    
for dataset in datasets:
    df = h5py.File(dataset, 'r')
    dump_name = os.path.basename(dataset)
    
    #distances = unyt_array(df[f'distances_{n}/array_data'], 'pc')
    distances = unyt_array(df['distances'][:], 'pc')
    #velocities = unyt_array(df[f'velocity_{n}/array_data'], 'km/s')
    velocities = transpose_unyt([np.mean(transpose_unyt(vels)) for vels in \
                                 zip(*[unyt_array(df[f'velocity_{n}'][:], 'km/s') for n in range (NUM_RAYS)])])
    #temperatures = unyt_array(df[f'temperature_{n}/array_data'], 'K')
    temperatures = transpose_unyt([np.mean(transpose_unyt(temps)) for temps in \
                                   zip(*[unyt_array(df[f'temperature_{n}'][:], 'K') for n in range (NUM_RAYS)])])
    #densities = unyt_array(df[f'distances_{n}/array_data'], 'pc')
    densities = transpose_unyt([np.mean(transpose_unyt(dens)) for dens in \
                                zip(*[unyt_array(df[f'density_{n}'][:], 'g/cm**3') for n in range (NUM_RAYS)])])
        
    velocities_diff = transpose_unyt([(velocities[i+1] - velocities[i]) / (distances[i+1] - distances[i]) \
                                    for i in range (len(velocities) - 1)])
    temperatures_diff = transpose_unyt([(temperatures[i+1] - temperatures[i]) / (distances[i+1] - distances[i]) \
                                    for i in range (len(temperatures) - 1)])

    plt.figure()
    plt.title(get_title(dump_name, star_type))
    accl_vel = np.abs(velocities[1:] * velocities_diff).to('m*s**(-2)')
    accl_temp = np.abs((GAMMA/(GAMMA-1))*(kb/(MU*mh)) * temperatures_diff).to('m*s**(-2)')
    plt.plot(distances[1:], accl_vel, label= 'velocity')
    plt.plot(distances[1:], accl_temp, label= 'temperature')
    plt.ylabel("$\dot{v} (m s^{-2})$")
    plt.yscale('log')
    plt.ylim(1e-15, 1e3)
    plt.xlabel("Radius (pc)")
    plt.xlim(0.0, 400.0)
    plt.legend(loc= 'upper right')
    plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dump_name)}_accelerations.png"))
    plt.close()

    
    plt.figure()
    plt.title(get_title(dump_name, star_type))
    energy_vel = (densities * velocities**2).to('erg*m**(-3)')
    energy_temp = ((1/(GAMMA-1))*(kb/(MU*mh)) * temperatures).to('erg*m**(-3)')
    plt.plot(distances, energy_vel, energy_temp)
    plt.ylabel("$E (erg)")
    plt.yscale('log')
    plt.xlabel("Radius (pc)")
    plt.xlim(0.0, 400.0)
    plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dump_name)}_energies.png"))
    plt.close()
    '''

