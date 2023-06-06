import os
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
from unyt import Myr, yr, pc, kb, mh


GAMMA = 5/3
MU = 1.6

#INDIR = "Sobolev/Ray_profiles"
INDIR = "Sobolev/Ray_profiles/Packed_rays"
OUTDIR = "Sobolev/Ray_graphs"
NUM_RAYS = 10
FIELDS = ['density', 'temperature', 'pressure', 'entropy', 'metallicity3', 'H2_p0_fraction', 'velocity']
CUSTOM_LIMS = (1e-29, 1e-21)


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass

#datasets = glob.glob(os.path.join(INDIR, "DD*_sobolev.h5"))
datasets = glob.glob(os.path.join(INDIR, "DD*_packed.h5"))
datasets = sorted(datasets)
#datasets.sort(key = lambda x: get_time_z(os.path.basename(x), star_type)[0])

'''
for dataset in datasets:
    dump_name = os.path.basename(dataset)
    ds = h5py.File(dataset, 'r')
    for field in FIELDS:    
        field_dict = get_field_dict(field)
        plt.figure()
        plt.title(get_title(dump_name, star_type))
        for n in range (NUM_RAYS):
            #plt.plot(ds[f"distances_{n}/array_data"], ds[f"{field}_{n}/array_data"])
            plt.plot(ds["distances"][:], ds[f"{field}_{n}"][:])
        label = ' '.join(np.char.capitalize(field.split('_')))
        plt.ylabel(f"{label} ({field_dict['units']})")
        plt.xlabel("Radius (pc)")
        plt.xlim(0.0, 400.0)
        if (field_dict['log'] == True):
            plt.yscale('log')
        if (field == 'temperature'):
            plt.axhline(y= 0.0)
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
    plt.xlabel("Radius (pc)")
    plt.xlim(0.0, 400.0)
    plt.legend(loc= 'upper right')
    plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dump_name)}_accelerations.png"))
    plt.close()

    '''
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

