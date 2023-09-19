import glob
import h5py
import numpy as np
import sys
import os
import csv

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from unyt import unyt_quantity, unyt_array
from unyt import kb, mh
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from yt.extensions.sam_mods.graph_funcs import *




NUM_RAYS = 20

INDIR_RAYS = "Sobolev/Ray_profiles"
OUTDIR_NU = "Sobolev/Viscous_params"

MU = 1.6
GAMMA = 5/3
CONV_CONST = kb / (MU * mh)


LENGTHS = unyt_array([5.0, 10.0, 15.0, 20.0, 25.0, 30.0], 'pc')
COLORS = ['red', 'orange', 'yellow', 'green', 'blue', 'purple']


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


ray_files = sorted(glob.glob(os.path.join(INDIR_RAYS, "DD*_packed.h5")))
for i, ray_file in enumerate(ray_files):
    
    df = h5py.File(ray_file, 'r')
    dump_name = os.path.basename(ray_file)
    distances = unyt_array(df['distances/array_data'], 'pc')
    arr_len = len(distances)
    time = get_time_z(dump_name, star_type)[0].to('Myr')
    
    v_para = [transpose_unyt(x) for x in zip(*[unyt_array(df[f'velocity_para_{n}/array_data'], 'km/s') for n in range (NUM_RAYS)])]
    v_rad = [transpose_unyt(x) for x in zip(*[unyt_array(df[f'velocity_rad_{n}/array_data'], 'km/s') for n in range (NUM_RAYS)])]
    v_phi = [transpose_unyt(x) for x in zip(*[unyt_array(df[f'velocity_phi_{n}/array_data'], 'km/s') for n in range (NUM_RAYS)])]
    temps = [transpose_unyt(x) for x in zip(*[unyt_array(df[f'temperature_{n}/array_data'], 'K') for n in range (NUM_RAYS)])]
    del_phi = np.ptp(df['angle_list/array_data'])

    v_para_ave = transpose_unyt([np.nanmean(vels) for vels in v_para])
    v_rad_ave = transpose_unyt([np.nanmean(vels) / (distances[j] * del_phi) for j, vels in enumerate(v_rad)])
    v_phi_ave = transpose_unyt([np.nanmean(vels) / (distances[j] * del_phi) for j, vels in enumerate(v_phi)])
    t_ave = transpose_unyt([np.nanmean(temperature) for temperature in temps])

    v_para_std = transpose_unyt([np.nanstd(vels) / (distances[j] * del_phi) for j, vels in enumerate(v_para)])
    v_rad_std = transpose_unyt([np.nanstd(vels) / (distances[j] * del_phi) for j, vels in enumerate(v_rad)])
    v_phi_std = transpose_unyt([np.nanstd(vels) / (distances[j] * del_phi) for j, vels in enumerate(v_phi)])
    t_std = transpose_unyt([np.nanstd(temperature) for temperature in temps])

    v_para_diff = transpose_unyt([np.nanmean(transpose_unyt([v_para[j+1][n] - v_para[j][n] for n in range (NUM_RAYS)])) / (distances[j+1] - distances[j]) for j in range (arr_len - 1)])
    v_rad_diff = transpose_unyt([np.nanmean(transpose_unyt([v_rad[j+1][n] - v_rad[j][n] for n in range (NUM_RAYS)])) / (distances[j+1] - distances[j]) for j in range (arr_len - 1)])
    v_phi_diff = transpose_unyt([np.nanmean(transpose_unyt([v_phi[j+1][n] - v_phi[j][n] for n in range (NUM_RAYS)])) / (distances[j+1] - distances[j]) for j in range (arr_len - 1)])
    t_diff = transpose_unyt([np.nanmean(transpose_unyt([temps[j+1][n] - temps[j][n] for n in range (NUM_RAYS)])) / (distances[j+1] - distances[j]) for j in range (arr_len - 1)])

    cooling_times = transpose_unyt([np.nanmean(transpose_unyt(times)) for times in zip(*[unyt_array(df[f'cooling_time_{n}/array_data'], 's') for n in range (NUM_RAYS)])])
    densities = transpose_unyt([np.nanmean(transpose_unyt(dens)) for dens in zip(*[unyt_array(df[f'density_{n}/array_data'], 'g/cm**3') for n in range (NUM_RAYS)])])
    
    numer = (1/(GAMMA-1)) * ((t_ave[:-1] / cooling_times[:-1]) + v_rad_ave[:-1] * t_std[:-1] + v_para_ave[:-1] * t_diff) + t_ave[:-1] * (v_rad_std[:-1] + v_rad_ave[:-1] + v_para_diff)
    denom = v_rad_std[:-1]**2 + v_para_diff**2 + 0.5 * v_rad_diff**2 + 0.5 * (v_phi_std[:-1] - v_phi_ave[:-1])**2 + 0.5 * (v_rad_diff + v_para_std[:-1])**2 + 0.5 * v_phi_diff**2 
    
    nu = (0.5 * CONV_CONST * (numer / denom)).to('m**2/s')
    #mu = (densities[:-1] * nu).to('Pa*s')
    #shear = denom.to('s**(-2)')

    plt.figure()
    plt.title(get_title(dump_name, star_type))
    plt.plot(distances[:-1], nu)
    plt.ylabel(r"$\nu ~(m^{2} ~s^{-1})$")
    plt.xlabel("Radius (pc)")
    plt.xlim(0.0, 400.0)
    plt.yscale('log')
    plt.ylim(1e26, 1e38)
    plt.savefig(os.path.join(OUTDIR_NU, f"DD{get_dump_num(dump_name)}_nu.png"))
    plt.close()

    fig1 = plt.figure()
    fig2 = plt.figure()
    ax1 = fig1.add_subplot()
    ax2 = fig2.add_subplot()
    ax1.set_title(get_title(dump_name, star_type))
    ax2.set_title(get_title(dump_name, star_type))
    for j, length in enumerate(LENGTHS):
        index_span = np.argwhere((distances - length) > 0)[0].item()
        index_start = int(index_span / 2)
        index_end = arr_len - int(index_span / 2) - 1
        distances_len = distances[index_start:index_end]
        averages = [0.0 for distance in distances_len]
        integrals = [0.0 for distance in distances_len]
        for k, index in enumerate(range(index_start, index_end)):
            averages[k] = np.mean(nu[index-int(index_span/2):index+int(index_span/2)])
            for l in range (index - int(index_span / 2), index + int(index_span / 2)):
                integrals[k] += nu[l]
        ax1.plot(distances_len, length / averages, color= COLORS[j], label=f"{length}")
        ax2.plot(distances_len, length / integrals, color= COLORS[j], label=f"{length}")
    ax1.set_ylabel(r"$ Re/v ~(s ~m^{-1})$")
    ax2.set_ylabel(r"$ Re/v ~(s ~m^{-1})$")
    ax1.set_xlabel("Radius (pc)")
    ax2.set_xlabel("Radius (pc)")
    ax1.set_xlim(0.0, 400.0)
    ax2.set_xlim(0.0, 400.0)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.set_ylim(1e-20, 1e-9)
    ax2.set_ylim(1e-20, 1e-9)
    ax1.legend(loc='lower right')
    ax2.legend(loc='lower right')
    fig1.savefig(os.path.join(OUTDIR_NU, f"DD{get_dump_num(dump_name)}_nus_ave.png"))
    fig2.savefig(os.path.join(OUTDIR_NU, f"DD{get_dump_num(dump_name)}_nus_int.png"))
    plt.close()
