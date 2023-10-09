import os
import yt
import sys
import glob
import h5py
import numpy as np
import scipy.signal as sig
from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc



NUM_RAYS = 20
FILT_WINDOW = 11

INDIR = "Sobolev/Ray_profiles"
OUTDIR = "Sobolev/Bary_vs_dm"

DATASET_NUMS = [98, 100, 110, 120]
#DATASET_NUMS = [99, 100, 101, 103, 106]
DATASET_NUMS = [137, 139, 141, 145, 151]
COLORS = ['green', 'blue', 'orange', 'magenta', 'cyan', 'brown']

DENS_CUTOFF = 1e-27


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


datasets = sorted([os.path.join(INDIR, f"DD{num:04d}_packed.h5") for num in DATASET_NUMS])
fig1 = plt.figure()
ax1 = fig1.add_subplot()
for i, dataset in enumerate(datasets):
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot()
    
    dump_name = os.path.basename(dataset)
    time = get_time_z(dump_name, star_type)[0].to('Myr')
    ds = h5py.File(dataset, 'r')
    distances = unyt_array(ds['distances/array_data'], 'pc')
    bary_densities = transpose_unyt([np.mean(transpose_unyt(x)) for x in \
                                     zip(*[unyt_array(ds[f'density_{n}/array_data'], 'g/cm**3') for n in range (NUM_RAYS)])])

    dm_arrs = list(zip(*[unyt_array(ds[f'dark_matter_density_{n}/array_data'], 'g/cm**3') for n in range (NUM_RAYS)]))
    dm_densities = [0.0 for j in range (len(dm_arrs))]
    for j, arr in enumerate(dm_arrs):
        arr_sane = [j for j in arr if j > unyt_quantity(DENS_CUTOFF, 'g/cm**3')]
        if (len(arr_sane) == 0):
            dm_mean = unyt_quantity(0.0, 'g/cm**3')
        else :
            dm_mean = np.mean(transpose_unyt(arr_sane))
        dm_densities[j] = dm_mean
    dm_densities = transpose_unyt(dm_densities)
    
    ax2.plot(distances, dm_densities, label="DM", color='black')
    ax2.plot(distances, bary_densities, label="baryonic", color='red')

    ax2.set_title(star_type.upper())
    ax2.set_ylabel(r"$\rho ~(g ~cm^{-3})$")
    ax2.set_xlabel("Radius (pc)")
    ax2.set_xlim(0.0, 400.0)
    ax2.set_yscale('log')
    ax2.set_ylim(1e-27, 1e-18)
    ax2.legend(loc='upper right')
    fig2.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dump_name)}_bary_vs_dm_density.png"))

    ax1.plot(distances, sig.medfilt(dm_densities, kernel_size= FILT_WINDOW), color= COLORS[i], linestyle='--')
    ax1.plot(distances, bary_densities, color= COLORS[i], linestyle='-', label=f'{time:.2f}')

ax1.set_ylabel(r"$\rho ~(g ~cm^{-3})$")
ax1.set_xlabel("Radius (pc)")
ax1.set_xlim(0.0, 400.0)
ax1.legend(loc='upper right')
ax1.set_yscale('log')
ax1.set_ylim(1e-27, 1e-18)
fig1.savefig(os.path.join(OUTDIR, f"bary_vs_dm_density_combined.png"))
plt.close()
