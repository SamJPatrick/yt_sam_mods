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
INDIR = "Sobolev/Ray_profiles/Packed_rays"
OUTDIR = "Sobolev/Bary_vs_dm"
DATASET_NUMS = [120, 130, 140, 150, 160, 165]   # <----- MODIFY SELECTED DATA DUMPS HERE!!
#DATASET_NUMS = [98, 100, 110]
DENS_CUTOFF = 1e-28


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass


datasets = sorted([os.path.join(INDIR, f"DD{num:04d}_packed.h5") for num in DATASET_NUMS])
for dataset in datasets:
    plt.figure()
    plt.title(star_type.upper())
    dump_name = os.path.basename(dataset)
    ds = h5py.File(dataset, 'r')
    distances = unyt_array(ds['distances'][:], 'pc')
    dm_arrs = list(zip(*[unyt_array(ds[f'dark_matter_density_{n}'][:], 'g/cm**3') for n in range (NUM_RAYS)]))
    dm_densities = [0.0 for i in range (len(dm_arrs))]
    for i, arr in enumerate(dm_arrs):
        arr_sane = [i for i in arr if i > unyt_quantity(DENS_CUTOFF, 'g/cm**3')]
        if (len(arr_sane) == 0):
            dm_mean = unyt_quantity(0.0, 'g/cm**3')
        else :
            dm_mean = np.mean(transpose_unyt(arr_sane))
        dm_densities[i] = dm_mean
    dm_densities = transpose_unyt(dm_densities)
    bary_densities = transpose_unyt([np.mean(transpose_unyt(x)) for x in zip(*[unyt_array(ds[f'density_{n}'][:], 'g/cm**3') for n in range (NUM_RAYS)])])
    plt.plot(distances, dm_densities, label="DM", color='black')
    plt.plot(distances, bary_densities, label="baryonic", color='red')
    plt.ylabel(r"$\rho (g ~cm^{-3})$")
    plt.xlabel("Radius (pc)")
    #plt.xlim(1e-1, 1e4)
    #plt.xscale('log')
    plt.xlim(0.0, 400.0)
    plt.legend(loc='upper right')
    plt.yscale('log')
    plt.ylim(1e-27, 1e-18)
    plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dump_name)}_bary_vs_dm_density.png"))
    plt.close()
