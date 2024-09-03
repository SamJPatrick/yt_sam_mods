import os
import yt
import sys
import glob
import numpy as np
import h5py
from  yt.extensions.sam_mods.graph_funcs import *

import pdb

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc


LAST_OUTPUTS = {'pisn': ('DD0295_velocity_mag_fourier_spectra.h5', 'DD0295_1D_profile_radius_cell_mass.h5'), \
                'hn': ('DD0232_velocity_mag_fourier_spectra.h5', 'DD0232_1D_profile_radius_cell_mass.h5'), \
                'ccsn': ('DD0179_velocity_mag_fourier_spectra.h5', 'DD0179_1D_profile_radius_cell_mass.h5')}
PRE_DIR = "Britton_sim_data"
VEL_DIR = "Profiles/Velocity_profiles"
SIM_FILE = "Britton_sim_data/PISN/simulation.h5"
OUTDIR = "."
#INDIR = "Fourier_spectra_mag_hn"
INDIR = "."
COLORS = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'yellow', 'black']
STYLES = ['-', '--', '-.']
X_TICK_GAP = 50
#DATASET_GAP = 3
DATASET_GAP = 10


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = "pisn"
    pass

'''
dsfns = glob.glob(os.path.join(INDIR, "DD*_velocity_mag_fourier_spectra.h5"))
dsfns.sort()
dsfns.reverse()
print(dsfns)
plt.figure()
for i, dsfn in enumerate(dsfns):
    df = h5py.File(dsfn, 'r')
    x_data = 1 / df['x_data/array_data'][1:]
    y_glob = np.array(df['velocity_fourier_globmean/array_data'][1:1025]) / 1e3
    y_diag = np.array(df['velocity_fourier_diagonals/array_data'][1:1025]) / 1e3
    time = get_time_z(dsfn.split('/')[-1], star_type, sim_file= SIM_FILE)[0]
    try :
        plt.plot(x_data, y_glob, color= COLORS[i], alpha= 0.8, label= f't={time:.2f}')
        plt.plot(x_data, y_diag, color= COLORS[i], alpha= 0.8, label= f't={time:.2f}', linestyle='--')
    except IndexError:
        break
plt.legend(loc='upper left')
plt.ylabel("Velocity power (km s$^{-1}$)")
plt.ylim(1e0, 1e7)
plt.yscale('log')
plt.xlabel(r"Wavelength (pc$^{-1}$)")
#plt.xticks(x_data[::X_TICK_GAP], [f'{x:.1f}' for x in (1/x_data)[::X_TICK_GAP]], fontsize=8)
plt.xscale('log')
plt.xlim(1e-1, 1e2)
plt.savefig(os.path.join(OUTDIR, f"velocity_fourier_spectra_test.png"))
plt.close()
'''

'''
plt.figure()
for i, (star_type, filenames) in enumerate(LAST_OUTPUTS.items()):
    df = h5py.File(filenames[0], 'r')
    x_data = 1 / df['x_data/array_data'][1:]
    y_glob = np.array(df['velocity_fourier_globmean/array_data'][1:1025]) / 1e6
    y_diag = np.array(df['velocity_fourier_diagonals/array_data'][1:1025]) / 1e6
    try :
        plt.plot(x_data, y_glob, color= COLORS[i], alpha= 0.8, label= f'{star_type}')
        plt.plot(x_data, y_diag, color= COLORS[i], alpha= 0.8, label= f'{star_type}', linestyle='--')
    except IndexError:
        break
plt.legend(loc='upper left')
plt.ylabel("Velocity power (km s$^{-1}$)")
plt.ylim(1e-6, 1e5)
plt.yscale('log')
plt.xlabel(r"Wavelength (pc)")
#plt.xticks(x_data[::X_TICK_GAP], [f'{x:.1f}' for x in (1/x_data)[::X_TICK_GAP]], fontsize=8)
plt.xscale('log')
plt.xlim(1e-1, 1e2)
plt.savefig(os.path.join(OUTDIR, f"velocity_fourier_spectra_comb.png"))
plt.close()
'''


for star_type, filenames in LAST_OUTPUTS.items():

    specfn = filenames[0]
    velfn = filenames[1]
    
    df_mass = yt.load(os.path.join(PRE_DIR, star_type.upper(), "star_None_mass.h5"))
    used = df_mass.data[('data', 'used')][-1].astype(bool)
    ind_max_be = np.argmax(df_mass.data[('data', 'bonnor_ebert_ratio')][-1][used])
    rad_max_be = df_mass.data[('data', 'radius')][-1][used][ind_max_be].to('pc')
    ind_zero_be = np.argmin(np.abs(df_mass.data[('data', 'bonnor_ebert_ratio')][-1][used][ind_max_be:] - 1.0)) + ind_max_be
    rad_zero_be = df_mass.data[('data', 'radius')][-1][used][ind_zero_be].to('pc')
    t_ff = df_mass.data[('data', 'total_dynamical_time')][-1][used][ind_zero_be].to('Myr')

    try :
        df_vel = yt.load(os.path.join(PRE_DIR, star_type.upper(), VEL_DIR, velfn))
    except FileNotFoundError:
        continue
    used = df_vel.data[('data', 'used')].astype(bool)
    ind_be_vel = np.argwhere(df_vel.data[('data', 'radius')][used].to('pc') >= rad_zero_be)[0].item()
    ind_mid_vel = np.argmax(np.abs(df_vel.data[('data', 'radius')][used].to('pc') - rad_max_be))
    max_vel = np.max(np.abs(df_vel.data[('data', 'velocity_spherical_radius')][used][:ind_be_vel])).to('km/s')
    mid_vel_halo = np.abs(df_vel.data[('data', 'velocity_spherical_radius')][used][ind_mid_vel]).to('km/s')

    df_spec = h5py.File(specfn,'r')
    inv_arr = unyt_array(1 / df_spec['x_data/array_data'][1:], 'pc')
    vel_arr = unyt_array(df_spec['velocity_fourier_globmean/array_data'][1:1026] / 1e6, 'km/s')
    ind_rad_max = np.argmin(np.abs(inv_arr - rad_zero_be))
    ind_rad_min = np.argmin(np.abs(vel_arr - max_vel))
    ind_rad_mid = np.argmin(np.abs(inv_arr - rad_max_be))
    mid_vel_spec = vel_arr[ind_rad_mid]

    print(star_type.upper())
    print("Radius of integration min:", inv_arr[ind_rad_min])
    print("Radius of maximal m_be/m:", rad_max_be)
    print("Radius of integration max:", inv_arr[ind_rad_max])
    print("Radial velocity at (m_be/m)max:", mid_vel_halo)
    print("Spectral velocity at (m_be/m)max:", mid_vel_spec)

    t_ave = 0
    for index in range (ind_rad_max, ind_rad_min + 1):
        t_ave += inv_arr[index] / vel_arr[index]
    print("Destruct time:", t_ave.to('Myr'))
    print("Free-fall time:", t_ff.to('Myr'))
    print("Clump ratio:", t_ave / t_ff)
    print("Velocity ratio:", mid_vel_halo / mid_vel_spec)
    print("-------------")
