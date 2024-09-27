import os
import yt
import sys
import glob
import numpy as np
import h5py
from  yt.extensions.sam_mods.graph_funcs import *
from grid_figure import GridFigure

import pdb

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc


LAST_OUTPUTS = {'ccsn': 'DD0179_velocity_mag_fourier_spectra.h5',
                'hn': 'DD0232_velocity_mag_fourier_spectra.h5',
                'pisn': 'DD0295_velocity_mag_fourier_spectra.h5'}
PRE_DIR = "Britton_sim_data"
VEL_DIR = "Profiles/Velocity_profiles"
OUTDIR = "."
INDIR = "Fourier_spectra"
COLORS = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'yellow', 'black']


my_fig = GridFigure(1, 2, figsize=(13, 5),
                    left_buffer=0.10, right_buffer=0.06,
                    bottom_buffer=0.12, top_buffer=0.06,
                    vertical_buffer=0.00, horizontal_buffer=0.08)


for i, (star_type, filename) in enumerate(LAST_OUTPUTS.items()):

    df_four = h5py.File(os.path.join(INDIR, filename), 'r')
    x_data = unyt_array(1 / df_four['x_data/array_data'][1:], 'pc')
    y_glob = unyt_array(np.array(df_four['velocity_fourier_globmean/array_data'][1:1025]) / 1e6, 'km/s')
    y_diag = unyt_array(np.array(df_four['velocity_fourier_diagonals/array_data'][1:1025]) / 1e6, 'km/s')
    my_fig[0].plot(x_data, y_glob, color= COLORS[i], alpha= 0.8, label= f'{star_type}')
    #my_fig[0].plot(x_data, y_diag, color= COLORS[i], alpha= 0.8, label= f'{star_type}', linestyle='--')

    try :
        df_vel = yt.load(os.path.join(PRE_DIR, star_type.upper(), "star_None_mass_velocity_profiles.h5"))
        df_mass = yt.load(os.path.join(PRE_DIR, star_type.upper(), "star_None_mass_normal_profiles.h5"))
    except FileNotFoundError:
        continue

    #used = df_mass.data[('data', 'used')][-1].astype(bool)
    #ind_max_be = np.argmax(df_mass.data[('data', 'bonnor_ebert_ratio')][-1][used])
    #rad_max_be = df_mass.data[('data', 'radius')][-1][used][ind_max_be].to('pc')
    #ind_zero_be = np.argmin(np.abs(df_mass.data[('data', 'bonnor_ebert_ratio')][-1][used][ind_max_be:] - 1.0)) + ind_max_be
    #rad_zero_be = df_mass.data[('data', 'radius')][-1][used][ind_zero_be].to('pc')
    #t_ff = df_mass.data[('data', 'total_dynamical_time')][-1][used][ind_zero_be].to('Myr')
    #plt.axvline(x= rad_zero_be, color= COLORS[i], linestyle= '--')

    ind_max_be = np.argmax(df_mass.data[('data', 'bonnor_ebert_ratio')][-1])
    rad_max_be = df_mass.data[('data', 'radius')][-1][ind_max_be].to('pc')
    ind_zero_be = np.argmin(np.abs(df_mass.data[('data', 'bonnor_ebert_ratio')][-1][ind_max_be:] - 1.0)) + ind_max_be
    rad_zero_be = df_mass.data[('data', 'radius')][-1][ind_zero_be].to('pc')
    t_ff = df_mass.data[('data', 'total_dynamical_time')][-1][ind_zero_be].to('Myr')
    my_fig[0].axvline(x= rad_zero_be, color= COLORS[i], linestyle= '--')

    #vel_max = np.max(np.abs(df_vel.data[('data', 'velocity_spherical_radius')][-1][used][:ind_zero_be])).to('km/s')
    #vel_be_max = np.abs(df_vel.data[('data', 'velocity_spherical_radius')][-1][used][ind_max_be]).to('km/s')
    #plt.axhline(y= vel_max, color= COLORS[i], linestyle= '--')
    #ind_rad_max = np.argmin(np.abs(x_data - rad_zero_be))
    #ind_rad_min = np.argmin(np.abs(y_glob - vel_max))
    #ind_max_be_four = np.argmin(np.abs(x_data - rad_max_be))
    #mid_vel_spec = y_glob[ind_max_be_four].to('km/s')
    #print(f"Velocity ratio {star_type}:", vel_be_max / mid_vel_spec)

    vel_max = np.max(np.abs(df_vel.data[('data', 'velocity_spherical_radius')][-1][:ind_zero_be])).to('km/s')
    vel_be_max = np.abs(df_vel.data[('data', 'velocity_spherical_radius')][-1][ind_max_be]).to('km/s')
    my_fig[0].axhline(y= vel_max, color= COLORS[i], linestyle= '--')
    ind_rad_max = np.argmin(np.abs(x_data - rad_zero_be))
    ind_rad_min = np.argmin(np.abs(y_glob - vel_max))
    ind_max_be_four = np.argmin(np.abs(x_data - rad_max_be))
    mid_vel_spec = y_glob[ind_max_be_four].to('km/s')
    print(f"Velocity ratio {star_type}:", vel_be_max / mid_vel_spec)

    t_ave = 0
    for index in range (ind_rad_max, ind_rad_min + 1):
        t_ave += x_data[index] / y_glob[index]
    print(f"Clumping factor {star_type}:", t_ave / t_ff)

my_fig[0].legend(loc='upper left')
my_fig[0].set_ylabel("Velocity power (km s$^{-1}$)")
my_fig[0].set_ylim(1e-4, 1e5)
my_fig[0].set_yscale('log')
my_fig[0].set_xlabel(r"Wavelength (pc)")
my_fig[0].set_xscale('log')
my_fig[0].set_xlim(1e-1, 1e2)
my_fig[0].grid(visible=True, axis="both", zorder=1, linestyle=":", color="black", alpha=0.6)
my_fig[0].tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)



for i, star_type in enumerate(LAST_OUTPUTS.keys()):
    
    df_vel = yt.load(os.path.join(PRE_DIR, star_type.upper(), "star_None_mass_velocity_profiles.h5"))
    used = df_vel.data[('data', 'used')][-1].astype(bool)
    velocities = - df_vel.data[('data', 'velocity_spherical_radius')][-1][used].to('km/s')
    radii = df_vel.data[('data', 'radius')][-1][used].to('pc')
    my_fig[1].plot(radii, velocities, color= COLORS[i])

my_fig[1].legend(loc='upper left')
my_fig[1].axhline(y= 0.0, color= 'red')
my_fig[1].set_xlabel("Radius (pc)", labelpad=8)
my_fig[1].set_ylabel(r"Radial velocity (km s$^{-1}$)", labelpad=8)
my_fig[1].set_xlim(1e-3, 5e2)
my_fig[1].set_xscale('log')
#my_fig[1].xaxis.set_ticks(np.logspace(5e-3, 5e2, 10))
my_fig[1].grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
my_fig[1].tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)
my_fig[1].set_ylim(-0.5, 12.0)
#my_fig[1].set_ylim(-1.0, 3.0)

plt.savefig("velocities_test.png")
