import os
import glob
import re
import yt
from unyt import unyt_array
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import get_dump_num, get_time_z, get_title



NORMDIR = "Profiles/Normal_profiles"
VELDIR = "Profiles/Velocity_and_timescale_profiles"
OUTDIR = "Profiles/Accretion_vs_sigma_profiles"


fn_norm = glob.glob(os.path.join(os.getcwd(), NORMDIR, "DD*_profile_weight_field_cell_mass.h5"))
fn_norm.sort(key= lambda x: int(get_dump_num(os.path.basename(x))))
fn_vel = glob.glob(os.path.join(os.getcwd(), VELDIR, "DD*_1D_profile_radius_cell_mass.h5"))
fn_vel.sort(key= lambda x: int(get_dump_num(os.path.basename(x))))
assert len(fn_norm) == len(fn_vel), "Error, directories must have equal number of nodes"

for i in range(len(fn_norm)):
    
    name = os.path.basename(fn_norm[i])
    title = get_title(name)
    dump = f"DD{get_dump_num(name)}"
    normdf = yt.load(fn_norm[i])
    veldf = yt.load(fn_vel[i])
    
    used = normdf.profile.used
    radii = normdf.profile.x[used].to("pc")
    acc = normdf.profile[('data', 'accretion_rate')][used]
    acc_z = normdf.profile[('data', 'accretion_rate_z')][used]
    var_z = normdf.profile.standard_deviation[('data', 'metallicity3')][used]
    var_vr = veldf.profile.standard_deviation[('data', 'velocity_spherical_phi')][used]
    turb = var_v = veldf.profile.standard_deviation[('data', 'velocity_magnitude')][used]
    #p_turb = normdf.profile[('data', 'density')][used] * veldf.profile.standard_deviation[('data', 'velocity_magnitude')][used]**2

    plt.figure()
    plt.title(title)
    _, ax1 = plt.subplots()
    ax1.plot(radii, abs(acc), color='blue', label='all', linestyle= '-')
    ax1.plot(radii, abs(acc_z), color='blue', label='metals', linestyle= '-.')
    ax1.set_xlabel("Radius (pc)")
    ax1.set_xscale('log')
    ax1.set_xlim(1e-2, 5e2)
    ax1.set_ylabel(r"$\dot{M} ~(M_{\odot} ~yr^{-1})$", color='blue')
    ax1.set_yscale('log')
    ax1.set_ylim(1e-13, 2e-2)
    ax1.legend(loc='upper left')
    ax2 = ax1.twinx()
    ax2.plot(radii, var_v, color='red', label= r'$v$', linestyle= '-')
    ax2.plot(radii, var_vr, color='red', label= r'$v_r$', linestyle= '-.')
    ax2.set_ylabel(r"$\sigma$", color='red')
    ax2.set_yscale('log')
    ax2.set_ylim(1e1, 3e6)
    ax2.legend(loc='lower left')
    plt.savefig(os.path.join(OUTDIR, f"{dump}_accretion_vs_sigma.png"))
    plt.close()
    
    plt.figure()
    plt.title(title)
    _, ax1 = plt.subplots()
    ax1.plot(radii, abs(acc), color='blue', label='all', linestyle= '-')
    ax1.plot(radii, abs(acc_z), color='blue', label='metals', linestyle= '-.')
    ax1.set_xlabel("Radius (pc)")
    ax1.set_xscale('log')
    ax1.set_xlim(1e-2, 5e2)
    ax1.set_ylabel(r"$\dot{M} ~(M_{\odot} ~yr^{-1})$", color='blue')
    ax1.set_yscale('log')
    ax1.set_ylim(1e-13, 2e-2)
    ax1.legend(loc='upper left')
    ax2 = ax1.twinx()
    ax2.plot(radii, var_z, color='red', label= r'Z')
    ax2.set_ylabel(r"$\sigma$", color='red')
    ax2.set_yscale('log')
    ax2.set_ylim(1e-7, 3e-1)
    ax2.legend(loc='lower left')
    plt.savefig(os.path.join(OUTDIR, f"{dump}_accretion_vs_sigma_metal.png"))
    plt.close()
