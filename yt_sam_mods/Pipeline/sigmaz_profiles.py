import os
import glob
import yt
import sys
from unyt import unyt_array

import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt



SIM_FILE = "Britton_sim_data/CCSN/simulation.h5"
NORMDIR = "Britton_sim_data/CCSN/Profiles/Normal_profiles"
OUTDIR = '.'


try:
    star_mode = sys.argv[1]
except IndexError:
    star_mode = ""
    pass

fn_norms = glob.glob(os.path.join(os.getcwd(), NORMDIR, "DD*_profile_weight_field_cell_mass.h5"))
for norm_file in fn_norms:
    
    name = os.path.basename(norm_file)
    title = get_title(name, star_mode, SIM_FILE)
    dump = f"DD{get_dump_num(name)}"
    
    normdf = yt.load(norm_file)    
    used = normdf.profile.used
    radii = normdf.profile.x[used].to("pc")
    sig_z = normdf.profile.standard_deviation[('data', 'metallicity3')][used]
    metallicity = normdf.profile[('data', 'metallicity3')][used]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.title(title)
    ax1.plot(radii, sig_z, c='orange')
    ax1.set_xlabel("Radius (pc)")
    ax1.set_xscale('log')
    ax1.set_xlim(5e-1, 5e2)
    ax1.set_ylabel(r"$\sigma_Z$ (Z$_{\odot}$)")
    ax1.set_yscale('log')
    ax1.set_ylim(1e-8, 1e-1)
    ax1.yaxis.label.set_color('orange')
    ax1.tick_params(axis='y', colors='orange')
    ax2 = ax1.twinx()
    ax2.plot(radii, metallicity, c='blue')
    ax2.set_ylabel(r"Z (Z$_{\odot}$)")
    ax2.set_yscale('log')
    ax2.set_ylim(1e-7, 1e-1)
    ax2.yaxis.label.set_color('blue')
    ax2.tick_params(axis='y', colors='blue')
    fig.savefig(os.path.join(OUTDIR, f"{dump}_sigmaz_profile.png"))
    plt.close()
