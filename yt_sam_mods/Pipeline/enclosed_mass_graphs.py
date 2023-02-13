import os
import sys
import h5py
import glob
import yt
from unyt import Msun, pc
from unyt import unyt_array
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from yt.extensions.sam_mods.graph_funcs import get_dump_num, get_title, get_time_z


INDIR = "Profiles/Normal_profiles"
OUTDIR = "Profiles/Enclosed_mass_graphs"


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass

filenames = glob.glob(os.path.join(os.getcwd(), INDIR, "DD*_profile_weight_field_None.h5"))
for path in filenames:

    dump_name = os.path.basename(path)
    title = get_title(dump_name, star_type)
    dump = f"DD{get_dump_num(dump_name)}"
    normdf = yt.load(path)
    
    used = normdf.profile.used
    radii = normdf.profile.x[used].to("pc")
    tot_mass = normdf.profile[('data', 'matter_mass')][used].in_units('Msun')
    bary_mass = normdf.profile[('data', 'cell_mass')][used].in_units('Msun')
    metal_mass = normdf.profile[('data', 'metal3_mass')][used].in_units('Msun')
    if ((len(tot_mass) != len(bary_mass))
        or (len(bary_mass) != len(metal_mass))
        or (len(tot_mass) != len(metal_mass))):
        print("Error, profiles are of unequal length")
        continue

    tot_mass_acml = [sum(tot_mass[:i+1]) for i in range(len(tot_mass))]
    bary_mass_acml = [sum(bary_mass[:i+1]) for i in range(len(bary_mass))]
    metal_mass_acml = [sum(metal_mass[:i+1]) for i in range(len(metal_mass))]

    plt.figure()
    plt.title(title)
    plt.plot(radii, tot_mass_acml, label='total', color='black')
    plt.plot(radii, bary_mass_acml, label='baryonic', color='red')
    plt.plot(radii, metal_mass_acml, label='metals', color='blue')
    plt.xlabel("Radius (pc)")
    plt.ylabel(r"$Mass ~(M_{\odot})$")
    plt.xscale('log')
    plt.xlim(1e-2, 5e2)
    plt.yscale('log')
    plt.ylim(1e-9, 1e7)
    plt.legend(loc='upper left')
    plt.savefig(os.path.join(OUTDIR, f"{dump}_enclosed_masses.png"))
    plt.close()    

    '''
    bary_frac = bary_mass / tot_mass
    metal_frac = metal_mass / bary_mass

    plt.figure()
    plt.title(title)
    _, ax1 = plt.subplots()
    ax1.plot(radii, bary_frac, label='baryon fraction', color='red')
    ax1.set_xlabel("Radius (pc)")
    ax1.set_xscale('log')
    ax1.set_xlim(1e-2, 5e2)
    ax1.set_ylabel(r"$\Omega_b$", color='red')
    ax1.set_ylim(0.0, 1.0)
    ax2 = ax1.twinx()
    ax2.plot(radii, metal_frac, label='metal fraction', color='blue')
    ax2.set_ylabel(r"$Z$", color='blue')
    ax2.set_yscale('log')
    ax2.set_ylim(1e-7, 1e-1)
    plt.title(title)
    plt.savefig(os.path.join(OUTDIR, f"{dump}_enclosed_mass_fractions.png"))
    plt.close()
    '''
