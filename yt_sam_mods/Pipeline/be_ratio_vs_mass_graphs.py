import os
import yt
import sys
import numpy as np
from  yt.extensions.sam_mods.graph_funcs import *

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt



MASS_FILE = "star_None_mass.h5"
SIMULATION_FILE = "simulation.h5"
OUTDIR = "Profiles/Bonnor_ebert_ratios_mass"


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass

df_sim = yt.load(SIMULATION_FILE)
df_mass = yt.load(MASS_FILE)
for index, time in enumerate(df_mass.data[('data', 'time')].to('Myr')):
    
    index_sim = np.argwhere(df_sim.data[('data', 'time')].to('Myr') == time).item()
    dsfn = df_sim.data[('data', 'filename')].astype(str)[index_sim].split('/')[-1]
    mass = df_mass.data[('data', 'gas_mass_enclosed')][index].to('Msun')
    be_ratio = df_mass.data[('data', 'bonnor_ebert_ratio')][index].to('')

    plt.figure()
    plt.title(get_title(dsfn, star_type))
    plt.plot(mass, be_ratio)
    plt.axhline(y=1.0, color='red')
    plt.xlabel("$Mass (M_{\odot})$")
    plt.xscale('log')
    plt.xlim(1e-2, 1e5)
    plt.ylabel(r"$M_{tot} / M_{BE}$")
    plt.yscale('log')
    plt.ylim(1e-4, 1e1)
    plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dsfn)}_be_ratio_mass.png"))
    plt.close()
