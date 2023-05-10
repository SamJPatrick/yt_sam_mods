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
idx_max = np.argmax(df_mass.data[('data', 'bonnor_ebert_ratio')][-1])
#idx_cross = np.argwhere(df_mass.data[('data', 'bonnor_ebert_ratio')][-1] > 1)[0].item()
for index, time in enumerate(df_mass.data[('data', 'time')].to('Myr')):
    
    index_sim = np.argwhere(df_sim.data[('data', 'time')].to('Myr') == time).item()
    dsfn = df_sim.data[('data', 'filename')].astype(str)[index_sim].split('/')[-1]
    mass = df_mass.data[('data', 'gas_mass_enclosed')][-1].to('Msun')
    be_ratio = df_mass.data[('data', 'bonnor_ebert_ratio')][index].to('')
    mass_be = mass[idx_max].to('Msun')
    #mass_be_cross = mass[idx_cross].to('Msun')

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
    plt.axvline(x= mass_be, color='blue')
    #plt.axvline(x= mass_be_cross, color='orange')
    plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dsfn)}_be_ratio_mass.png"))
    plt.close()
