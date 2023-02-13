import os
import yt
import sys
import numpy as np
from  yt.extensions.sam_mods.graph_funcs import *

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt



RADIUS_FILE = "star_None_radius.h5"
SIMULATION_FILE = "simulation.h5"
OUTDIR = "Profiles/Bonnor_ebert_ratios_radius"


try :
    star_type = sys.argv[1]
except IndexError:
    star_type = ""
    pass

df_sim = yt.load(SIMULATION_FILE)
df_radius = yt.load(RADIUS_FILE)
for index, time in enumerate(df_radius.data[('data', 'time')].to('Myr')):
    
    index_sim = np.argwhere(df_sim.data[('data', 'time')].to('Myr') == time).item()
    dsfn = df_sim.data[('data', 'filename')].astype(str)[index_sim].split('/')[-1]
    radius = df_radius.data[('data', 'radius')][1:].to('pc')
    be_ratio = df_radius.data[('data', 'bonnor_ebert_ratio')][index].to('')

    plt.figure()
    plt.title(get_title(dsfn, star_type))
    plt.plot(radius, be_ratio)
    plt.axhline(y=1.0, color='red')
    plt.xlabel("$Radius (pc)$")
    plt.xscale('log')
    plt.xlim(1e-2, 5e2)
    plt.ylabel(r"$M_{tot} / M_{BE}$")
    plt.yscale('log')
    plt.ylim(1e-4, 1e1)
    plt.savefig(os.path.join(OUTDIR, f"DD{get_dump_num(dsfn)}_be_ratio_radius.png"))
    plt.close()
