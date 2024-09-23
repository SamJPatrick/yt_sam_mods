from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import os
import yt

plt.rcParams['font.size'] = 14

from grid_figure import GridFigure
from yt.extensions.p2p.stars import get_star_data
from yt.utilities.physical_constants import G
from yt.visualization.color_maps import yt_colormaps
from unyt import uvstack

def _int_fmt(t, pos):
    return f"{t:d}"

def _flt_fmt(t, pos):
    return np.format_float_positional(t, trim="-")



PRE_DIR = "Britton_sim_data"
STAR_TYPES = ['CCSN', 'HN', 'PISN']
COLORS = ['blue', 'orange', 'green']


if __name__ == "__main__":

    
    my_fig = GridFigure(1, 2, figsize=(12, 8),
                        left_buffer=0.08, right_buffer=0.02,
                        bottom_buffer=0.08, top_buffer=0.08,
                        vertical_buffer=0.08, horizontal_buffer=0.08)

    xticks = np.logspace(-2, 2, 5)

    
    for my_axes in my_fig:        
        my_axes.set_xscale('log')
        my_axes.set_ylim(1e-2, 5e1)
        my_axes.set_ylabel(r"$M_{tot} / M_{BE}$", labelpad=8)
        my_axes.set_yscale('log')
        my_axes.grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
        my_axes.axhline(y= 1.0, color= 'red')
        
    my_fig[0].set_xlim(1e-2, 5e2)
    my_fig[0].set_xlabel("Radius (pc)")
    for i, star_type in enumerate(STAR_TYPES):
        df_radius = yt.load(os.path.join(PRE_DIR, star_type.upper(), "star_None_radius.h5"))
        used = df_radius.data[('data', 'used')][-1].astype(bool)
        be_ratios = df_radius.data[('data', 'bonnor_ebert_ratio')][-1][used]
        radii = df_radius.data[('data', 'radius')][1:][used].to('pc').value
        my_fig[0].plot(radii, be_ratios, color= COLORS[i], label= star_type)
        my_fig[0].axvline(x= radii[np.argmax(be_ratios)], color= COLORS[i], linestyle= '--', alpha=0.6)

    my_fig[1].set_xlim(1e-2, 1e5)
    my_fig[1].set_xlabel(r"Mass (M$_{\odot}$)")
    for i, star_type in enumerate(STAR_TYPES):
        df_mass = yt.load(os.path.join(PRE_DIR, star_type.upper(), "star_None_mass.h5"))
        used = df_mass.data[('data', 'used')][-1].astype(bool)
        be_ratios = df_mass.data[('data', 'bonnor_ebert_ratio')][-1][used]
        masses = df_mass.data[('data', 'gas_mass_enclosed')][-1][used].to('Msun').value
        my_fig[1].plot(masses, be_ratios, color= COLORS[i], label= star_type)
        my_fig[1].axvline(x= masses[np.argmax(be_ratios)], color= COLORS[i], linestyle= '--', alpha=0.6)
    my_fig[1].legend(loc='upper left')


    plt.savefig("model_profiles.png")
