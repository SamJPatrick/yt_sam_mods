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
STAR_DIR = ['CCSN', 'HN', 'PISN']
VEL_DIR = "Profiles/Velocity_profiles"
VEL_FILE = ["DD0179_1D_profile_radius_cell_mass.h5",
            "DD0232_1D_profile_radius_cell_mass.h5",
            "DD0295_1D_profile_radius_cell_mass.h5"]





if __name__ == "__main__":

    
    my_fig = GridFigure(2, 3, figsize=(15, 10),
                        left_buffer=0.08, right_buffer=0.02,
                        bottom_buffer=0.08, top_buffer=0.08,
                        vertical_buffer=0.02, horizontal_buffer=0.08)

    xticks = np.logspace(-2, 2, 5)

    
    for my_axes in my_fig:
        
        my_axes.set_xscale('log')
        my_axes.set_xlim(1e-2, 5e2)
        my_axes.grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
                
    for i, my_axes in enumerate(my_fig.top_axes):
        
        #my_axes.tick_params(axis="x", top=False, bottom=True, labeltop=False, labelbottom=True)
        #my_axes.xaxis.set_label_position("top")
        my_axes.set_ylabel(r"$M_{tot} / M_{BE}$", labelpad=8)
        my_axes.set_yscale('log')
        my_axes.set_ylim(1e-2, 5e1)

        df_radius = yt.load(os.path.join(PRE_DIR, STAR_DIR[i], "star_None_radius.h5"))
        used = df_radius.data[('data', 'used')][-1].astype(bool)
        be_ratios = df_radius.data[('data', 'bonnor_ebert_ratio')][-1][used]
        radii = df_radius.data[('data', 'radius')][1:][used].to('pc').value
        my_axes.plot(radii, be_ratios)

        my_axes.axvline(x= radii[np.argmax(be_ratios)], color='blue')
        my_axes.axhline(y= 1.0, color= 'red')

        my_axes.xaxis.set_label_position("top")
        my_axes.xaxis.set_label_text(STAR_DIR[i], fontsize=18)
        my_axes.tick_params(labelbottom=False)

        #ax_top = my_axes.twiny()
        #ax_top.plot(radii, be_ratios, ls="", color="none")
        masses = df_radius.data[('data', 'gas_mass_enclosed')][-1][used]
        top_labels = [masses[np.argmin(np.abs(radii - num))].to("Msun").value.item() for num in xticks]
        print(xticks)
        print(top_labels)
        #plt.xticks(ticks= xticks, labels= [f'{num:.2e}' for num in top_labels])
        #ax_top.xaxis.set_major_locator(ax.xaxis.get_major_locator())
        #ax_top.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{k*x:g}"))


    for i, my_axes in enumerate(my_fig.bottom_axes):

        #my_axes.tick_params(axis="x", top=False, bottom=True, labeltop=True, labelbottom=True)
        #my_axes.xaxis.set_label_position("top")
        my_axes.set_xlabel("Radius (pc)", labelpad=8)
        my_axes.set_ylabel(r"Radial velocity (km s$^{-1}$)", labelpad=8)
        if (i == 2):
            my_axes.set_ylim(-1.0, 12.0)
        else :
            my_axes.set_ylim(-1.0, 3.0)

        df_vel = yt.load(os.path.join(PRE_DIR, STAR_DIR[i], VEL_DIR, VEL_FILE[i]))
        used = df_vel.data[('data', 'used')].astype(bool)
        velocities = - df_vel.data[('gas', 'radial_velocity')][used][20:].to('km/s')
        radii = df_vel.data[('data', 'radius')][used][20:].to('pc')
        my_axes.plot(radii, velocities)

        my_axes.axvline(x= radii[np.argmin(np.abs(velocities))], color='blue')
        my_axes.axhline(y= 0.0, color= 'red')
        
        
    plt.savefig("model_profiles.png")
