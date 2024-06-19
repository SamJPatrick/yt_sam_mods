from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter
import numpy as np
import os
import yt

from grid_figure import GridFigure
from yt.visualization.color_maps import yt_colormaps

from yt.extensions.sam_mods.graph_funcs import get_time_offset
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt


pyplot.rcParams['font.size'] = 12



def _int_fmt(t, pos):
    return f"{t:d}"

def _flt_fmt(t, pos):
    return np.format_float_positional(t, trim="-")

def get_xaxis_data(star_type, df_mass):
    time_offset = get_time_offset(star_type)
    otimes = df_mass[('data', 'time')]
    times = transpose_unyt([time.to('Myr') - time_offset for time in otimes if time >= time_offset])
    t_indicies = np.arange(np.argwhere(otimes >= time_offset)[0].item(), len(otimes))
    return times, t_indicies



if __name__ == "__main__":

    fields = ["density", "temperature", "metallicity3", "bonnor_ebert_ratio"]
    color_dict = {'ccsn': 'blue', 'hn': 'orange', 'pisn': 'green'}
    file_dict = {'ccsn': 'CCSN', 'hn': 'HN', 'pisn': 'PISN'}
    
    my_fig = GridFigure(2, 2, figsize=(11, 9),
                    left_buffer=0.10, right_buffer=0.10,
                    bottom_buffer=0.08, top_buffer=0.08,
                    vertical_buffer=0.02, horizontal_buffer=0.04)
        
    for i, my_axes in enumerate(my_fig):
        for star_type in color_dict.keys():
            filename = os.path.join("Britton_sim_data", file_dict[star_type], "star_None_mass.h5")
            pds = yt.load(filename)
            df_mass = pds.data
            r_index = np.argmax(df_mass[('data', 'bonnor_ebert_ratio')][-1])
            times, t_indicies = get_xaxis_data(star_type, df_mass)
            data = [df_mass[('data', fields[i])][t_index][r_index] for t_index in t_indicies]
            my_axes.plot(times, data, color= color_dict[star_type], alpha= 0.8, label= star_type)
        my_axes.set_yscale("log")
        my_axes.xaxis.set_label_text("time (Myr)")
        my_axes.grid(visible=True, axis="both", zorder=0, linestyle=":",color="black", alpha=0.6)
        my_axes.xaxis.set_ticks(np.arange(0, 131, 10))
        my_axes.set_xlim(0, 130)
        my_axes.tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)

    my_fig[0].yaxis.set_label_text("Density (g cm$^{-3}$)")
    my_fig[0].set_ylim(1e-27, 1e-18)
    my_fig[0].yaxis.set_ticks(np.logspace(-27, -18, 10), minor= True)

    my_fig[1].yaxis.set_label_text("Temperature (K)")
    my_fig[1].set_ylim(1e1, 1e5)
    my_fig[1].legend(loc= 'upper right')

    my_fig[2].yaxis.set_label_text("Metallicity ($Z_{\odot}$)")
    my_fig[2].set_ylim(1e-7, 1e0)
    my_fig[2].axhline(10**(-5.3), color='red', linestyle='--')
    
    my_fig[3].yaxis.set_label_text("M$_{\\rm gas, enc}$ / M$_{\\rm BE}$")
    my_fig[3].set_ylim(1e-4, 1e2)
    my_fig[3].axhline(y= 1.0, color='red', linestyle='--')    
    
    for my_axes in my_fig.left_axes:
        my_axes.tick_params(axis="y", left=True, direction="inout", which="both")
        my_axes.tick_params(axis="y", right=True, direction="in", which="both")
    for my_axes in my_fig.right_axes:
        my_axes.tick_params(axis="y", right=True, direction="inout", which="both",
                            labelright=True)
        my_axes.tick_params(axis="y", left=True, direction="in", which="both",
                            labelleft=False)
        my_axes.yaxis.set_label_position("right")
    for my_axes in my_fig.top_axes:
        my_axes.tick_params(axis="x", top=True, bottom=True,
                            labeltop=True, labelbottom=False)
        my_axes.xaxis.set_label_position("top")

    pyplot.savefig("model_profiles.png")
