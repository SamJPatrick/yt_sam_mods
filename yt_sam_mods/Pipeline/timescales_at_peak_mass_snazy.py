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

    field_dict = {"turbulent_sound_crossing_time": 'orange', "sound_crossing_time": 'red',
                  "total_dynamical_time": 'green', "cooling_time": 'blue'}
    file_dict = {'ccsn': 'CCSN', 'hn': 'HN', 'pisn': 'PISN'}
    
    my_fig = GridFigure(3, 1, figsize=(11, 9),
                        left_buffer=0.10, right_buffer=0.02,
                        bottom_buffer=0.08, top_buffer=0.02,
                        vertical_buffer=0.02)
        
    for i, my_axes in enumerate(my_fig):
        star_type = list(file_dict.keys())[i]
        filename = os.path.join("Britton_sim_data", file_dict[star_type], "star_None_mass.h5")
        pds = yt.load(filename)
        df_mass = pds.data
        r_index = np.argmax(df_mass[('data', 'bonnor_ebert_ratio')][-1])
        times, t_indicies = get_xaxis_data(star_type, df_mass)
        for field, color in field_dict.items():
            data = transpose_unyt([df_mass[('data', field)][t_index][r_index].to('Myr') for t_index in t_indicies])
            my_axes.plot(times, data, color= color, alpha= 0.8, label= field)
        my_axes.set_xlim(0, 130)
        my_axes.set_yscale("log")
        my_axes.yaxis.set_label_text(f"Timescales {star_type} (Myr)")
        my_axes.set_ylim(1e0, 1e6)
        my_axes.xaxis.set_ticks(np.arange(0, 131, 10))
        my_axes.grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)
        my_axes.tick_params(axis="x", direction="inout", which="both", top=True, bottom=True)
        
    my_fig[0].legend(loc= 'upper right')
    my_fig[2].xaxis.set_label_text("time (Myr)")
    my_fig[0].xaxis.set_ticklabels(['' for i in range (len(np.arange(0, 131, 10)))])
    my_fig[1].xaxis.set_ticklabels(['' for i in range (len(np.arange(0, 131, 10)))])
    
    pyplot.savefig("model_profiles.png")
