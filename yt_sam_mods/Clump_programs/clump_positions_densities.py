import numpy as np
import os
import yt

from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter
pyplot.rcParams['font.size'] = 14

from grid_figure import GridFigure
from yt.extensions.sam_mods.clump_funcs import Node, sanatise_clumps
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from unyt import unyt_quantity, unyt_array
from unyt import mh


MU = 1.6

RADIUS_DIR = "Britton_sim_data"
CLUMP_DIR = "Britton_sim_data/Clump_data/Clumps"
FILENAME = "star_None_radius.h5"
FILE_DICT = {'ccsn': 'DD0179_clump_info_ccsn.h5',
             'hn': 'DD0232_clump_info_hn.h5',
             'pisn': 'DD0295_clump_info_pisn.h5'}
STAR_TYPES = ['CCSN', 'HN', 'PISN']
CLUMP_FILES = ['DD0179_clump_info_ccsn.h5', 'DD0232_clump_info_hn.h5', 'DD0295_clump_info_pisn.h5']


def _int_fmt(t, pos):
    return f"{t:d}"

def _flt_fmt(t, pos):
    return np.format_float_positional(t, trim="-")


my_fig = GridFigure(2, 3, figsize=(16, 9),
                    left_buffer=0.08, right_buffer=0.02,
                    bottom_buffer=0.08, top_buffer=0.08,
                    vertical_buffer=0.04, horizontal_buffer=0.08)


for i, my_axes in enumerate(my_fig.top_axes):
    
    df_radius = yt.load(os.path.join(RADIUS_DIR, STAR_TYPES[i], FILENAME))
    radii = df_radius.data[('data', 'radius')][1:].to('pc')
    densities = (df_radius.data[('data', 'density')][-1] / (MU * mh)).to('cm**(-3)')
    index_be = np.argmax(df_radius.data[('data', 'bonnor_ebert_ratio')][-1])
    radius_be = radii[index_be]
    my_axes.plot(radii, densities)
    
    my_axes.set_ylabel(r"Gas density (cm$^{-3}$)")
    my_axes.axvline(x= radius_be)
    my_axes.set_xlabel(r"Radius (pc)")

    my_axes.xaxis.set_label_position("top")
    my_axes.xaxis.set_label_text(STAR_TYPES[i], fontsize=18)
    my_axes.tick_params(labelbottom=False)


for i, my_axes in enumerate(my_fig.bottom_axes):

    df = yt.load(os.path.join(CLUMP_DIR, CLUMP_FILES[i]))
    #sane_leaves = sanatise_clumps(df.leaves)
    sane_leaves = df.leaves

    positions = np.array([clump[('clump', 'com')] for clump in sane_leaves])
    masses = np.array([clump[('clump', 'cell_mass')] for clump in sane_leaves])
    pos_start, com, mean_inter, mean_com = Node.set_threshold(positions, masses)
    master_clump = Node(pos_start)
    master_clump.search_first_right()
    node_list = np.array(master_clump.produce_tree_list())
    
    densities = transpose_unyt([clump[('clump', 'max_number_density')] for clump in sane_leaves])
    distances = unyt_array([np.linalg.norm(clump[('clump', 'com')] - com) for clump in sane_leaves], 'pc')
    densities_slct = transpose_unyt([clump[('clump', 'max_number_density')] for clump in sane_leaves])[node_list]
    distances_slct = unyt_array([np.linalg.norm(clump[('clump', 'com')] - com) for clump in sane_leaves], 'pc')[node_list]

    dots_all = my_axes.scatter(distances, densities, c='blue')
    dots_slct = my_axes.scatter(distances_slct, densities_slct, c='red')
    my_axes.set_xlabel("Distance from COM (pc)")
    my_axes.set_ylabel("Clump density (cm$^{-3}$)")
    my_axes.legend([dots_slct], [f'{x:.3f} pc' for x in [mean_inter]], loc='upper right', fontsize=10)
    my_axes.axvline(x= mean_com, linestyle='--')


for my_axes in my_fig:

    my_axes.set_xlim(5e-1, 5e1)
    my_axes.set_xscale('log')
    my_axes.set_ylim(1e0, 1e9)
    my_axes.set_yscale('log')
    my_axes.grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)


pyplot.savefig("clump_positions_densities_test.png")
