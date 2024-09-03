import numpy as np
import os
import yt

from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
plt.rcParams['font.size'] = 14

from grid_figure import GridFigure
from yt.extensions.sam_mods.clump_funcs import Node, sanatise_clumps
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from unyt import unyt_quantity, unyt_array
from unyt import mh


MU = 1.6

PRE_DIR = "Britton_sim_data"
POST_DIR = "Profiles/Normal_profiles"
STAR_TYPES = ['CCSN', 'HN', 'PISN']
CLUMP_DIR = "Britton_sim_data/Clump_data/Clumps"
CLUMP_FILES = ['DD0179_clump_info_ccsn.h5',
               'DD0232_clump_info_hn.h5',
               'DD0295_clump_info_pisn.h5']
FILENAMES = ["DD0179_profile_weight_field_cell_mass.h5",
             "DD0232_profile_weight_field_cell_mass.h5",
             "DD0295_profile_weight_field_cell_mass.h5"]    


def _int_fmt(t, pos):
    return f"{t:d}"

def _flt_fmt(t, pos):
    return np.format_float_positional(t, trim="-")


my_fig = GridFigure(2, 3, figsize=(16, 9),
                    left_buffer=0.06, right_buffer=0.08,
                    bottom_buffer=0.08, top_buffer=0.08,
                    vertical_buffer=0.04, horizontal_buffer=0.10)


for i, my_axes in enumerate(my_fig.top_axes):
    
    df = yt.load(os.path.join(PRE_DIR, STAR_TYPES[i], POST_DIR, FILENAMES[i]))
    used = df.profile.used
    radii = df.profile.x[used].to("pc")
    sig_z = df.profile.standard_deviation[('data', 'metallicity3')][used]
    metallicity = df.profile[('data', 'metallicity3')][used].to('Zsun').value

    my_axes.plot(radii, sig_z, c='orange')
    my_axes.set_xlabel("Radius (pc)")
    my_axes.set_ylabel(r"$\sigma_Z$ (Z$_{\odot}$)")
    my_axes.set_ylim(1e-6, 1e-1)
    my_axes.set_yscale('log')
    my_axes.yaxis.label.set_color('orange')
    my_axes.tick_params(axis='y', colors='orange')

    my_axes_dual = my_axes.twinx()
    my_axes_dual.plot(radii, metallicity, c='blue')
    my_axes_dual.set_ylabel(r"Z (10$^{-3}$ Z$_{\odot}$)")
    my_axes_dual.set_ylim(0.0, 1e-3)
    my_axes_dual.yaxis.label.set_color('blue')
    #plt.yticks(ticks = np.linspace(0, 1e-3, 6), labels = [f'{num:.1e}' for num in np.linspace(0, 1e-3, 6)])
    plt.yticks(ticks= np.linspace(0, 1e-3, 11), labels= [f'{num:.1f}' for num in np.linspace(0, 1.0, 11)])
    my_axes_dual.tick_params(axis='y', colors='blue')

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
    
    metallicities = transpose_unyt([clump[('clump', 'max_metallicity')] for clump in sane_leaves])
    distances = unyt_array([np.linalg.norm(clump[('clump', 'com')] - com) for clump in sane_leaves], 'pc')
    metallicities_slct = transpose_unyt([clump[('clump', 'max_metallicity')] for clump in sane_leaves])[node_list]
    distances_slct = unyt_array([np.linalg.norm(clump[('clump', 'com')] - com) for clump in sane_leaves], 'pc')[node_list]

    dots_all = my_axes.scatter(distances, metallicities, c='blue')
    dots_slct = my_axes.scatter(distances_slct, metallicities_slct, c='red')
    my_axes.set_xlabel("Distance from COM (pc)")
    my_axes.set_ylabel(r"Clump metallicty ( 10$^{-3}$ Z$_{\odot}$)")
    my_axes.set_ylim(0.0, 1e-3)
    #plt.yticks(ticks = np.linspace(0, 1e-3, 6), labels = [f'{num:.1e}' for num in np.linspace(0, 1e-3, 6)])
    plt.yticks(ticks= np.linspace(0, 1e-3, 11), labels= [f'{num:.1f}' for num in np.linspace(0, 1.0, 11)])
    my_axes.legend([dots_slct], [f'{x:.3f} pc' for x in [mean_inter]], loc='upper right', fontsize=10)
    my_axes.axvline(x= mean_com, linestyle='--')


for my_axes in my_fig:

    my_axes.set_xlim(5e-1, 5e1)
    my_axes.set_xscale('log')
    my_axes.grid(visible=True, axis="both", zorder=0, linestyle=":", color="black", alpha=0.6)


plt.savefig("clump_positions_metallicities_test.png")
