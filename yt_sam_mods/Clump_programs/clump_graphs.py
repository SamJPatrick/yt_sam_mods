import yt
import numpy as np
from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from fragplot import TreePlot

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt



CLUMP_FILE_NAME = "Clumps/DD0295_clump_info.h5"
TREE_NAME = "clump_tree_test.png"
MATRIX_NAME = "correlations_ave.png"
PLOT_NAME = "histo"

MATRIX_FIELDS = ['alpha', 'ratio', 'volume', 'cell_mass', 'jeans_mass',
                 'mean_metallicity', 'mean_temperature', 'mean_number_density', 'mean_H2_fraction']
#MATRIX_FIELDS = ['cell_mass', 'min_temperature', 'mean_temperature', 'max_temperature',
#                 'min_metallicity', 'mean_metallicity', 'max_metallicity']
PLOT_FIELDS = ['mean_number_density', 'min_number_density', 'max_number_density']

N_BINS_LOG = 20
N_BINS_LIN = 10

ERR_TOL = 0.2
CMB_TOL = unyt_quantity(10.0, 'K')
TEMP_CRIT_H2 = unyt_quantity(200.0, 'K')
N_CRIT_H2 = unyt_quantity(1.0e4, 'cm**(-3)')
VOL_THRESH = unyt_quantity(1.0e-3, 'pc**3')
Z_CRIT = unyt_quantity(1.0e-4, '')
CMB_Z0 = unyt_quantity(2.725, 'K')



gclumps= lambda field, alist: transpose_unyt([clump[('clump', field)] for clump in alist])


def my_node(clump):

    child_clumps = clump.children
    if (child_clumps is None):
        label = ''
    else :
        density = clump[('clump', 'min_number_density')].to('cm**(-3)').value
        label = f'{density:.3g}'

    if (clump[('clump', 'fragmentation_instability')] > 1.0):
        shape = 'square'
    else :
        shape = 'circle'

    my_kwargs = {'label': label, 'fontsize': 16, 'shape': shape}
    return my_kwargs
    

def sanatise_clumps(clump_list):
    
    order_mask = np.argsort([clump[('clump', 'cell_mass')] for clump in clump_list])
    
    vol_mask = (gclumps('volume', clump_list)[order_mask] > VOL_THRESH)
    dens_mask = (gclumps('max_number_density', clump_list)[order_mask] > N_CRIT_H2)
    temp_mask = (gclumps('mean_temperature', clump_list)[order_mask] < TEMP_CRIT_H2)
    cool_mask = (gclumps('min_temperature', clump_list)[order_mask] < temp_cmb + CMB_TOL)
    err_mask = (gclumps('max_number_density', clump_list)[order_mask] < \
                10 * (1 + ERR_TOL) * gclumps('mean_number_density', clump_list)[order_mask])
    sanity_mask = err_mask & vol_mask & temp_mask & cool_mask #& dens_mask
    
    sane_list = np.array(clump_list)[order_mask][sanity_mask]
    duff_frac = (sum(sanity_mask) / len(sanity_mask)) * 100
    print(f"Percentage of clumps remaining is {duff_frac:.2f}%")
    return sane_list


def create_tree_plot(df):

    p = TreePlot(df, dot_kwargs={'rankdir': 'TB', 'size': 18}, node_function= my_node)
    p.size_field = "mass"
    p.color_field = "max_metallicity"
    p.node_cmap = 'turbo'
    p.save(TREE_NAME)


def create_matrix_plot(sane_list):
        
    corr_matrix = np.zeros((len(MATRIX_FIELDS), len(MATRIX_FIELDS)))
    for i, field_a in enumerate(MATRIX_FIELDS):
        for j, field_b in enumerate(MATRIX_FIELDS):
            corr_matrix[i][j] = np.corrcoef(gclumps(field_a, sane_list), gclumps(field_b, sane_list))[0][1]

    plt.figure()
    plt.matshow(np.abs(corr_matrix))
    for i in range (len(MATRIX_FIELDS)):
        for j in range (len(MATRIX_FIELDS)):
            plt.text(i, j, f'{corr_matrix[i][j]:.3f}', ha='center', va='center', color='black', fontsize=5)
    plt.xticks(np.arange(len(MATRIX_FIELDS)), MATRIX_FIELDS, fontsize= 5, rotation='vertical')
    plt.yticks(np.arange(len(MATRIX_FIELDS)), MATRIX_FIELDS, fontsize= 5)
    plt.savefig(MATRIX_NAME)
    plt.close()


def create_histogram_plots(sane_list, temp_cmb):

    n_bins = N_BINS_LIN
    sep_data  = transpose_unyt(sorted(gclumps(PLOT_FIELDS[0], sane_list)))
    sep_values = np.hstack([unyt_quantity(0.0, sep_data.units), \
                        np.linspace(sep_data[0], sep_data[-1], n_bins + 1), \
                        np.max(gclumps(PLOT_FIELDS[2], sane_list))])
    plot_values = np.hstack([np.min(gclumps(PLOT_FIELDS[0], sane_list)), \
                        np.linspace(sep_data[0], sep_data[-1], n_bins + 1), \
                        np.max(gclumps(PLOT_FIELDS[2], sane_list))])

    if ('temperature' in PLOT_FIELDS[0]):
        plt_cmb = (temp_cmb - sep_data[0]) * (n_bins / (sep_data[-1] - sep_data[0])) + 1
    n_clumps_max = 0
    n_clumps = np.zeros((3, n_bins + 3), dtype=int)

    for i, field in enumerate(PLOT_FIELDS):
        data  = transpose_unyt(sorted(gclumps(field, sane_list)))
        for datum in data:
            for j in range (n_bins + 3):
                if (datum > sep_values[j] and datum < sep_values[j+1]):
                    n_clumps[i][j] += 1
                    if (n_clumps[i][j] > n_clumps_max):
                        n_clumps_max = n_clumps[i][j]
                    break

    for i, field in enumerate(PLOT_FIELDS):
        plt.figure()
        barlist = plt.bar(np.arange(n_bins + 3), n_clumps[i], width=1.0, align='edge')
        barlist[np.argmax(n_clumps[i])].set_color('red')
        val_median = np.median(gclumps(PLOT_FIELDS[i], sane_list))
        if (val_median < 1.0 or (sep_data[-1] - sep_data[0]) > 100.0):
            plt.xticks(np.arange(n_bins + 3), [f'{num:.1e}' for num in plot_values.value], rotation=20, fontsize=8)
        else :
            plt.xticks(np.arange(n_bins + 3), [f'{num:.0f}' for num in plot_values.value])
        plt.xlabel(plot_values.units, horizontalalignment='center', y=1.0)
        plt.ylabel("N")
        plt.ylim(0, (n_clumps_max // 5) * 5 + 5)
        plt_median = (val_median - sep_data[0]) * (n_bins / (sep_data[-1] - sep_data[0])) + 1
        plt.axvline(x= plt_median, color='red')
        if ('temperature' in field):
            plt.axvline(x= plt_cmb, color='green', linestyle='--')
        ax = plt.gca()
        ax.get_xticklabels()[-2].set_color('red')
        ax.get_xticklabels()[-2].set_weight('bold')
        ax.get_xticklabels()[1].set_color('red')
        ax.get_xticklabels()[1].set_weight('bold')
        plt.savefig(f"{PLOT_NAME}_{field}.png")
        plt.close()


df = yt.load(CLUMP_FILE_NAME)
temp_cmb = CMB_Z0 * (df.current_redshift + 1.0)
clump_list = df.leaves
sane_list = sanatise_clumps(clump_list)
#create_matrix_plot(sane_list)
create_histogram_plots(sane_list, temp_cmb)
#create_tree_plot(df)
