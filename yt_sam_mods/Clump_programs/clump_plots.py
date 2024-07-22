import yt
import numpy as np
from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
from yt.extensions.sam_mods.fragplot import TreePlot

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


N_BINS_LOG = 10
N_BINS_LIN = 10

ERR_TOL = 0.2
N_CRIT = unyt_quantity(5.0, 'cm**(-3)')
N_CELLS = 8
VOL_MAX = unyt_quantity(8.0e3, 'pc**3')
Z_CRIT = unyt_quantity(1.0e-4, '')
CMB_Z0 = unyt_quantity(2.725, 'K')



gclumps = lambda field, alist: transpose_unyt([clump[('clump', field)] for clump in alist])


def plt_median(field, sub_list, sep_data, log=False):

    min_value = min(sep_data)
    max_value = max(sep_data)
    if (log == False):
        median = np.median(gclumps(field, sub_list))
        plt_val = N_BINS_LIN * (median - min_value) / (max_value - min_value) + 1
    else :
        log_mean = np.mean(np.log10(gclumps(field, sub_list)))
        plt_val = N_BINS_LOG * (log_mean - np.log10(min_value)) / \
            (np.log10(max_value) - np.log10(min_value)) + 1
    return plt_val


def label_densities(clump):

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

    cond = (clump[('clump', 'volume')] < VOL_MAX and clump[('clump', 'total_cells')] > N_CELLS and clump[('clump', 'mean_number_density')] > N_CRIT)
        
    if (cond):
        my_kwargs = {'label': label, 'fontsize': 16, 'shape': shape}
    else :
        my_kwargs = {'width': 0}
    return my_kwargs


def label_ids(clump):

    label = clump[('clump', 'clump_id')]
    if (clump[('clump', 'fragmentation_instability')] > 1.0):
        shape = 'square'
    else :
        shape = 'circle'
    my_kwargs = {'label': label, 'fontsize': 16, 'shape': shape}
    return my_kwargs






def create_tree_plot(df, tree_name, function= label_densities):

    p = TreePlot(df, dot_kwargs={'rankdir': 'TB', 'size': 18}, node_function= function)
    p.size_field = "mass"
    p.color_field = "max_metallicity"
    p.node_cmap = 'turbo'
    p.save(tree_name)


def sanatise_clumps(clump_list):

    order_mask = np.argsort([clump[('clump', 'cell_mass')] for clump in clump_list])

    vol_min = (np.array([clump[('clump', 'total_cells')] for clump in clump_list])[order_mask] > N_CELLS)
    vol_max = (gclumps('volume', clump_list)[order_mask] < VOL_MAX)
    vol_mask = vol_min & vol_max
    dens_mask = (gclumps('mean_number_density', clump_list)[order_mask] > N_CRIT)
    err_mask = (gclumps('max_number_density', clump_list)[order_mask] < \
                10 * (1 + ERR_TOL) * gclumps('mean_number_density', clump_list)[order_mask])
    sanity_mask = err_mask & vol_mask & dens_mask

    sane_list = np.array(clump_list)[order_mask][sanity_mask]
    good_frac = (sum(sanity_mask) / len(sanity_mask)) * 100
    bad_clumps = len(sanity_mask) - sum(sanity_mask)
    print(f"Percentage of clumps remaining is {good_frac:.2f}% with {bad_clumps} bad clumps caught")
    return sane_list


def create_matrix_plot(sane_list, matrix_fields, matrix_name):

    corr_matrix = np.zeros((len(matrix_fields), len(matrix_fields)))
    for i, field_a in enumerate(matrix_fields):
        for j, field_b in enumerate(matrix_fields):
            corr_matrix[i][j] = np.corrcoef(gclumps(field_a, sane_list), gclumps(field_b, sane_list))[0][1]

    plt.figure()
    plt.matshow(np.abs(corr_matrix))
    for i in range (len(matrix_fields)):
        for j in range (len(matrix_fields)):
            plt.text(i, j, f'{corr_matrix[i][j]:.3f}', ha='center', va='center', color='black', fontsize=5)
    plt.xticks(np.arange(len(matrix_fields)), matrix_fields, fontsize= 5, rotation='vertical')
    plt.yticks(np.arange(len(matrix_fields)), matrix_fields, fontsize= 5)
    plt.savefig(matrix_name)
    plt.close()


def create_histogram_plots(sane_list, plot_fields, histogram_name, temp_cmb, sep_frag= True):

    sep_data  = np.sort(gclumps(plot_fields[0], sane_list))
    if (sep_data[-1]/sep_data[0] < 100.0):
        n_bins = N_BINS_LIN
        plot_values = np.hstack([np.min(gclumps(plot_fields[1], sane_list)), \
                                 np.linspace(sep_data[0], sep_data[-1], n_bins + 1), \
                                 np.max(gclumps(plot_fields[2], sane_list))])
    else :
        n_bins = N_BINS_LOG
        plot_values = np.hstack([np.min(gclumps(plot_fields[1], sane_list)), \
                                 np.logspace(np.log10(sep_data[0].value), np.log10(sep_data[-1].value), n_bins + 1), \
                                 np.max(gclumps(plot_fields[2], sane_list))])

    if ('temperature' in plot_fields[0]):
        plt_cmb = (temp_cmb - sep_data[0]) * (n_bins / (sep_data[-1] - sep_data[0])) + 1

    if (sep_frag):
        stable_clumps = [clump for clump in sane_list if clump[('clump', 'fragmentation_instability')] <= 1.0]
        unstable_clumps = [clump for clump in sane_list if clump[('clump', 'fragmentation_instability')] > 1.0]
        clump_tuple = (stable_clumps, unstable_clumps)
    else :
        clump_tuple = (sane_list,)
    
    n_clumps_max = 0
    n_clumps = np.zeros((len(plot_fields), len(clump_tuple), n_bins + 2), dtype=int)
    for i, field in enumerate(plot_fields):
        for k, clump_type in enumerate(clump_tuple):
            data  = np.sort(gclumps(field, clump_type))
            n_clumps[i][k][0] = np.sum(data < plot_values[1])
            n_clumps[i][k][-1] = np.sum(data > plot_values[-2])
            for datum in data:
                for j in range (1, n_bins + 1):
                    if (datum >= plot_values[j] and datum <= plot_values[j+1]):
                        n_clumps[i][k][j] += 1
                        if (n_clumps[i][k][j] > n_clumps_max):
                            n_clumps_max = n_clumps[i][k][j]
                        break

    for i, field in enumerate(plot_fields):
        plt.figure()
        if (sep_frag):
            stable_barlist = plt.bar(np.arange(n_bins + 2) + 0.5, n_clumps[i][0], width=0.5, color='blue', align='edge')
            unstable_barlist = plt.bar(np.arange(n_bins + 2), n_clumps[i][1], width=0.5, color='red', align='edge')
        else :
            barlist = plt.bar(np.arange(n_bins + 2), n_clumps[i][0], width=1.0, color='blue', align='edge')
            #barlist[np.argmax(n_clumps[i][0])].set_color('red')
        if (sep_data[-1]/sep_data[0] > 100.0):
            log = True
        else :
            log = False
        if (log == True or sep_data[-1] < 1.0):
            plt.xticks(np.arange(n_bins + 3), [f'{num:.1e}' for num in plot_values.value], rotation=20, fontsize=8)
        else :
            plt.xticks(np.arange(n_bins + 3), [f'{num:.0f}' for num in plot_values.value])
        plt.xlabel(f'{plot_fields[0].split("_", 1)[1]} ({sep_data.units})', horizontalalignment='center', y=1.0)
        plt.ylabel("N")
        plt.ylim(0, int((n_clumps_max // 5) * 5 + 5))
        plt.axvline(x= plt_median(plot_fields[i], clump_tuple[0], sep_data, log= log), color='blue')
        if (sep_frag):
            plt.axvline(x= plt_median(plot_fields[i], clump_tuple[1], sep_data, log= log), color='red')
        if ('temperature' in field):
            plt.axvline(x= plt_cmb, color='green', linestyle='--')
        ax = plt.gca()
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.grid(True, which='both', axis='y')
        ax.get_xticklabels()[-2].set_color('red')
        ax.get_xticklabels()[-2].set_weight('bold')
        ax.get_xticklabels()[1].set_color('red')
        ax.get_xticklabels()[1].set_weight('bold')
        if (sep_frag):
            plt.savefig(f"{histogram_name}_{field}_sep.png")
        else :
            plt.savefig(f"{histogram_name}_{field}_nosep.png")
        plt.close()


def create_single_histogram(sane_list, plot_field, histogram_name, sep_frag= True):

    data  = np.sort(gclumps(plot_field, sane_list))
    if (data[-1]/data[0] < 100.0):
        n_bins = N_BINS_LIN
        plot_values = unyt_array(np.linspace(data[0], data[-1], n_bins + 1), data.units)
    else :
        n_bins = N_BINS_LOG
        plot_values = unyt_array(np.logspace(np.log10(data[0].value), np.log10(data[-1].value), n_bins + 1), data.units)

    if (sep_frag):
        stable_clumps = [clump for clump in sane_list if clump[('clump', 'fragmentation_instability')] <= 1.0]
        unstable_clumps = [clump for clump in sane_list if clump[('clump', 'fragmentation_instability')] > 1.0]
        clump_tuple = (stable_clumps, unstable_clumps)
    else :
        clump_tuple = (sane_list,)
    
    n_clumps_max = 0
    n_clumps = np.zeros((len(clump_tuple), n_bins), dtype=int)
    for k, clump_type in enumerate(clump_tuple):
        for datum in gclumps(plot_field, clump_type):
            for j in range (0, n_bins):
                if (datum >= plot_values[j] and datum <= plot_values[j+1]):
                    n_clumps[k][j] += 1
                    if (n_clumps[k][j] > n_clumps_max):
                        n_clumps_max = n_clumps[k][j]
                    break

    plt.figure()
    if (sep_frag):
        stable_barlist = plt.bar(np.arange(n_bins) + 0.5, n_clumps[0], width=0.5, color='blue', align='edge')
        unstable_barlist = plt.bar(np.arange(n_bins), n_clumps[1], width=0.5, color='red', align='edge')
    else :
        barlist = plt.bar(np.arange(n_bins), n_clumps[0], width=1.0, color='blue', align='edge')
        #barlist[np.argmax(n_clumps[0])].set_color('red')
    if (data[-1]/data[0] > 100.0):
        log = True
    else :
        log = False
    if (log == True or data[-1] < 1.0):
        plt.xticks(np.arange(n_bins + 1), [f'{num:.1e}' for num in plot_values.value], rotation=20, fontsize=8)
    else :
        plt.xticks(np.arange(n_bins + 1), [f'{num:.0f}' for num in plot_values.value])
    plt.xlabel(f'{plot_field} ({data.units})', horizontalalignment='center', y=1.0)
    plt.ylabel("N")
    plt.ylim(0, int((n_clumps_max // 5) * 5 + 5))
    plt.axvline(x= plt_median(plot_field, clump_tuple[0], data, log= log), color='blue')
    if (sep_frag):
        plt.axvline(x= plt_median(plot_field, clump_tuple[1], data, log= log), color='red')
    ax = plt.gca()
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.grid(True, which='both', axis='y')
    if (sep_frag):
        plt.savefig(f"{histogram_name}_{plot_field}_sep.png")
    else :
        plt.savefig(f"{histogram_name}_{plot_field}_nosep.png")
    plt.close()

        
def create_graph(sane_list, plot_fields, graph_name, temp_cmb):
    
    plt.figure()
    plt.scatter(gclumps(plot_fields[0], sane_list), np.abs(gclumps(plot_fields[1], sane_list)),
                c= gclumps(plot_fields[2], sane_list), cmap= 'turbo', vmin= 1e-3, vmax= 5e1, norm=matplotlib.colors.LogNorm())
    #plt.xlabel(f'{plot_fields[0]}')
    #plt.xlabel(r'$Mass (M_{\odot})$')
    plt.xlabel(r'$Volume ~(pc^{-3})$')
    #plt.xlabel(r'$\alpha_{\rho}$')
    #plt.ylabel(r'$\alpha_{\rho}$')
    #plt.ylabel(r'$\alpha_{Z}$')
    plt.ylabel(r'$Z/Z_{\odot}$')
    #plt.ylabel(r'$\tau_{dyn}/\tau_{cool}$')
    #plt.ylabel(r'$Temperature ~(K)$')
    #plt.ylim(0, 8)
    #plt.ylim(0, 60)
    #plt.ylim(1e-3, 5e1)
    #plt.ylim(30, 70)
    #plt.ylim(30, 800)
    plt.ylim(1e-7, 1e-3)
    #plt.xlim(0, 60)
    #plt.xlim(1e-3, 5e1)
    plt.xlim(5e-12, 1e1)
    plt.xscale('log')
    plt.yscale('log')
    #plt.colorbar(label= "r_factor")
    #plt.colorbar(label= r"$Volume ~(pc^{-3})$")
    plt.colorbar(label= r"$\tau_{cool}/\tau_{dyn}$")
    #plt.axhline(y= 1.0)
    #plt.axhline(y= temp_cmb, color='green', linestyle='--')
    plt.axhline(y= 10**(-5.3), color='green', linestyle='--')
    plt.axhline(y= 1.35e-5, color='red', linestyle='--')
    #plt.colorbar(label= r"$\alpha_{\rho}$")
    plt.savefig(f"{graph_name}.png")
    plt.close()
  


if __name__ == '__main__':
    

    clump_file_name = "Clumps/DD0232_clump_info_hn.h5"
    tree_name = "clump_tree_hn_232_ids.png"
    matrix_name = "correlations_intr_hn_232"
    histogram_name = "histo_hn_232"
    graph_name = "volume_vs_metal_vs_ratio_232_hn_leaf"
    sep_frag = False
        
    matrix_fields = ['alpha_density', 'alpha_metal', 'volume', 'cell_mass', 'r_factor',
                     'jeans_mass', 'fragmentation_instability',
                     'mean_metallicity', 'mean_temperature', 'mean_number_density', 'mean_H2_fraction']
    matrix_fields = ['mean_temperature',
                     'mean_number_density',
                     'mean_H2_fraction',
                     'min_metallicity', 'mean_metallicity', 'max_metallicity']
    histo_fields = ['temperature', 'number_density', 'metallicity']
    plot_fields = ['volume', 'mean_metallicity', 'fragmentation_instability']
    
    df = yt.load(clump_file_name)
    temp_cmb = CMB_Z0 * (df.current_redshift + 1.0)
    
    leaf_clumps = df.leaves
    sane_leaves = sanatise_clumps(leaf_clumps)
    all_clumps = list(df.tree)
    sane_clumps = sanatise_clumps(all_clumps)
    
    #create_tree_plot(df, tree_name, function= label_ids)
    #create_graph(sane_leaves, plot_fields, graph_name, temp_cmb)
    #create_matrix_plot(sane_clumps, matrix_fields, matrix_name + "_all.png")
    #create_matrix_plot(sane_leaves, matrix_fields, matrix_name + "_leaves.png")
    create_single_histogram(sane_leaves, 'cell_mass', "histo_hn_232", sep_frag= sep_frag)
    for field in histo_fields:
        to_plot = [f'mean_{field}', f'min_{field}', f'max_{field}']
        create_histogram_plots(sane_leaves, to_plot, histogram_name, temp_cmb, sep_frag= sep_frag)

