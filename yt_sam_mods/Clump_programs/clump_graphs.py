import yt
import numpy as np
from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt



CLUMP_FILE_NAME = "DD0295_clump_0.h5"
PLOT_NAME = "clump_correlations_dens.png"

PLOT_FIELDS = {'max_metallicity': False, 'mean_metallicity': False, 'min_temperature': False, 'mean_temperature': False, \
               'max_number_density': True, 'mean_number_density': True, 'stability': True, 'cell_mass': True, 'volume': True}

#MATRIX_FIELDS = ['mean_temperature', 'max_temperature', 'min_temperature', 'volume']
#MATRIX_FIELDS = ['mean_metallicity', 'max_metallicity', 'min_metallicity', 'cell_mass']
MATRIX_FIELDS = ['mean_number_density', 'max_number_density', 'min_number_density', 'cell_mass']

#MATRIX_FIELDS = ['mean_number_density', 'mean_metallicity', 'mean_temperature']
#MATRIX_FIELDS = ['max_number_density', 'max_metallicity', 'min_temperature']
#MATRIX_FIELDS = ['distance', 'alpha', 'ratio', 'volume', 'a_max']
#MATRIX_FIELDS = ['distance', 'alpha', 'ratio', 'volume', 'cell_mass', 'stability', 'max_metallicity', \
#                 'mean_metallicity', 'min_metallicity', 'max_temperature', 'mean_temperature', 'min_temperature', \
#                 'max_number_density', 'mean_number_density', 'min_number_density']

N_BINS_LOG = 20
N_BINS_LIN = 10



df = yt.load(CLUMP_FILE_NAME)
order_mask = [np.argwhere(df.data[('clump', 'cell_mass')] == mass).item() \
              for mass in sorted(df.data[('clump', 'cell_mass')])]
sanity_mask = (df.data[('clump', 'mean_temperature')][order_mask].to('K').value < 40) & \
    (df.data[('clump', 'max_number_density')][order_mask] < 10 * df.data[('clump', 'mean_number_density')][order_mask]) & \
    (df.data[('clump', 'volume')][order_mask].to('pc**3').value > 1e-6)
duff_frac = (len(np.nonzero(sanity_mask)) / len(sanity_mask)) * 100
print(f"Percentage of clumps eliminated is {duff_frac:.2f}%")

temp_diff = df.data[('clump', 'max_temperature')][order_mask][sanity_mask] - \
    df.data[('clump', 'min_temperature')][order_mask][sanity_mask]
stability = df.data[('clump', 'jeans_mass')][order_mask][sanity_mask] / \
    df.data[('clump', 'cell_mass')][order_mask][sanity_mask]
coms = df.data[('clump', 'com')][order_mask][sanity_mask]
mean_com = unyt_array([np.mean(coord) for coord in list(zip(*coms.value.tolist()))], 'pc')
distance = [np.linalg.norm(unyt_array(com, 'pc') - mean_com) for com in coms]

corr_matrix = np.zeros((len(MATRIX_FIELDS), len(MATRIX_FIELDS)))
for i, field in enumerate(MATRIX_FIELDS):
    if (field == 'temp_diff'):
        arr1 = temp_diff
    elif (field == 'stability'):
        arr1 = stability
    elif (field == 'distance'):
        arr1 = distance
    else :
        arr1 = df.data[('clump', field)][order_mask][sanity_mask]
    for j, field in enumerate(MATRIX_FIELDS):
        if (field == 'temp_diff'):
            arr2 = temp_diff
        elif (field == 'stability'):
            arr2 = stability
        elif (field == 'distance'):
            arr2 = distance
        else :
            arr2 = df.data[('clump', field)][order_mask][sanity_mask]
        corr_matrix[i][j] = np.corrcoef(arr1, arr2)[0][1]

plt.figure()
plt.matshow(np.abs(corr_matrix))
for i in range (len(MATRIX_FIELDS)):
    for j in range (len(MATRIX_FIELDS)):
        plt.text(i, j, f'{corr_matrix[i][j]:.3f}', ha='center', va='center', color='black', fontsize=5)
plt.xticks(np.arange(len(MATRIX_FIELDS)), MATRIX_FIELDS, fontsize= 5, rotation='vertical')
plt.yticks(np.arange(len(MATRIX_FIELDS)), MATRIX_FIELDS, fontsize= 5)
plt.savefig(PLOT_NAME)
plt.close()

'''
x_axis_metal = None
x_axis_temps = None
x_axis_number = None
for field in PLOT_FIELDS.keys():

    if (field == 'stability'):
        data = transpose_unyt(sorted(stability[order_mask][sanity_mask]))
    elif (field == 'temp_diff'):
        data = transpose_unyt(sorted(temp_diff[order_mask][sanity_mask]))
    else :
        data  = transpose_unyt(sorted(df.data[('clump', field)][order_mask][sanity_mask]))

    if (PLOT_FIELDS[field]):
        n_bins = N_BINS_LOG
    else :
        n_bins = N_BINS_LIN
    if (PLOT_FIELDS[field]):
        x_axis = unyt_array(np.logspace(np.log10(data[0]), np.log10(data[-1]), n_bins), \
                            data.units)
    else :
        x_axis = np.linspace(data[0], data[-1], n_bins)
        
    if ('temperature' in field):
        if (x_axis_temps is None):
            x_axis_temps = x_axis
        x_axis = x_axis_temps   
    if ('metal' in field):
        if (x_axis_metal is None):
            x_axis_metal = x_axis
        x_axis = x_axis_metal
    if ('number' in field):
        if (x_axis_number is None):
            x_axis_number = x_axis
        x_axis = x_axis_number
        
    n_clumps = np.zeros(n_bins)
    for datum in data:
        for i in range (n_bins):
            if (i == (n_bins - 1)):
                n_clumps[-1] +=1
                break
            if (datum >= x_axis[i] and datum <= x_axis[i+1]):
                n_clumps[i] += 1
                break

    plt.figure()
    barlist = plt.bar(np.arange(n_bins), n_clumps, width=1.0, align='edge')
    barlist[np.argmax(n_clumps)].set_color('red')
    if ('metal' in field or PLOT_FIELDS[field]):
        plt.xticks(np.arange(n_bins), [f'{x:.2E}' for x in x_axis.value], fontsize=5, rotation='vertical')
    elif ('temperature' in field):
        plt.xticks(np.arange(n_bins), [f'{x:.1f}' for x in x_axis.value], fontsize=7, rotation=30)
    else :
        plt.xticks(np.arange(n_bins), x_axis.value, fontsize=7)
    if (('temperature' in field) or ('metal' in field)):
        plt.ylim(0, 15)
    elif ('density' in field):
        plt.ylim(0, 20)
    plt.xlabel(x_axis.units, horizontalalignment='center', y=1.0)
    plt.ylabel("N")
    #print(field)
    #print(x_axis)
    #print(n_clumps)
    if (PLOT_FIELDS[field]):
    #    x_ave = np.exp(np.average(np.log(x_axis), weights= n_clumps)) * (N_BINS/ np.log10(x_axis[-1])) - 1
    #    print(x_ave)
    #    plt.axvline(x= x_ave, color= 'green')
        pass
    else :
        plt.axvline(x= np.average(x_axis, weights= n_clumps) * (n_bins/x_axis[-1]) - 1, color='red')
    plt.savefig(f"clump_{field}.png")
    plt.close()
'''
