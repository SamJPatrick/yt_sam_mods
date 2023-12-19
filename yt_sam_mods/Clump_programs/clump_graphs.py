import yt
import numpy as np
from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt



CLUMP_FILE_NAME = "Clumps/DD0295_clump_info.h5"
PLOT_NAME = "clump_correlations.png"
#MATRIX_FIELDS = ['distance', 'alpha', 'ratio', 'volume', 'cell_mass', 'jeans_mass',
#                 'max_metallicity', 'min_temperature', 'mean_number_density']
MATRIX_FIELDS = ['volume', 'min_temperature', 'mean_temperature', 'max_temperature']
#PLOT_FIELDS = ['mean_temperature', 'min_temperature', 'max_temperature']
PLOT_FIELDS = ['mean_number_density', 'min_number_density', 'max_number_density']
N_BINS_LOG = 20
N_BINS_LIN = 10

ERR_TOL = 0.2
TEMP_THRESH = unyt_quantity(200.0, 'K')
CMB_TOL = unyt_quantity(10.0, 'K')
VOL_THRESH = unyt_quantity(1e-6, 'pc**3')
CMB_Z0 = unyt_quantity(2.725, 'K')


df = yt.load(CLUMP_FILE_NAME)
temp_cmb = CMB_Z0 * (df.current_redshift + 1.0)
order_mask = [np.argwhere(df.data[('clump', 'cell_mass')] == mass).item() \
              for mass in sorted(df.data[('clump', 'cell_mass')])]
sanity_mask = (df.data[('clump', 'volume')][order_mask] > VOL_THRESH) & \
    (df.data[('clump', 'max_number_density')][order_mask] < 10 * (1 + ERR_TOL) * df.data[('clump', 'mean_number_density')][order_mask]) & \
    (df.data[('clump', 'mean_temperature')][order_mask] < TEMP_THRESH) & \
    (df.data[('clump', 'min_temperature')][order_mask] < temp_cmb + CMB_TOL)
duff_frac = (len(np.nonzero(sanity_mask)) / len(sanity_mask)) * 100
print(f"Percentage of clumps eliminated is {duff_frac:.2f}%")

coms = df.data[('clump', 'com')][order_mask][sanity_mask]
mean_com = unyt_array([np.mean(coord) for coord in list(zip(*coms.value.tolist()))], 'pc')
distance = [np.linalg.norm(unyt_array(com, 'pc') - mean_com) for com in coms]

corr_matrix = np.zeros((len(MATRIX_FIELDS), len(MATRIX_FIELDS)))
for i, field in enumerate(MATRIX_FIELDS):
    if (field == 'distance'):
        arr1 = distance
    else :
        arr1 = df.data[('clump', field)][order_mask][sanity_mask]
    for j, field in enumerate(MATRIX_FIELDS):
        if (field == 'distance'):
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

n_bins = N_BINS_LIN
sep_data  = transpose_unyt(sorted(df.data[('clump', PLOT_FIELDS[0])][order_mask][sanity_mask]))
sep_values = np.hstack([unyt_quantity(0.0, sep_data.units), np.linspace(sep_data[0], sep_data[-1], n_bins + 1),
                        np.max(transpose_unyt(df.data[('clump', PLOT_FIELDS[2])]))]) * unyt_quantity(1.0, sep_data.units)
print(sep_values)
if ('temperature' in PLOT_FIELDS[0]):
    plt_cmb = (temp_cmb - sep_data[0]) * (n_bins / (sep_data[-1] - sep_data[0])) + 1
n_clumps_max = 0
n_clumps = np.zeros((3, n_bins + 3), dtype=int)

for i, field in enumerate(PLOT_FIELDS):
    data  = transpose_unyt(sorted(df.data[('clump', field)][order_mask][sanity_mask]))
    for datum in data:
        for j in range (n_bins + 3):
            if (datum > sep_values[j] and datum <= sep_values[j+1]):
                n_clumps[i][j] += 1
                if (n_clumps[i][j] > n_clumps_max):
                    n_clumps_max = n_clumps[i][j]
                break

for i, field in enumerate(PLOT_FIELDS):
    plt.figure()
    barlist = plt.bar(np.arange(n_bins + 3), n_clumps[i], width=1.0, align='edge')
    barlist[np.argmax(n_clumps[i])].set_color('red')
    val_median = np.median(df.data[('clump', PLOT_FIELDS[i])])
    if (val_median < 1.0 or (sep_data[-1] - sep_data[0]) > 100.0):
        plt.xticks(np.arange(n_bins + 3), [f'{num:.1e}' for num in sep_values.value], rotation=20, fontsize=8)
    else :
        plt.xticks(np.arange(n_bins + 3), [f'{num:.0f}' for num in sep_values.value])
    plt.xlabel(sep_values.units, horizontalalignment='center', y=1.0)
    plt.ylabel("N")
    plt.ylim(0, (n_clumps_max // 5) * 5 + 5)
    plt_median = (val_median - sep_data[0]) * (n_bins / (sep_data[-1] - sep_data[0])) + 1
    plt.axvline(x= plt_median, color='red')
    if ('temperature' in field):
        plt.axvline(x= plt_cmb, color='green', linestyle='--')
    plt.savefig(f"clump_{field}.png")
    plt.close()
