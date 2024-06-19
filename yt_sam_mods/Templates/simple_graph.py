import yt
import numpy as np
from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt



FILE_NAME = "DD0098_profile_weight_field_cell_mass.h5"

df = yt.load(FILE_NAME)
x_data = df.data[('data', 'radius')]
y_data  = df.data[('data', 'density')]
plt.figure()
plt.plot(x_data, y_data)
plt.xlabel(r'$r (pc)$')
plt.ylabel(r"$\rho (g cm^{-3})$")
plt.xscale('log')
plt.yscale('log')
plt.savefig("pop3_pisn_profile.png")
plt.close()

