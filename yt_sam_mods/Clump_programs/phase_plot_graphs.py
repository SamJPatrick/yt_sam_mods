import os
import yt
import sys
import numpy as np

from yt.extensions.sam_mods.graph_funcs import *
from yt.extensions.sam_mods.plots import *
from grid_figure import GridFigure

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import unyt
from unyt import unyt_quantity, unyt_array
from unyt import Myr, yr, pc

'''
CMB_Z0 = unyt_quantity(2.725, 'K')
INFILE = "Profiles/DD0295_profile_temperature_number_density.h5"
OUTDIR = "Profiles"
X_AXIS = 'number_density'
Y_AXIS = 'temperature'
df = yt.load(INFILE)
used = df.data[('data', 'used')].astype(bool).value
x_axis = df.data[('data', X_AXIS)][used]
y_axis = df.data[('data', Y_AXIS)][used]
'''

PROF_DIR = "Profiles"
OUT_DIR = '.'

xfield = 'number_density'
yfield = 'temperture'
input_filename = os.path.join(PROF_DIR, "DD0179_profile_number_density_temperature.h5")
output_filename = os.path.join(OUT_DIR, "t_vs_n_vs_M_phase_plot_ccsn.png")
my_fig = GridFigure(1, 1, figsize=(6, 4.5),
                    top_buffer = 0.03, bottom_buffer = 0.13,
                    left_buffer = 0.16, right_buffer = 0.19)

my_axes = my_fig[0]
xscale = 'log'
yscale = 'log'

field = 'cell_mass'
units = 'Msun'
cmap = plt.get_cmap('dusk')
#cmap = yt_colormaps['turbo']
clabel = "M [M$_{\odot}$]"

xlim = (1e-4, 1e9)
xmajor = np.logspace(-4, 9, 14)
#xminor = np.logspace(-4, 9, 16)
xminor = None
xlabel = "n [cm$^{-3}$]"

ylim = (1e1, 1e4)
ymajor = np.logspace(1, 4, 4)
#yminor = np.logspace(1, 3, 30)
yminor = None
ylabel = "T [K]"

make_phase_plot(
    my_fig, my_axes, input_filename,
    field, units, cmap, clabel,
    xlim, xmajor, xminor, xscale, xlabel,
    ylim, ymajor, yminor, yscale, ylabel,
    output_filename)
