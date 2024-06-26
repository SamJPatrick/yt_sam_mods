import numpy as np
import yt

from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.projection_image import multi_image


DENS_LIMS = unyt_array([1e-27, 1e-21], 'g/cm**3')
METAL_LIMS = unyt_array([1e-7, 1e0], 'Zsun')
TEMP_LIMS = unyt_array([1e1, 1e8], 'K')
TEMP_LIMS_SLICE = unyt_array([1e1, 1e6], 'K')
#RANGE = unyt_array([-150, 150], 'pc')
RANGE = [-150, 150]

'''
filenames = ['DD0150_density_proj.h5', 'DD0160_density_proj.h5', 'DD0175_density_proj.h5',
             'DD0195_density_proj.h5', 'DD0225_density_proj.h5', 'DD0254_density_proj.h5']

panels = []
panels.append({"filename": filenames[0], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "turbo", "label": "$\\rho$ [g/cm$^{3}$]"})
panels.append({"filename": filenames[1], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "turbo", "label": "$\\rho$ [g/cm$^{3}$]"})
panels.append({"filename": filenames[2], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "turbo", "label": "$\\rho$ [g/cm$^{3}$]"})
panels.append({"filename": filenames[3], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "turbo", "label": "$\\rho$ [g/cm$^{3}$]"})
panels.append({"filename": filenames[4], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "turbo", "label": "$\\rho$ [g/cm$^{3}$]"})
panels.append({"filename": filenames[5], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "turbo", "label": "$\\rho$ [g/cm$^{3}$]"})

my_fig = multi_image(panels, "paper_projections_pisn.pdf", figsize=(8, 8), fontsize=16, dpi=200,
                     n_columns=3, bg_color="white", text_color="black", show_cbar_label= False,
                     tick_range = RANGE, tick_label = 'pc', n_ticks = 7, cbar_orientation='horizontal',
                     bottom_buffer=0.1, top_buffer=0.0, left_buffer=0.12, right_buffer=0.12)
my_fig = multi_image(panels, "paper_projections_pisn.png", figsize=(8, 8), fontsize=16, dpi=200,
                     n_columns=3, bg_color="white", text_color="black", show_cbar_label= False,
                     tick_range = RANGE, tick_label = 'pc', n_ticks = 7, cbar_orientation='horizontal',
                     bottom_buffer=0.1, top_buffer=0.0, left_buffer=0.12, right_buffer=0.12)
'''

'''
filenames = ['DD0115_density_slice.h5', 'DD0115_temperature_slice.h5', 'DD0115_metallicity3_slice.h5',
             'DD0161_density_slice.h5', 'DD0161_temperature_slice.h5', 'DD0161_metallicity3_slice.h5',
             'DD0150_density_slice.h5', 'DD0150_temperature_slice.h5', 'DD0150_metallicity3_slice.h5']

panels = []
panels.append({"filename": filenames[0], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "turbo", "label": "$\\rho$ [g/cm$^{3}$]"})
panels.append({"filename": filenames[1], "quantity": ("gas", "temperature"),
               "range": [*TEMP_LIMS], "cmap": "inferno", "label": "T [K]"})
panels.append({"filename": filenames[2], "quantity": ("gas", "metallicity3"),
               "range": [*METAL_LIMS], "cmap": "cool", "label": "Z [Z$_{\odot}$]"})
panels.append({"filename": filenames[3], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "turbo", "label": "$\\rho$ [g/cm$^{3}$]"})
panels.append({"filename": filenames[4], "quantity": ("gas", "temperature"),
               "range": [*TEMP_LIMS], "cmap": "inferno", "label": "T [K]"})
panels.append({"filename": filenames[5], "quantity": ("gas", "metallicity3"),
               "range": [*METAL_LIMS], "cmap": "cool", "label": "Z [Z$_{\odot}$]"})
panels.append({"filename": filenames[6], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "turbo", "label": "$\\rho$ [g/cm$^{3}$]"})
panels.append({"filename": filenames[7], "quantity": ("gas", "temperature"),
               "range": [*TEMP_LIMS], "cmap": "inferno", "label": "T [K]"})
panels.append({"filename": filenames[8], "quantity": ("gas", "metallicity3"),
               "range": [*METAL_LIMS], "cmap": "cool", "label": "Z [Z$_{\odot}$]"})

my_fig = multi_image(panels, "paper_slices.pdf", figsize=(8, 8), fontsize=16, dpi=200,
                     n_columns=3, bg_color="white", text_color="black", show_cbar_label= False,
                     tick_range = RANGE, tick_label = 'pc', n_ticks = 7, cbar_orientation='horizontal',
                     bottom_buffer=0.1, top_buffer=0.0, left_buffer=0.12, right_buffer=0.12)
my_fig = multi_image(panels, "paper_slices.png", figsize=(8, 8), fontsize=16, dpi=200,
                     n_columns=3, bg_color="white", text_color="black", show_cbar_label= False,
                     tick_range = RANGE, tick_label = 'pc', n_ticks = 7, cbar_orientation='horizontal',
                     bottom_buffer=0.1, top_buffer=0.0, left_buffer=0.12, right_buffer=0.12)
'''



filenames = ["DD0102_temperature_slice.h5", "DD0110_temperature_slice.h5",
             "DD0104_density_proj.h5", "DD0130_density_proj.h5"]

panels = []
panels.append({"filename": filenames[0], "quantity": ("gas", "temperature"),
               "range": [*TEMP_LIMS_SLICE], "cmap": "turbo", "label": "T [K]"})
panels.append({"filename": filenames[1], "quantity": ("gas", "temperature"),
               "range": [*TEMP_LIMS_SLICE], "cmap": "turbo", "label": "T [K]"})
panels.append({"filename": filenames[2], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "algae", "label": "$\\rho$ [g/cm$^{3}$]"})
panels.append({"filename": filenames[3], "quantity": ("gas", "density"),
               "range": [*DENS_LIMS], "cmap": "algae", "label": "$\\rho$ [g/cm$^{3}$]"})

my_fig = multi_image(panels, "paper_multiplot.pdf", figsize=(8, 8), fontsize=16, dpi=200,
                     n_columns=2, bg_color="white", text_color="black", show_cbar_label= False,
                     tick_range = RANGE, tick_label = 'pc', n_ticks = 7, cbar_orientation='vertical',
                     bottom_buffer=0.1, top_buffer=0.0, left_buffer=0.12, right_buffer=0.12)
my_fig = multi_image(panels, "paper_multiplot.png", figsize=(8, 8), fontsize=16, dpi=200,
                     n_columns=2, bg_color="white", text_color="black", show_cbar_label= False,
                     tick_range = RANGE, tick_label = 'pc', n_ticks = 7, cbar_orientation='vertical',
                     bottom_buffer=0.1, top_buffer=0.0, left_buffer=0.12, right_buffer=0.12)

