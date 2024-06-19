import numpy as np
import yt

from yt.extensions.sam_mods.projection_image import multi_image

filenames = ["DD0102_temperature_slice.h5", "DD0110_temperature_slice.h5",
             "DD0104_density_proj.h5", "DD0130_density_proj.h5"]

panels = []
panels.append({"filename": filenames[0], "quantity": ("gas", "temperature"),
               "range": ["min", "max"], "cmap": "turbo", "label": "T [K]"})
panels.append({"filename": filenames[1], "quantity": ("gas", "temperature"),
               "range": ["min", "max"], "cmap": "turbo", "label": "T [K]"})
panels.append({"filename": filenames[2], "quantity": ("gas", "density"),
               "range": ["min", "max"], "cmap": "algae", "label": "$\\rho$ [g/cm$^{3}$]"})
panels.append({"filename": filenames[3], "quantity": ("gas", "density"),
               "range": ["min", "max"], "cmap": "algae", "label": "$\\rho$ [g/cm$^{3}$]"})

my_fig = multi_image(panels, "projections_noicom.pdf", figsize=(8, 8), fontsize=16, dpi=200,
                     n_columns=2, bg_color="white", text_color="black",
                     bottom_buffer=0.1, top_buffer=0.0, left_buffer=0.12, right_buffer=0.12)
my_fig = multi_image(panels, "projections_noicom.png", figsize=(8, 8), fontsize=16, dpi=200,
                     n_columns=2, bg_color="white", text_color="black",
                     bottom_buffer=0.1, top_buffer=0.0, left_buffer=0.12, right_buffer=0.12)
