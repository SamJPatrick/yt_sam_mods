import numpy as np
import yt

from yt_p2p.projection_image import *

if __name__ == "__main__":

    ### to make projections and save the data
    # ds = yt.load('DD0046/DD0046')

    # p = yt.ProjectionPlot(ds, "x", [("gas", "density"),
    #                                 ("gas", "temperature"),
    #                                 ("gas", "metallicity")],
    #                       weight_field=("gas", "density"))
    # p.frb.save_as_dataset(filename="projections/gas.h5")

    # p = yt.ParticleProjectionPlot(ds, "x", ("all", "particle_mass"))
    # p.frb.save_as_dataset(filename="projections/particles.h5")

    ### to load in projections and plot them
    
    gas_fn = "projections/gas.h5"
    dm_fn = "projections/particles.h5"

    panels = []
    panels.append({"filename": gas_fn, "quantity": ("gas", "density"),
                   "range": ["min", "max"],
                   "cmap": "algae", "label": "$\\rho$ [g/cm$^{3}$]"})
    panels.append({"filename": gas_fn, "quantity": ("gas", "temperature"),
                   "range": ["min", "max"],
                   "cmap": "gist_heat", "label": "T [K]"})
    panels.append({"filename": gas_fn, "quantity": ("gas", "metallicity"),
                   "range": ["min", "max"],
                   "cmap": "kamae", "label": "Z [Z$_{\\odot}$]"})
    panels.append({"filename": dm_fn, "quantity": ("gas", "particle_mass"),
                   "range": ["min", "max"],
                   "cmap": "turbo", "label": "$\\rho_{dm}$"})


    my_fig = multi_image(panels, "projections.pdf", figsize=(8, 8), fontsize=16, dpi=200,
                         n_columns=2, bg_color="white", text_color="black",
                         bottom_buffer=0.1, top_buffer=0.0,
                         left_buffer=0.12, right_buffer=0.12)
