import numpy as np
import os
import yt
import h5py
yt.enable_parallelism()
import ytree

from yt.frontends.enzo.data_structures import EnzoDataset
from yt.data_objects.level_sets.api import Clump, find_clumps, add_validator

from ytree.analysis import AnalysisPipeline

from yt.extensions.sam_mods.misc import return_sphere, align_sphere, transpose_unyt, modify_grackle
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.tree_analysis_operations import yt_dataset, node_profile, garbage_collect, delattrs
from yt.extensions.sam_mods import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc

from grid_figure import GridFigure



BIN_DENSITY = 20


def all_profiles(ap, outdir='.'):
    
    x_bin_field = ("index", "radius")
    profile_fields = [("gas", "cell_mass")]
    weight_field = None
    pkwargs = {"accumulation": False, "bin_density": BIN_DENSITY}
    fpath = os.path.join(os.getcwd(), outdir)
    
    _fields = {"velocity_magnitude": {"units": "km/s", "log": False},
               "velocity_spherical_radius": {"units": "km/s", "log": False},
               "velocity_spherical_theta": {"units": "km/s", "log": False},
               "velocity_spherical_phi": {"units": "km/s", "log": False},
               "tangential_velocity_magnitude": {"units": "km/s", "log": False},
               "sound_speed": {"units": "km/s", "log": False},
               "vortical_time": {"units": "yr", "log": True},
               "cooling_time": {"units": "yr", "log": True},
               "dynamical_time": {"units": "yr", "log": True},
               "total_dynamical_time": {"units": "yr", "log": True}}
    my_fields = [("gas", field) for field in _fields]
    for field in my_fields:
        my_kwargs = pkwargs.copy()
        my_kwargs["logs"] = {x_bin_field: True, field: _fields[field[1]]["log"]}
        my_kwargs["units"] = {x_bin_field: "pc", field: _fields[field[1]]["units"]}
        ap.add_operation(node_profile, [x_bin_field, field], profile_fields, weight_field,
                         profile_kwargs=my_kwargs, output_dir=fpath)

    my_kwargs = pkwargs.copy()
    my_kwargs["logs"] = {x_bin_field: True}
    my_kwargs["units"] = {x_bin_field: "pc"}
    for field in my_fields:
        my_kwargs["logs"][field] = _fields[field[1]]["log"]
        my_kwargs["units"][field] = _fields[field[1]]["units"]
    ap.add_operation(node_profile, [x_bin_field], my_fields, ("gas", "cell_mass"),
                     profile_kwargs=my_kwargs, output_dir=fpath)

    # 1D mass vs. radius profiles to get circular velocity
    mass_fields = [("gas", "cell_mass"),
                   ("gas", "dark_matter_mass"),
                   ("gas", "matter_mass")]
    my_kwargs = pkwargs.copy()
    my_kwargs["logs"] = {x_bin_field: True}
    my_kwargs["units"] = {x_bin_field: "pc"}
    ap.add_operation(node_profile, [x_bin_field], mass_fields, weight_field,
                     profile_kwargs=my_kwargs, output_dir=fpath)




if __name__ == "__main__":

    output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"

    es = yt.load(sim_path)
    a = ytree.load(tree_path)
    if "icom_gas2_position_x" in a.field_list:
        a.add_vector_field("icom_gas2_position")

    ap = AnalysisPipeline()
    ap.add_operation(yt_dataset, data_dir, add_fields= True)
    #ap.add_operation(yt_dataset, data_dir, es)
    ap.add_operation(modify_grackle)
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)
   
    ap.add_recipe(all_profiles, outdir='Profiles/Velocity_and_timescale_profiles_nocmb')
   
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
