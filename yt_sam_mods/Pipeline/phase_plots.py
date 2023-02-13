import numpy as np
import os
import yt
import h5py
yt.enable_parallelism()
import ytree

from yt.frontends.enzo.data_structures import EnzoDataset
from yt.data_objects.level_sets.api import Clump, find_clumps, add_validator

from ytree.analysis import AnalysisPipeline

from yt.extensions.sam_mods.misc import return_sphere, align_sphere, transpose_unyt
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.tree_analysis_operations import yt_dataset, garbage_collect, delattrs
from yt.extensions.sam_mods import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc




TINY_RADIUS = 1e-3



def metallicity_distribution(node, outdir='.', nbins=120):

    ds = node.ds
    sphere = node.sphere

    metal_dstrb = my_profile(sphere, ('gas', 'metallicity3'), ('gas', 'cell_volume'),\
                             extrema={('gas', 'metallicity3'): (1e-17, 1e0), ('gas', 'cell_volume'): unyt_array((1e7, 1e1), 'pc**3')}, \
                             n_bins= nbins, weight_field= None)
    plt = yt.ProfilePlot.from_profiles(metal_dstrb)
    plt.set_unit(('gas', 'cell_volume'), 'pc**3')
    plt.annotate_text(xpos= 1e-14, ypos= 1e3, text= 'z={:.2f}'.format(ds.current_redshift)) #coord_system= 'axis'

    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_metallicity_distribution")
    plt.save(name= fpath, suffix= ".png")
    

def metal_phase_plots(node, outdir='.'):

    ds = node.ds
    sphere = node.sphere

    plt = yt.PhasePlot(sphere, ("gas", "metallicity3_min7"), ("gas", "density"), [("gas", "cell_mass")], weight_field=None)
    plt.set_xlim(1e-7, 1e0)
    plt.set_unit(('gas', 'density'), 'g * cm**(-3)')
    plt.set_ylim(1e-29, 1e-15)
    plt.set_unit(('gas', 'cell_mass'), 'Msun')
    plt.set_zlim(('gas', 'cell_mass'), 1e-2, 1e2)
    decorate_plot(node, plt)
    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_metallicity_density_mass.png")
    plt.save(name= fpath)
    
    plt = yt.PhasePlot(sphere, 'radius', ('gas', 'metallicity3_min7'), [('gas', 'density')], weight_field=None)
    plt.set_unit('radius', 'pc')
    plt.set_log('radius', True)
    plt.set_xlim(1e-2, 5e2)
    plt.set_ylim(1e-7, 1e0)
    plt.set_unit(('gas', 'density'), 'g * cm**(-3)')
    plt.set_zlim(('gas', 'density'), 1e-27, 1e-15)
    decorate_plot(node, plt)
    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_radius_metallicity_density.png")
    plt.save(name= fpath)

    plt = yt.PhasePlot(sphere, ('gas', 'metallicity3_min7'), ('gas', 'density'), [('gas', 'cooling_ratio')], weight_field=('gas', 'cell_mass'))
    plt.set_xlim(1e-7, 1e0)
    plt.set_ylim(1e-29, 1e-15)
    plt.set_zlim(('gas', 'cooling_ratio'), 1e-2, 1e3)
    plt.set_log(('gas', 'cooling_ratio'), True)
    decorate_plot(node, plt)
    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_metallicity_density_ratio.png")
    plt.save(name= fpath)
    


    
if __name__ == "__main__":

    #output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
    #sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
    #data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    #tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"

    output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/minihalo_analysis"
    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/simulation.h5"
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust"
    tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"

    es = yt.load(sim_path)
    a = ytree.load(tree_path)
    if "icom_gas2_position_x" in a.field_list:
        a.add_vector_field("icom_gas2_position")
    
    ap = AnalysisPipeline()
    ap.add_operation(yt_dataset, data_dir, add_fields=True)
    #ap.add_operation(yt_dataset, data_dir, es)
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)
    
    ap.add_operation(metal_phase_plots, outdir= 'Phase_plots')
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
