import statistics
import numpy as np
import os
import yt
import h5py
yt.enable_parallelism()
import ytree

from yt.frontends.enzo.data_structures import EnzoDataset
from yt.data_objects.level_sets.api import Clump, find_clumps, add_validator

from ytree.analysis import AnalysisPipeline

from yt.extensions.p2p.misc import return_sphere, align_sphere, transpose_unyt
from yt.extensions.p2p.profiles import my_profile
from yt.extensions.p2p.plots import *
from yt.extensions.p2p.tree_analysis_operations import yt_dataset, garbage_collect, delattrs
from yt.extensions.p2p import add_p2p_fields

from yt import derived_field
from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc


def bonnor_ebert_rate(node, outdir='.'):

    BE_PATH = "Profiles/Bonnor_ebert_old/be_time_mass_radius.h5"

    ds = node.ds
    sphere = node.sphere

    be_info = h5py.File(os.path.join(os.getcwd(), BE_PATH), 'r')
    be_times = unyt_array(list(be_info['time/array_data']), 'Myr')
    time = ds.current_time.in_units('Myr')
    i, time_diff = min(enumerate(np.abs(be_times - time)), key= lambda x: x[1])
    be_mass = unyt_quantity(be_info['mass/array_data'][i], 'Msun')
    be_radius = unyt_quantity(be_info['radius/array_data'][i], 'pc')

    be_sphere = ds.sphere(sphere.center, be_radius)
    sigma = np.std(be_sphere[('gas', 'velocity_magnitude')])
    cs = be_sphere.quantities.weighted_average_quantity(('gas', 'sound_speed'), ('gas', 'cell_mass'))
    tvort = be_sphere.quantities.weighted_average_quantity(('gas', 'vortical_time'), ('gas', 'cell_mass'))
    be_rate = (4 * be_mass * (sigma / cs)**2 / tvort).in_units('Msun/yr')

    dx = ds.index.get_smallest_dx()
    ds.add_field(('gas', 'mass_flux'), function= multiply_fields(ds, ('gas', 'density'), ('gas', 'radial_velocity')), \
                 sampling_type='local', units='kg*s**(-1)*m**(-2)', take_log= False)
    in_bin = (sphere['radius'] < be_radius) & (sphere['radius'] >= be_radius - dx)
    mass_flux = np.sum(sphere[('gas', 'mass_flux')][up_bin] * \
                       sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
    mdot = (4 * np.pi * mass_flux * be_radius**2).in_units('Msun/yr')

    # yes, the tuples have to be there because h5py is absolutely moronic >.>
    fpath = os.path.join(os.getcwd(), outdir, f"bonnor_ebert_rates.h5")
    if (os.path.exists(fpath)):
        rate_info = h5py.File(fpath, 'a')
        rate_info['be_rate'].resize((len(rate_info['be_rate']) + 1,))
        rate_info['be_rate'][-1] = be_rate
        rate_info['mdot'].resize((len(rate_info['mdot']) + 1,))
        rate_info['mdot'][-1] = mdot
        rate_info['time'].resize((len(rate_info['time']) + 1,))
        rate_info['time'][-1] = time
    else :
        rate_info = h5py.File(fpath, 'w')
        rate_info.create_dataset('be_rate', (1,), maxshape= (len(be_times),), dtype= 'f')
        rate_info['be_rate'][0] = be_rate
        rate_info.create_dataset('mdot', (1,), maxshape= (len(be_times),), dtype= 'f')
        rate_info['mdot'][0] = mdot
        rate_info.create_dataset('time', (1,), maxshape= (len(be_times),), dtype= 'f')
        rate_info['time'][0] = time




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
    ap.add_operation(yt_dataset, data_dir, es)
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)
    
    ap.add_operation(bonnor_ebert_rate, outdir='.')
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
