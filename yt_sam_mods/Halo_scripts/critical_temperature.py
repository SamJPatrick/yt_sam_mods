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




def critical_temperature(node, outdir='.'):

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
    vsigma = np.std(be_sphere[('gas', 'velocity_magnitude')])
    cs = be_sphere.quantities.weighted_average_quantity(('gas', 'sound_speed'), ('gas', 'cell_mass'))
    mu = be_sphere.quantities.weighted_average_quantity(('gas', 'mean_molecular_mass'), ('gas', 'cell_mass'))
    temperature = be_sphere.quantities.weighted_average_quantity(('gas', 'temperature'), ('gas', 'cell_mass'))
    tvort = be_sphere.quantities.weighted_average_quantity(('gas', 'vortical_time'), ('gas', 'cell_mass'))
    tcool = be_sphere.quantities.weighted_average_quantity(('gas', 'cooling_time'), ('gas', 'cell_mass'))

    dx = ds.index.get_smallest_dx()
    ds.add_field(('gas', 'mass_flux'), function= multiply_fields(ds, ('gas', 'density'), ('gas', 'radial_velocity')), \
                 sampling_type='local', units='kg*s**(-1)*m**(-2)', take_log= False)
    in_bin = (sphere['radius'] < be_radius) & (sphere['radius'] >= be_radius - dx)
    mass_flux = np.sum(sphere[('gas', 'mass_flux')][up_bin] * \
                       sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
    mdot = (4 * np.pi * mass_flux * be_radius**2).in_units('Msun/yr')
    vel_r = np.sum(sphere[('gas', 'radial_velocity')][up_bin] * \
                   sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
    density = np.sum(sphere[('gas', 'density')][up_bin] * \
                     sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])

    temp_crit = vsigma * ((1/tvort) + (mdot/be_mass)) / ((0.696 * np.pi * density * vel_r * G / cs**4) - ((2/3) / tcool))
    temp_crit = (temp_crit * (mu*mh/kb)).in_units('K').value
    temperature = temperature.in_units('K').value
    
    fpath = os.path.join(os.getcwd(), outdir, f"critical_temperature.h5")
    if (os.path.exists(fpath)):
        temp_info = h5py.File(fpath, 'a')
        temp_info['temperature'].resize((len(temp_info['temperature']) + 1,))
        temp_info['temperature'][-1] = temperature
        temp_info['temp_crit'].resize((len(temp_info['temp_crit']) + 1,))
        temp_info['temp_crit'][-1] = temp_crit
        temp_info['time'].resize((len(temp_info['time']) + 1,))
        temp_info['time'][-1] = time
    else :
        temp_info = h5py.File(fpath, 'w')
        temp_info.create_dataset('temperature', (1,), maxshape= (len(be_times),), dtype= 'f')
        temp_info['temperature'][0] = temperature
        temp_info.create_dataset('temp_crit', (1,), maxshape= (len(be_times),), dtype= 'f')
        temp_info['temp_crit'][0] = temp_crit
        temp_info.create_dataset('time', (1,), maxshape= (len(be_times),), dtype= 'f')
        temp_info['time'][0] = time




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
    
    ap.add_operation(critical_temperature, outdir='.')

    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
