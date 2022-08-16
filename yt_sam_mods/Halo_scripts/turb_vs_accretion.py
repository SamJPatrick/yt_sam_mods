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




def turb_vs_accretion(node, outdir='.', nbins=120, take_log= False):

    ds = node.ds
    sphere = node.sphere
    
    dx = ds.index.get_smallest_dx()        
    profile = my_profile(sphere, 'radius', [('gas', 'vturb')], logs={'radius': take_log}, n_bins= nbins, weight_field= ('gas', 'cell_mass'))
    radii = profile.x_bins[1:].in_units('pc')
    vturb_shell_prof = profile[('gas', 'vturb')].in_units('km/s')
    
    spheres = [ds.sphere(sphere.center, radius) for radius in radii]
    shells = [spheres[0]] + [spheres[i] - spheres[i-1] for i in range(1, len(spheres))]
    vsigma_shell_prof = [np.std(shell[('gas', 'velocity_magnitude')]) for shell in shells]
    vsigma_shell_prof = transpose_unyt(vsigma_shell_prof).in_units('km/s')
    vsigma_tot_prof = [np.std(sphere[('gas', 'velocity_magnitude')]) for sphere in spheres]
    vsigma_tot_prof = transpose_unyt(vsigma_tot_prof).in_units('km/s')
    vturb_tot_prof = [sphere.quantities.weighted_average_quantity(('gas', 'vturb'), ('gas', 'cell_mass')) for sphere in spheres]
    vturb_tot_prof = transpose_unyt(vturb_tot_prof).in_units('km/s')
    
    mdot_shell_prof = [0.0 for i in range (nbins)]
    mdot_tot_prof = [0.0 for i in range (nbins)]
    edot_shell_prof = [0.0 for i in range (nbins)]
    edot_tot_prof = [0.0 for i in range (nbins)]
 
    for i in range (nbins):
        
        up_bin = (sphere['radius'] < radii[i]) & (sphere['radius'] >= radii[i] - dx)
        low_bin = (sphere['radius'] < radii[i-1]) & (sphere['radius'] >= radii[i-1] - dx)
        
        if (i == 0):
            up_m_flux = np.sum(sphere[('gas', 'mass_flux')][up_bin] * \
                               sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            up_e_flux = np.sum(sphere[('gas', 'energy_flux')][up_bin] * \
                               sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            m_dot_prof[i] = up_m_flux * radii[i]**2
            e_dot_prof[i] = up_e_flux * radii[i]**2            
            continue
        
        if (len(sphere['radius'][up_bin]) == 0):
            up_m_flux = unyt_quantity(0.0, 'kg*s**(-1)*m**(-2)')
            up_e_flux = unyt_quantity(0.0, 'kg*s**(-3)')
        else :
            up_m_flux = np.sum(sphere[('gas', 'mass_flux')][up_bin] * \
                               sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            up_e_flux = np.sum(sphere[('gas', 'energy_flux')][up_bin] * \
                               sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])

        if (len(sphere['radius'][low_bin]) == 0):
            low_m_flux = unyt_quantity(0.0, 'kg*s**(-1)*m**(-2)')
            low_e_flux = unyt_quantity(0.0, 'kg*s**(-3)')
        else :
            low_m_flux = np.sum(sphere[('gas', 'mass_flux')][low_bin] * \
                                sphere[('gas', 'cell_mass')][low_bin]) / np.sum(sphere[('gas', 'cell_mass')][low_bin])
            low_e_flux = np.sum(sphere[('gas', 'energy_flux')][low_bin] * \
                                sphere[('gas', 'cell_mass')][low_bin]) / np.sum(sphere[('gas', 'cell_mass')][low_bin])

        mdot_shell_prof = up_m_flux * radii[i]**2 - low_m_flux * radii[i-1]**2
        mdot_tot_prof[i] = up_m_flux * radii[i]**2
        edot_shell_prof[i] = up_e_flux * radii[i]**2 - low_e_flux * radii[i-1]**2
        edot_tot_prof[i] = up_e_flux * radii[i]**2

    mdot_shell_prof = transpose_unyt(mdot_shell_prof).in_units('Msun/yr')
    mdot_tot_prof = transpose_unyt(mdot_tot_prof).in_units('Msun/yr')
    edot_shell_prof = transpose_unyt(edot_shell_prof).in_units('erg/s')
    edot_tot_prof = transpose_unyt(edot_tot_prof).in_units('erg/s')
    
    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_turb_vs_accretion_shell.h5")
    vsigma_shell_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='vsigma')
    vturb_shell_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='vturb')
    mdot_shell_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='mass_dot')
    edot_shell_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='energy_dot')
    radii.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='radius')

    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_turb_vs_accretion_tot.h5")
    vsigma_tot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='vsigma')
    vturb_tot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='vturb')
    mdot_tot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='mass_dot')
    edot_tot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='energy_dot')
    radii.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='radius')



def perform_test(node):

    ds = node.ds
    sp = ds.sphere(unyt_array([2000, 2000, 2000], 'pc'), unyt_quantity(200, 'pc'))
    y = sp.quantities.weighted_average_quantity(('gas', 'mass_flux'), ('gas', 'cell_mass'))
    print(y)
    x = sp.quantities.weighted_average_quantity(('gas', 'energy_flux'), ('gas', 'cell_mass'))
    print(x)


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
    ap.add_operation(yt_dataset, data_dir)
    #ap.add_operation(yt_dataset, data_dir, es)
    ap.add_operation(perform_test)
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)
    
    ap.add_operation(turb_vs_accretion, outdir= 'Profiles/Turb_vs_accretion', nbins= 120)

    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
