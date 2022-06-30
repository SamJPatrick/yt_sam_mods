import statistics
import numpy as np
import os
import yt
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




def entropy_analysis(node, outdir='.', nbins=120):

    ds = node.ds
    sphere = node.sphere
    
    dx = ds.index.get_smallest_dx()    
    profile = my_profile(sphere, 'radius', \
                         [('gas', 'cooling_time'), ('gas', 'vortical_time'), \
                          ('gas', 'temperature'), \
                          ('gas', 'radial_velocity'), ('gas', 'velocity_spherical_phi')], \
                         logs={'radius': False}, n_bins= nbins, weight_field= ('gas', 'cell_mass'))
    mass_prof = my_profile(sphere, 'radius', ('gas', 'cell_mass'), logs={'radius': False}, n_bins= nbins, weight_field= None)
    radii = profile.x_bins[1:].in_units('pc')
    spheres = [ds.sphere(sphere.center, radius) for radius in radii]
    shells = [spheres[0]] + [spheres[i] - spheres[i-1] for i in range(1, len(spheres))]
    mass_tot_prof = [sum(sphere.quantities.total_mass()) for sphere in spheres]
    vsigma_prof = [np.std(shell[('gas', 'velocity_magnitude')])**2 for shell in shells]
    vsigma_prof = transpose_unyt(vsigma_prof).in_units('erg/g')
    #vsigma_tot_prof = [np.std(sphere[('gas', 'velocity_magnitude')])**2 for sphere in spheres]
    #vsigma_tot_prof = transpose_unyt(vsigma_tot_prof).in_units('erg/g')
    
    m_dot_prof = [0.0 for i in range (nbins)]
    #m_grd_prof = [0.0 for i in range (nbins)]
    e_dot_prof = [0.0 for i in range (nbins)]
    #vrad_grd_prof = [0.0 for i in range (nbins)]
    #vrad_dot_prof = [0.0 for i in range (nbins)]
    #vphi_dot_prof = [0.0 for i in range (nbins)]
    #temp_grd_prof = [0.0 for i in range (nbins)]
    #e_tot_dot_prof = [0.0 for i in range (nbins)]
    #m_tot_dot_prof = [0.0 for i in range (nbins)]
    
    for i in range (nbins):
        up_bin = (sphere['radius'] < radii[i]) & (sphere['radius'] >= radii[i] - dx)
        low_bin = (sphere['radius'] < radii[i-1]) & (sphere['radius'] >= radii[i-1] - dx)
        
        if (i == 0):
            up_m_flux = np.sum(sphere[('gas', 'mass_flux')][up_bin] * \
                               sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            up_e_flux = np.sum(sphere[('gas', 'energy_flux')][up_bin] * \
                               sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            #up_vrad = np.sum(sphere[('gas', 'radial_velocity')][up_bin] * \
            #              sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            #up_vphi = np.sum(sphere[('gas', 'velocity_spherical_phi')][up_bin] * \
            #              sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            #up_temp = np.sum(sphere[('gas', 'temperature')][up_bin] * \
            #              sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            #up_mass = np.sum(sphere[('gas', 'cell_mass')][up_bin]) / len(sphere['radius'][up_bin])
            
            m_dot_prof[i] = up_m_flux * radii[i]**2
            #m_grd_prof[i] = up_mass / radii[0]
            e_dot_prof[i] = up_e_flux * radii[i]**2
            #vrad_grd_prof[i] = up_vrad / radii[0]
            #vrad_dot_prof[i] = profile[('gas', 'radial_velocity')][0] * vrad_grd_prof[i]
            #vphi_dot_prof[i] = profile[('gas', 'radial_velocity')][0] * (up_vphi / radii[0])
            #temp_grd_prof[i] = up_temp / radii[0]
            #m_tot_dot_prof[i] = up_m_flux * radii[i]**2
            #e_tot_dot_prof[i] = up_e_flux * radii[i]**2
            
            continue
        
        if (len(sphere['radius'][up_bin]) == 0):
            up_m_flux = unyt_quantity(0.0, 'kg*s**(-1)*m**(-2)')
            up_e_flux = unyt_quantity(0.0, 'kg*s**(-3)')
            #up_vrad = unyt_quantity(0.0, 'm*s**(-1)')
            #up_vphi = unyt_quantity(0.0, 'm*s**(-1)')
            #up_temp = unyt_quantity(0.0, 'K')
            #up_mass = unyt_quantity(0.0, 'Msun')
        else :
            up_m_flux = np.sum(sphere[('gas', 'mass_flux')][up_bin] * \
                               sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            up_e_flux = np.sum(sphere[('gas', 'energy_flux')][up_bin] * \
                               sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            #up_vrad = np.sum(sphere[('gas', 'radial_velocity')][up_bin] * \
            #              sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            #up_vphi = np.sum(sphere[('gas', 'velocity_spherical_phi')][up_bin] * \
            #              sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            #up_temp = np.sum(sphere[('gas', 'temperature')][up_bin] * \
            #              sphere[('gas', 'cell_mass')][up_bin]) / np.sum(sphere[('gas', 'cell_mass')][up_bin])
            #up_mass = np.sum(sphere[('gas', 'cell_mass')][up_bin]) / len(sphere['radius'][up_bin])
            
        if (len(sphere['radius'][low_bin]) == 0):
            low_m_flux = unyt_quantity(0.0, 'kg*s**(-1)*m**(-2)')
            low_e_flux = unyt_quantity(0.0, 'kg*s**(-3)')
            #low_vrad = unyt_quantity(0.0, 'm*s**(-1)')
            #low_vphi = unyt_quantity(0.0, 'm*s**(-1)')
            #low_temp = unyt_quantity(0.0, 'K')
            #low_mass = unyt_quantity(0.0, 'Msun')
        else :
            low_m_flux = np.sum(sphere[('gas', 'mass_flux')][low_bin] * \
                               sphere[('gas', 'cell_mass')][low_bin]) / np.sum(sphere[('gas', 'cell_mass')][low_bin])
            low_e_flux = np.sum(sphere[('gas', 'energy_flux')][low_bin] * \
                                sphere[('gas', 'cell_mass')][low_bin]) / np.sum(sphere[('gas', 'cell_mass')][low_bin])
            #low_vrad = np.sum(sphere[('gas', 'radial_velocity')][low_bin] * \
            #               sphere[('gas', 'cell_mass')][low_bin]) / np.sum(sphere[('gas', 'cell_mass')][low_bin])
            #low_vphi = np.sum(sphere[('gas', 'velocity_spherical_phi')][low_bin] * \
            #              sphere[('gas', 'cell_mass')][low_bin]) / np.sum(sphere[('gas', 'cell_mass')][low_bin])
            #low_temp = np.sum(sphere[('gas', 'temperature')][low_bin] * \
            #              sphere[('gas', 'cell_mass')][low_bin]) / np.sum(sphere[('gas', 'cell_mass')][low_bin])
            #low_mass = np.sum(sphere[('gas', 'cell_mass')][low_bin]) / len(sphere['radius'][low_bin])

        m_dot_prof[i] = up_m_flux * radii[i]**2 - low_m_flux * radii[i-1]**2
        #m_grd_prof[i] = (up_mass - low_mass) / (radii[i] - radii[i-1])
        e_dot_prof[i] = up_e_flux * radii[i]**2 - low_e_flux * radii[i-1]**2
        #vrad_grd_prof[i] = (up_vrad - low_vrad) / (radii[i] - radii[i-1])
        #vrad_dot_prof[i] = profile[('gas', 'radial_velocity')][i] * vrad_grd_prof[i]
        #vphi_dot_prof[i] = profile[('gas', 'radial_velocity')][i] * (up_vphi - low_vphi) / (radii[i] - radii[i-1])
        #temp_grd_prof[i] = (up_temp - low_temp) / (radii[i] - radii[i-1])
        #m_tot_dot_prof[i] = up_m_flux * radii[i]**2
        #e_tot_dot_prof[i] = up_e_flux * radii[i]**2

    m_dot_prof = transpose_unyt(m_dot_prof).in_units('Msun/yr')
    #m_grd_prof = transpose_unyt(m_grd_prof).in_units('Msun/pc')
    e_dot_prof = transpose_unyt(e_dot_prof).in_units('erg/s')
    #vrad_grd_prof = transpose_unyt(vrad_grd_prof).in_units('s**(-1)')
    #vrad_dot_prof = transpose_unyt(vrad_dot_prof).in_units('m*s**(-2)')
    #vphi_dot_prof = transpose_unyt(vphi_dot_prof).in_units('m*s**(-2)')
    #temp_grd_prof = transpose_unyt(temp_grd_prof).in_units('K/m')
    #m_tot_dot_prof = transpose_unyt(m_tot_dot_prof).in_units('Msun/yr')
    #e_tot_dot_prof = transpose_unyt(e_tot_dot_prof).in_units('erg/s')

    mu = sphere.quantities.weighted_average_quantity(('gas', 'mean_molecular_weight'), ('gas', 'cell_mass'))
    temp_prof = ((kb/(mu*mh)) * profile[('gas', 'temperature')]).in_units('erg/g')
    #temp_grd_prof = (kb/(mu*mh)) * temp_grd_prof
    mass_prof = mass_prof[('gas', 'cell_mass')].in_units('Msun')
    #vrad_prof = profile[('gas', 'radial_velocity')].in_units('km/s')
    #vphi_prof = profile[('gas', 'velocity_spherical_phi')].in_units('km/s')
    tcool_prof = profile[('gas', 'cooling_time')].in_units('yr')
    tvort_prof = profile[('gas', 'vortical_time')].in_units('yr')
    mdot_turb_prof = (mass_prof / tvort_prof).in_units('Msun/yr')
    mdot_cool_prof = (mass_prof / tcool_prof).in_units('Msun/yr')

    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_shell_mrate_profiles.h5")
    m_dot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='m_dot')
    mdot_turb_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='turb_mrate')
    mdot_cool_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='cool_mrate')
    radii.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='radius')
    
    #egrv_dot_prof = mass_prof * (vrad_prof * \
    #                             ( (G / radii) * (m_grd_prof - mass_tot_prof / radii) \
    #                               + temp_grd_prof + (temp_prof / radii) ) \
    #                             + temp_prof * vrad_grd_prof)
    #ekin_dot_prof = m_dot_prof * ((1/2) * vrad_prof**2 + (1/3) * vphi_prof**2) + \
    #    mass_prof * (vrad_prof * vrad_dot_prof + (2/3) * vphi_prof * vphi_dot_prof)
    #eplus_dot_prof = ekin_dot_prof + e_dot_prof #+ egrv_dot_prof
    #alpha_prof = (temp_prof * (mdot_cool_prof + m_dot_prof) + (2/3)* (mdot_cool_prof * temp_prof + eplus_dot_prof)) / \
    #    (temp_prof * (mdot_cool_prof + m_dot_prof) - vsigma_prof * (mdot_turb_prof + m_dot_prof))
    #coolttheo_prof = - ((5/2) * mass_prof * temp_prof) / (eplus_dot_prof + (3/2) * m_dot_prof * temp_prof)
    #turbttheo_prof = mass_prof / ((2/5) * ((m_dot_prof * temp_prof - eplus_dot_prof) / vsigma_prof) - m_dot_prof)

    temp_term = temp_prof * (mdot_cool_prof + m_dot_prof)
    vsigma_term  = vsigma_prof * (mdot_turb_prof + m_dot_prof)
    alpha_prof = ((2.0/3.0) * e_dot_prof - temp_term) / (vsigma_term - temp_term)
    err1 = ((2.0/3.0) * e_dot_prof - vsigma_term) / ((2.0/3.0) * e_dot_prof + vsigma_term)
    err2 = (temp_term - (2.0/3.0) * e_dot_prof) / (temp_term + (2.0/3.0) * e_dot_prof)
    err3 = (temp_term - vsigma_term) / (temp_term + vsigma_term)
    coolttheo_prof = mass_prof / ((2/3) * (e_dot_prof / vsigma_prof) - m_dot_prof)
    turbttheo_prof = mass_prof / ((2/3) * (e_dot_prof / temp_prof) - m_dot_prof)
    
    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_shell_theory_profiles.h5")
    alpha_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='alpha')
    err1.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='error_1')
    err2.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='error_2')
    err3.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='error_3')
    turbttheo_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='turb_time_theory')
    coolttheo_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='cool_time_theory')
    radii.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='radius')

    '''
    angmom_tot_prof = [sp.quantities.total_quantity(('gas', 'specific_angular_momentum_magnitude')) for sp in spheres]
    angmom_tot_prof = transpose_unyt(angmom_tot_prof).in_units('m**2/s')
    temp_tot_prof = (kb/(mu*mh)) * [sp.quantities.weighted_average_quantity(('gas', 'temperature'), ('gas', 'cell_mass')) for sp in spheres]
    temp_tot_prof = transpose_unyt(temp_tot_prof).in_units('erg/g')
    tcool_tot_prof = [sp.quantities.total_quantity(('gas', 'cooling_time')) for sp in spheres]
    tcool_tot_prof = transpose_unyt(tcool_tot_prof).in_units('yr')
    tvort_tot_prof = [sp.quantities.total_quantity(('gas', 'vortical_time')) for sp in spheres]
    tvort_tot_prof = transpose_unyt(tcool_tot_prof).in_units('yr')
    mdot_turb_tot_prof = (mass_tot_prof / tvort_tot_prof).in_units('Msun/yr')
    mdot_cool_tot_prof = (mass_tot_prof / tcool_tot_prof).in_units('Msun/yr')

    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_total_mrate_profiles.h5")
    m_tot_dot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='m_dot')
    mdot_turb_tot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='turb_mrate')
    mdot_cool_tot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='cool_mrate')
    radii.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='radius')
        
    ds.add_field(('gas', 'specific_kinetic_energy'), function=_specific_kinetic_energy, sampling_type='local', units='erg/g', take_log=False)
    profile = my_profile(sphere, 'radius', [('gas', 'specific_kinetic_energy'), ('gas', 'specific_thermal_energy')], \
                         n_bins= nbins, logs={'radius': False})
    
    pe_prof = mass_prof * G / radii
    ke_prof = profile[('gas', 'specific_kinetic_energy')].in_units('erg/g')
    te_prof = profile[('gas', 'specific_thermal_energy')].in_units('erg/g')
    turb_prof = (mu*mh/kb) * vsigma_prof

    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_shell_energy_profiles.h5")
    pe_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='potential')
    ke_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='kinetic')
    te_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='thermal')
    turb_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='turbulent')
    radii.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='radius')

    pe_tot_prof = (5 * mass_tot_prof**2) / (3 * radii)
    ke_tot_prof = [sp.quantities.total_quantity(('gas', 'specific_kinetic_energy')) * mass_tot_prof for sp in spheres]
    ke_tot_prof = transpose_unyt(ke_tot_prof).in_units('erg')
    te_tot_prof = [sp.quantities.total_quantity(('gas', 'specific_thermal_energy')) * mass_tot_prof for sp in spheres]
    te_tot_prof = transpose_unyt(te_tot_prof).in_units('erg')
    
    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_total_energy_profiles.h5")
    pe_tot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='potential')
    ke_tot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='kinetic')
    te_tot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='thermal')
    vsigma_tot_prof.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='turbulent')
    radii.write_hdf5(fpath, info={'time':ds.current_time, 'dump': str(ds)}, group_name='radius')
    '''


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
    
    ap.add_operation(entropy_analysis, outdir= 'Profiles/Entropy_simple', nbins= 120)
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
