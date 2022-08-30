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




def plot_velocity_profiles(node, indir= '.', outdir='.'):

    ds = node.ds    
    my_fig = GridFigure(1, 1, figsize=(6, 4.5),
                        top_buffer = 0.14, bottom_buffer = 0.12,
                        left_buffer = 0.09, right_buffer = 0.02)
    my_axes = my_fig[0]
    my_axes.set_xscale('log')

    fields = ["velocity_magnitude",
              "tangential_velocity_magnitude",
              "velocity_spherical_radius",
              "sound_speed"]
    colors = ["red", "green", "blue", "purple", "#544D4D", "black"]
    # http://www.colourlovers.com/palette/1329926/graphic_artist.
    # colors = ["#B00029", "#90B004", "#19849C", "#851370", "#544D4D", "black"]
    for i, field in enumerate(fields):
        filename = os.path.join(indir, f"{str(ds)}_2D_profile_radius_{field}_None.h5")
        plot_profile_distribution(my_axes, filename, 'cell_mass',
                                  x_units="pc", y_units='km/s', alpha_scale=0.7,
                                  pkwargs=dict(color=colors[i], linewidth=1))

    fn = os.path.join(indir, f"{str(ds)}_1D_profile_radius_cell_mass.h5")
    pds = yt.load(fn)
    pradius = pds.profile.x.to("pc")
    vsigma = pds.profile.standard_deviation['data', 'velocity_magnitude'].to("km/s")
    my_axes.plot(pradius[vsigma > 0], vsigma[vsigma > 0], alpha=0.9,
                 linewidth=1, color=colors[4], zorder=998)

    field = "matter_mass"
    fn = os.path.join(indir, f"{str(ds)}_1D_profile_radius_None.h5")
    mds = yt.load(fn)
    radius = mds.profile.x.to("pc")
    mass = mds.profile[field]
    dfil = mass > 0
    v_sp = np.sqrt(G * mass[dfil].cumsum() / radius[dfil]).to("km/s")
    my_axes.plot(radius[dfil], v_sp, alpha=0.9, linewidth=1,
                 color=colors[5], zorder=997)

    ylim = (-5, 13)
    ymajor = np.arange(-5, 16, 5.)
    yminor = np.arange(-5, 15, 1.)
    my_axes.yaxis.set_label_text("v [km / s]")
    my_axes.yaxis.labelpad = -3
    draw_major_grid(my_axes, 'y', ymajor,color='black', linestyle='-',
                    linewidth=1, alpha=0.2)
    ty = mirror_yticks(my_axes, ylim, ymajor, yminor=yminor)
    xlim = (1e-1, 2e2)
    tx = twin_unit_axes(my_axes, xlim, "r", "pc", top_units="AU")
    labels = ["|v|", "v$_{\\rm tan}$", "v$_{\\rm rad}$",
              "c$_{\\rm s}$", "$\\sigma$", "v$_{\\rm c}$"]
    dist = [True]*4 + [False]*2
    plot_items = list(zip(colors, labels, dist))
    plot_profile_distribution_legend(my_axes, plot_items, alpha_scale=0.7)
    fpath = os.path.join(os.getcwd(), outdir, f"{str(ds)}_velocity_profiles.png")
    pyplot.savefig(fpath)
    pyplot.close()



def plot_timescale_profiles(node, indir='.', outdir='.'):

    my_fig = GridFigure(1, 1, figsize=(6, 4.5),
                        top_buffer = 0.13, bottom_buffer = 0.12,
                        left_buffer = 0.14, right_buffer = 0.02,
                        horizontal_buffer = 0.05, vertical_buffer = 0)
    my_axes = my_fig[0]
    my_axes.set_xscale('log')
    my_axes.set_yscale('log')

    xlim = (5e-3, 3e2)
    tx = twin_unit_axes(my_axes, xlim, "r", "pc", top_units="AU")
    fields = ["turbulent_sound_crossing_time",
              "sound_crossing_time",
              "total_dynamical_time",
              "cooling_time",
              "vortical_time"]
    units = ["yr",
             "yr",
             #f"{np.sqrt(2)}*yr",
             "yr",
             "yr",
             #f"{1/(2*np.pi)}*yr",
             "yr"]
    labels = ["turbulent sound-crossing",
              "sound-crossing",
              "free-fall",
              "cooling",
              "mixing"]
    colors = ["orange", "red", "green", "blue"]

    dsfn = node.ds_filename.split('/')[-1]
    filename = os.path.join(indir, f"{str(dsfn)}_1D_profile_radius_cell_mass.h5")
    df = yt.load(filename)
    used = df.profile.used
    x_data = df.profile.x[used].to("pc")
    
    profile_dict = {}
    for field, unit, label, color in zip(fields, units, labels, colors):
        if (field == "sound_crossing_time"):
            cs = df.profile[('data', 'sound_speed')][used]
            y_sc = (2 * x_data / cs)
            profile_dict[field] = y_sc.to(unit)
        elif (field == "turbulent_sound_crossing_time"):
            #cs = df.profile[('data', 'sound_speed')][used]
            vt = df.profile.standard_deviation[('data', 'velocity_magnitude')][used]
            #v = np.sqrt(cs**2 + vt**2)
            y_sct = (2 * x_data / vt)
            profile_dict[field] = y_sct.to(unit)
        else:
            profile_dict[field] = df.profile[('data', field)][used].to(unit)

        my_axes.plot(x_data, profile_dict[field], color=color,
                     alpha=0.7, linewidth=1.5, label=label)

    ylim = (1e4, 1e9)
    ymajor = np.logspace(2, 10, 5)
    yminor = np.logspace(1, 9, 5)
    draw_major_grid(my_axes, 'y', ymajor,
                color='black', linestyle='-',
                linewidth=1, alpha=0.2)
    ty = mirror_yticks(my_axes, ylim, ymajor, yminor=yminor)

    bkgcolor = ['blue', 'green', 'yellow']
    for i in range (len(x_data) - 1):
        if (profile_dict['cooling_time'][i] < profile_dict['total_dynamical_time'][i]):
            pyplot.axvspan(x_data[i], x_data[i+1], facecolor= bkgcolor[0], alpha=0.5)
        elif (profile_dict['total_dynamical_time'][i] < profile_dict['turbulent_sound_crossing_time'][i] and \
              profile_dict['total_dynamical_time'][i] < profile_dict['sound_crossing_time'][i]):
            pyplot.axvspan(x_data[i], x_data[i+1], facecolor= bkgcolor[1], alpha=0.5)
        elif (profile_dict['total_dynamical_time'][i] < profile_dict['sound_crossing_time'][i]):
            pyplot.axvspan(x_data[i], x_data[i+1], facecolor= bkgcolor[2], alpha=0.5)
    
    my_axes.yaxis.set_label_text("t [yr]")
    my_axes.legend()
    fpath = os.path.join(os.getcwd(), outdir, f"{str(dsfn)}_timescale_profiles.png")
    pyplot.savefig(fpath)
    pyplot.close()



def create_state_data(node, indir='.', outdir='.'):
   
    dsfn = node.ds_filename.split('/')[-1]
    filename = os.path.join(indir, f"{str(dsfn)}_1D_profile_radius_cell_mass.h5")
    df = yt.load(filename)
    used = df.profile.used
    x_data = df.profile.x[used].to("pc")

    profile_dict = {}
    fields = ["turbulent_sound_crossing_time", "sound_crossing_time", "total_dynamical_time", \
              "cooling_time", "vortical_time"]
    for field in fields:
        if (field == "sound_crossing_time"):
            cs = df.profile[('data', 'sound_speed')][used]
            y_sc = (2 * x_data / cs)
            profile_dict[field] = y_sc.to('yr')
        elif (field == "turbulent_sound_crossing_time"):
            #cs = df.profile[('data', 'sound_speed')][used]
            vt = df.profile.standard_deviation[('data', 'velocity_magnitude')][used]
            #v = np.sqrt(cs**2 + vt**2)
            y_sct = (2 * x_data / vt)
            profile_dict[field] = y_sct.to('yr')
        else:
            profile_dict[field] = df.profile[('data', field)][used].to('yr')

    cell_mass = df.profile.weight[used].in_units('Msun')
    state_dict = {'frag':unyt_array([0.0 for i in range (len(x_data))], 'Msun'),
                  'support':unyt_array([0.0 for i in range (len(x_data))], 'Msun'),
                  'collapse':unyt_array([0.0 for i in range (len(x_data))], 'Msun')}
    for i in range (len(x_data) - 1):
        if (profile_dict['cooling_time'][i] < profile_dict['total_dynamical_time'][i]):
            state_dict['frag'][i] = cell_mass[i]
        elif (profile_dict['total_dynamical_time'][i] < profile_dict['turbulent_sound_crossing_time'][i] and \
              profile_dict['total_dynamical_time'][i] < profile_dict['sound_crossing_time'][i]):
            state_dict['collapse'][i] = cell_mass[i]
        elif (profile_dict['total_dynamical_time'][i] < profile_dict['sound_crossing_time'][i]):
            state_dict['support'][i] = cell_mass[i]

    gas_mass = unyt_array(np.zeros(len(state_dict.keys())), 'Msun')
    gas_frac = unyt_array(np.zeros(len(state_dict.keys())), '')
    mean_radius = unyt_array(np.zeros(len(state_dict.keys())), 'pc')
    radius_frac = unyt_array(np.zeros(len(state_dict.keys())), '')
    for i, state in enumerate(state_dict.keys()):
        gas_mass[i] = sum(state_dict[state])
        gas_frac[i] = gas_mass[i] / sum(cell_mass)
        mean_radius[i] = sum([x_data[i] * mass for i, mass in enumerate(state_dict[state])]) / gas_mass[i]
        radius_frac[i] = mean_radius[i] / x_data[-1]
    
    state_fn = os.path.join(outdir, f"{str(dsfn)}_state_info.h5")
    gas_mass.write_hdf5(state_fn, dataset_name='gas mass')
    gas_frac.write_hdf5(state_fn, dataset_name='gas fraction')
    mean_radius.write_hdf5(state_fn, dataset_name='mean radius')
    radius_frac.write_hdf5(state_fn, dataset_name='radius fraction')







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
    ap.add_operation(yt_dataset, data_dir, add_fields= False, load_data=False)
   
    ap.add_operation(plot_velocity_profiles, indir='Profiles/Velocities_and_timescales', outdir='Profiles/Velocities')
    ap.add_operation(plot_timescale_profiles, indir='Profiles/Velocities_and_timescales', outdir='Profiles/Timescales')
    ap.add_operation(create_state_data, indir='Profiles/Velocities_and_timescales', outdir='State_data')
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
