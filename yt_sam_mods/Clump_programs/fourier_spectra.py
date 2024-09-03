import os
import yt
import ytree
yt.enable_parallelism()

import numpy as np
import scipy.fft as fft

from ytree.analysis import AnalysisPipeline

from yt.extensions.sam_mods.misc import *
from yt.extensions.sam_mods.profiles import my_profile
from yt.extensions.sam_mods.plots import *
from yt.extensions.sam_mods.tree_analysis_operations import yt_dataset, \
    node_profile, garbage_collect, delattrs
from yt.extensions.sam_mods import add_p2p_fields

from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh, pc

#import pdb



VEL_UNIT = 'km/s'
BIN_DENSITY = 20
CUBE_SIZE = unyt_quantity(50.0, 'pc') #old: 500 pc
LEVEL_N = 11 # old: 10
OUTDIR = "."


def coord_generator(cell_num):

    for i in range (0, cell_num + 1):
        for j in range (0, cell_num - i):
            k = cell_num - i - j - 1
            yield [i, j, k]


def fourier_spectra(node, outdir='.'):

    ds = node.ds
    sp = node.sphere
    cube_size = ds.quan(CUBE_SIZE.to('pc').value * 2, 'pc')
    n_cells = 2**LEVEL_N
    cell_size = cube_size/n_cells
    offset_vec = ds.arr([cube_size/2, cube_size/2, cube_size/2], 'pc')
    cube = ds.arbitrary_grid(sp.center - offset_vec, sp.center + offset_vec, 
                             dims=[n_cells, n_cells, n_cells])

    inv_cube = np.abs(fft.fftn(cube[('gas', 'velocity_magnitude')].to(VEL_UNIT)))
    yfs_mean = unyt_array(np.zeros(n_cells), VEL_UNIT)
    yfs_diag = unyt_array(np.zeros(n_cells), VEL_UNIT)
    for cell_num in range (1, n_cells):
        #coord_list = []
        #for i in range (0, cell_num + 1):
        #    for j in range (0, cell_num - i):
        #        k = cell_num - i - j - 1
        #        coord_list.append([i, j, k])
        yfs_mean[cell_num] = np.mean([inv_cube[i][j][k] for i, j, k in coord_generator(cell_num)])
        yfs_diag[cell_num] = inv_cube[cell_num][cell_num][cell_num]
    freqs = unyt_array(fft.rfftfreq(n_cells, cell_size.value), 1/cube_size.units)

    fpath = os.path.join(outdir, f"{str(ds)}_velocity_mag_fourier_spectra.h5")
    yfs_mean.write_hdf5(fpath, group_name= f'velocity_fourier_globmean')
    yfs_diag.write_hdf5(fpath, group_name= f'velocity_fourier_diagonals')
    freqs.write_hdf5(fpath, group_name= 'x_data')
    




if __name__ == "__main__":

    #output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/minihalo_analysis"
    #sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/simulation.h5"
    #data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust"
    #tree_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"

    output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust/minihalo_analysis"
    sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust/simulation.h5"
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust"
    tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/hyper_512_collapse_solar_dust/merger_trees/target_halos/target_halos.h5"

    #output_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/minihalo_analysis"
    #sim_path = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/simulation.h5"
    #data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    #tree_path= "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo/merger_trees/target_halos/target_halos.h5"

    #es = yt.load(sim_path)
    a = ytree.load(tree_path)
    if "icom_gas2_position_x" in a.field_list:
        a.add_vector_field("icom_gas2_position")

    ap = AnalysisPipeline()
    ap.add_operation(yt_dataset, data_dir, add_fields=True)
    #ap.add_operation(yt_dataset, data_dir, es)
    #ap.add_operation(modify_grackle)
    ap.add_operation(return_sphere)
    ap.add_operation(align_sphere)
    
    ap.add_operation(fourier_spectra, outdir= OUTDIR)
    
    ap.add_operation(delattrs, ["sphere", "ds"], always_do=True)
    ap.add_operation(garbage_collect, 60, always_do=True)

    tree = a[0]
    # N_MISSED = 19
    # tree_mod = list(tree['prog'])[N_MISSED:]
    # for node in ytree.parallel_trees(tree_mod):
    #     ap.process_target(node)
    for node in ytree.parallel_tree_nodes(tree, group="prog"):
        ap.process_target(node)
