import numpy as np
import os
import sys
import yt
import ytree

from yt.data_objects.level_sets.api import Clump, find_clumps
from yt.extensions.sam_mods.fields import add_p2p_fields
from yt.extensions.sam_mods.clump_info_items import *
from yt.extensions.sam_mods.clump_validators import *
from yt.extensions.sam_mods.misc import iterate_center_of_mass
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt

from unyt import unyt_quantity, unyt_array
from unyt import G, Msun, kb, mh



RECENTER_RADIUS = unyt_quantity(40.0, 'pc')
SEARCH_RADIUS = unyt_quantity(20.0, 'pc')
INNER_RADIUS = unyt_quantity(1000.0, 'AU')
CLUMP_STEP = 5.0

OUTPUT_DIR = "Clumps"
PROJ_SUBDIR = "Projections"


CLUMPING_FIELD = ('gas', 'density')
INFO_FIELDS = ["min_number_density", "mean_number_density", "max_number_density", "min_metallicity", "mean_metallicity", "max_metallicity", \
               "min_temperature", "mean_temperature", "max_temperature", "jeans_mass", "volume", "mass", "metal_mass", "com", "a_max", \
               "ratio", "alpha"]
PLOT_FIELDS = [('gas', 'number_density'), ('gas', 'temperature')]


if __name__ == "__main__":
    
    data_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
    tree_path = "merger_trees/target_halos/target_halos.h5"
    sim_path = "DD0295/DD0295"
    
    a = ytree.load(os.path.join(data_dir, tree_path))
    ds = yt.load(os.path.join(data_dir, sim_path))
    add_p2p_fields(ds)

    center = ds.arr(a[0]["icom_gas2_position"], 'code_length')
    sp0 = ds.sphere(center, RECENTER_RADIUS)
    new_center = ds.arr(transpose_unyt(sp0.quantities.max_location(("gas", "density"))[1:]))
    sp = ds.sphere(new_center, SEARCH_RADIUS)
    del sp0
    
    c_min = 10**np.floor(np.log10(sp[CLUMPING_FIELD]).min())
    c_max = 10**np.floor(np.log10(sp[CLUMPING_FIELD]).max() + 1)
    yt.mylog.info(f"Density extrema: {c_min} - {c_max}.")
    master_clump = Clump(sp, CLUMPING_FIELD)
    master_clump.add_validator("future_bound", use_thermal_energy=True, truncate=True, include_cooling=True)
    ds.add_field(('gas', 'metal_laplacian'), function=generate_vol_laplacian(ds, ('gas', 'metallicity3')), \
                 sampling_type='local', units='pc', take_log=False)
    for field in INFO_FIELDS:
        master_clump.add_info_item(field)
    find_clumps(master_clump, c_min, c_max, CLUMP_STEP)
    fn = master_clump.save_as_dataset(filename=os.path.join(OUTPUT_DIR, f'{str(ds)_clump_info.h5}'), fields=["density"])
    leaf_clumps = master_clump.leaves

    '''
    pdir = os.path.join(OUTPUT_DIR, PROJ_SUBDIR)
    units = 'pc'
    dx = ds.index.get_smallest_dx()
    for i, sphere in enumerate(
            iterate_center_of_mass(sp, INNER_RADIUS, com_kwargs={'use_gas': True, 'use_particles': False})):
        if (sphere.radius < 10 * dx):
            break
        width = 2 * sphere.radius
        region = ds.box(sphere.center-(1.05 * width / 2), sphere.center+(1.05 * width / 2))
        if width.to('pc') < 0.01:
            units = 'AU'
        for ax in 'xyz':
            p = yt.ProjectionPlot(ds, ax, PLOT_FIELDS, weight_field=('gas', 'density'), \
                                  center=sphere.center, width=width, data_source=region)
            p.set_axes_unit(units)
            p.set_cmap(('gas', 'number_density'), "turbo")
            p.annotate_clumps(leaf_clumps)
            p.save(os.path.join(pdir, "%03d" % i))
    '''
