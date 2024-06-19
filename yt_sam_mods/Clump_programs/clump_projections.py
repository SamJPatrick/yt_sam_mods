import os
import yt
yt.enable_parallelism()
import numpy as np

#from yt.extensions.sam_mods.misc import *
#from yt.extensions.sam_mods.profiles import my_profile
#from yt.extensions.sam_mods.tree_analysis_operations import *
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
#from yt.extensions.sam_mods.graph_funcs import get_field_dict
#from yt.extensions.sam_mods import add_p2p_fields
from unyt import unyt_quantity, unyt_array


#DATA_DIR = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/pisn_solo"
DATA_DIR = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_collapse_solar_dust/"
#DUMP_NUM = "DD0295/DD0295"
DUMP_NUM = "DD0179/DD0179"
#CLUMP_DIR = "Clumps/DD0295_clump_info_pisn.h5"
CLUMP_DIR = "Clumps/DD0179_clump_info_ccsn.h5"
OUTDIR = "Clump_projections_ccsn_179"

PLOT_FIELDS = [('gas', 'density')]
DENSITY_LIMS = (1e-25, 1e-15)
VOL_THRESH = unyt_quantity(5e-12, 'pc**3')
ENLARGEMENT = 2

ds = yt.load(os.path.join(DATA_DIR, DUMP_NUM))
df = yt.load(CLUMP_DIR)
clumps = [clump for clump in df.tree if clump[('clump', 'volume')] > VOL_THRESH]
for clump in clumps:
    position = ds.arr(unyt_array(clump[('clump', 'com')], 'pc'))
    width = clump[('clump', 'volume')]**(1/3) * ENLARGEMENT
    region = ds.box(position - width, position + width)
    n = clump[('clump', 'clump_id')]
    for ax in 'xyz':
        p = yt.ProjectionPlot(ds, ax, PLOT_FIELDS, weight_field=('gas', 'density'),
                              center= position, width= width, data_source=region)
        p.set_axes_unit('pc')
        p.set_cmap(('gas', 'density'), "turbo")
        p.set_zlim(('gas', 'density'), *DENSITY_LIMS)
        children = clump.children
        if (children != None):
            children.append(clump)
            p.annotate_clumps(children)
        p.save(os.path.join(OUTDIR, f'{str(ds)}_clump_{n}_{ax}.png'))
