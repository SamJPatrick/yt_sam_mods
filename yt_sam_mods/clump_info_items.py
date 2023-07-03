import numpy as np
from yt.data_objects.level_sets.api import add_clump_info
from unyt import G, kb, mh, pc
from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.field_utils import *




Z_CRIT = 5e-3
SMOL_SP_INIT_FACT = 1.10 * np.sqrt(3.0)
SMOL_SP_INC_FACT = 1.5




def _total_volume(clump):
    volume = clump.data.quantities.total_quantity(('gas', 'cell_volume')).in_units('pc**3')
    return ("Total clump volume: %.3e pc^3", volume)
add_clump_info("volume", _total_volume)

def _total_mass(clump):
    mass = clump.data.quantities.total_mass()[0].in_units('Msun')
    return ("Total clump mass: %.3e Msun", mass)
add_clump_info("mass", _total_mass)

def _metal_mass(clump):
    metal_mass = clump.data.quantities.total_quantity(('gas', 'metal3_mass')).in_units('Msun')
    return ("Clump mass in metals: %.3e Msun", metal_mass)
add_clump_info("metal_mass", _metal_mass)

def _center_of_mass(clump):
    com = clump.data.quantities.center_of_mass(True, False).in_units('pc')
    return ("Clump center of mass: %.5f, %.5f, %.5f pc", (com[0], com[1], com[2]))
add_clump_info("com", _center_of_mass) 




def _min_number_density(clump):
    min_n = clump.data["gas", "number_density"].min().in_units("cm**-3")
    return "Min number density: %.6e cm^-3.", min_n
add_clump_info("min_number_density", _min_number_density)

def _mean_number_density(clump):
    mean_n = sum(clump.data["gas", "number_density"] * clump.data["gas", "cell_volume"]) / \
        sum(clump.data["gas", "cell_volume"])
    return "Mean number density: %.6e cm^-3.", mean_n
add_clump_info("mean_number_density", _mean_number_density)

def _max_number_density(clump):
    max_n = clump.data["gas", "number_density"].max().in_units("cm**-3")
    return "Max number density: %.6e cm^-3.", max_n
add_clump_info("max_number_density", _max_number_density)

def _min_metallicity(clump):
    min_Z = clump.data["gas", "metallicity3"].min()
    return "Min metallicity: %.3e Zsolar", min_Z
add_clump_info("min_metallicity", _min_metallicity)

def _mean_metallicity(clump):
    mean_Z = sum(clump.data["gas", "metallicity3"] * clump.data["gas", "cell_mass"]) / \
        sum(clump.data["gas", "cell_mass"])
    return "Mean metallicity: %.3e Zsolar", mean_Z
add_clump_info("mean_metallicity", _mean_metallicity)

def _max_metallicity(clump):
    max_Z = clump.data["gas", "metallicity3"].max()
    return "Max metallicity: %.3e Zsolar", max_Z
add_clump_info("max_metallicity", _max_metallicity)

def _min_temperature(clump):
    min_temp = clump.data["gas", "temperature"].min()
    return "Min temperature: %.1f K", min_temp
add_clump_info("min_temperature", _min_temperature)

def _mean_temperature(clump):
    mean_temp = sum(clump.data["gas", "temperature"] * clump.data["gas", "cell_mass"]) / \
        sum(clump.data["gas", "cell_mass"])
    return "Mean temperature: %.1f K", mean_temp
add_clump_info("mean_temperature", _mean_temperature)

def _max_temperature(clump):
    max_temp = clump.data["gas", "temperature"].max()
    return "Max temperature: %.1f K", max_temp
add_clump_info("max_temperature", _max_temperature)

def _jeans_mass(clump):
    temperature = clump.data.quantities.weighted_average_quantity(
        ("gas", "temperature"), ("gas", "cell_mass"))
    density = clump.data.quantities.weighted_average_quantity(
        ("gas", "density"), ("index", "cell_volume"))
    mu = clump.data.quantities.weighted_average_quantity(
        ("gas", "mean_molecular_weight"), ("gas", "cell_mass"))
    MJ_constant = (((5.0 * kb) / (G * mh))**(1.5)) * (3.0 / (4.0 * np.pi))**(0.5)
    u = MJ_constant * ((temperature / mu)**(1.5)) * (density**(-0.5))
    return "Jeans mass: %.6e Msun.", u.in_units("Msun")
add_clump_info("jeans_mass", _jeans_mass)




def _amax(clump):

    com = clump.info['com'][1]
    clump_mass = clump.info['mass'][1]

    x_x = np.sum(clump.data[('gas', 'cell_mass')] * \
                 ( ((clump.data[('gas', 'y')] - com[1]).in_units('pc'))**2 + \
                   ((clump.data[('gas', 'z')] - com[2]).in_units('pc'))**2 ) ) / clump_mass
    y_y = np.sum(clump.data[('gas', 'cell_mass')] * \
                 ( ((clump.data[('gas', 'x')] - com[0]).in_units('pc'))**2 + \
                   ((clump.data[('gas', 'z')] - com[2]).in_units('pc'))**2 ) ) / clump_mass
    z_z = np.sum(clump.data[('gas', 'cell_mass')] * \
                 ( ((clump.data[('gas', 'x')] - com[0]).in_units('pc'))**2 + \
                   ((clump.data[('gas', 'y')] - com[1]).in_units('pc'))**2 ) ) / clump_mass
    x_y = - np.sum(clump.data[('gas', 'cell_mass')] * \
                   (clump.data[('gas', 'x')] - com[0]).in_units('pc') * \
                   (clump.data[('gas', 'y')] - com[1]).in_units('pc')) / clump_mass        
    x_z = - np.sum(clump.data[('gas', 'cell_mass')] * \
                   (clump.data[('gas', 'x')] - com[0]).in_units('pc') * \
                   (clump.data[('gas', 'z')] - com[2]).in_units('pc')) / clump_mass
    y_z = - np.sum(clump.data[('gas', 'cell_mass')] * \
                   (clump.data[('gas', 'y')] - com[1]).in_units('pc') * \
                   (clump.data[('gas', 'z')] - com[2]).in_units('pc')) / clump_mass
    I = np.array([[x_x, x_y, x_z], [x_y, y_y, y_z], [x_z, y_z, z_z]])
    
    eigvalues, eigvectors = np.linalg.eig(I)
    a_max = max(clump.data[('gas', 'dx')].min().to('pc'), unyt_quantity(np.sqrt(max(eigvalues)), 'pc'))
    return ("Clump characteristic length scale: %.5f pc", a_max)
add_clump_info("a_max", _amax)

def _ratio(clump):

    ds = clump.data.ds
    com = clump.info['com'][1]
    a_max = clump.info['a_max'][1]
    clump_volume = clump.info['volume'][1]
    
    factor = SMOL_SP_INIT_FACT
    min_area = (24 * clump.data[('gas', 'dx')].min()**2).in_units('pc**2')
    clump_area = min_area
    i = 0
    while (clump_area <= min_area):
        if (i > 10):
            break
        i += 1
        try :
            smol_sp = ds.sphere(com, factor * a_max)
            surface = ds.surface(smol_sp, ('gas', 'metallicity3'), Z_CRIT)
            clump_area = max(surface.surface_area.in_units('pc**2'), clump_area)
        except YTSphereTooSmall:
            factor = factor * SMOL_SP_INC_FACT
            continue
    ratio = (clump_volume / (clump_area * a_max)) * (3 * np.sqrt(2/5))
    return ("Clump volume to surface area ratio: %.5f", ratio)
add_clump_info("ratio", _ratio)

def _alpha(clump):
    a_max = clump.info['a_max'][1]
    div2 = clump.data.quantities.total_quantity(('gas', 'metal_laplacian'))
    alpha = div2 / a_max
    return ("Clump dimensionless metallicity gradient over boundary: %.3e", alpha)
add_clump_info("alpha", _alpha)
