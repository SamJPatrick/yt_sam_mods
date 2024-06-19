from functools import reduce
import yt
from unyt import unyt_quantity, unyt_array
from unyt import kb, G, mh, Msun, pc
import numpy as np



def multiply_fields(ds, *quantities):
    def _field_product(field, data):
        return reduce(lambda x, y: x * y, [data[quantity] for quantity in quantities])
    return _field_product


def generate_vol_laplacian(ds, quantity):
    hassian = []
    grad_fields = ds.add_gradient_fields(quantity)
    for grad_field in grad_fields[:-1]:
        hassian.append(ds.add_gradient_fields(grad_field)[:-1])
    def _vol_laplacian(field, data):
        div2 = sum([data[hassian[i][i]] * data[('gas', 'cell_volume')] for i in range (0, 3)])
        return div2
    return _vol_laplacian


def generate_field_overlap(ds, field_1, field_2):
    grads_1 = np.array(ds.add_gradient_fields(field_1))
    grads_2 = np.array(ds.add_gradient_fields(field_2))
    grad_prod = np.dot(grads_1[:-1], grads_2[:-1].reshape(2,3)) / (grads_1[3] * grads_2[3])
    def _field_overlap(field, data):
        overlap = data[grad_prod] * data[('gas', 'cell_volume')]
        return overlap
    return _field_overlap


def generate_divergence(ds, quantity):
    for ax in 'xyz':
        ds.add_field((f'{quantity[0]}', f'{quantity[1]}_vel_{ax}'), \
                      function= multiply_fields(ds, ('gas', f'velocity_{ax}'), quantity), \
                      sampling_type= 'local', units='', take_log=False)
    divs = [ds.add_gradient_fields((f'{quantity[0]}', f'{quantity[1]}_vel_{ax}'))[i] \
            for i, ax in enumerate('xyz')]
    def _divergence(field, data):
        return sum([data[div_comp] for div_comp in divs])
    return _divergence




def get_jeans_global(sp):
    temp_ave = sp.quantities.weighted_average_quantity(('gas', 'temperature'), ('gas', 'cell_mass'))
    mu_ave = sp.quantities.weighted_average_quantity(('gas', 'mean_molecular_weight'), ('gas', 'cell_mass'))
    rho_ave = sp.quantities.weighted_average_quantity(('gas', 'density'), ('gas', 'cell_mass'))
    mj_const = ((5.0 * kb) / (G * mh))**1.5 * (3.0 / (4.0 * np.pi))**0.5
    jeans_mass =  mj_const * (temp_ave / mu_ave)**1.5 * rho_ave**(-0.5)
    return jeans_mass


def get_vkep(sp):
    mass = sp.quantities.total_quantity(('gas', 'cell_mass'))
    v_kep = (((G * mass) / sp.radius)**(1/2)).in_units('km/s')
    return v_kep
