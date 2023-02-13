import statistics
import re
import unyt
from unyt import unyt_array, unyt_quantity
import numpy as np




def yaml_to_unyt_quantity(yaml_quantity):
    exp = re.match(r'((-)?[0-9]*(\.[0-9]+)?(e(\+|-)?[0-9]*)?) ([A-Za-z]*)$', yaml_quantity)
    quantity = float(exp.group(1))
    units = exp.group(6)
    if (units == 'unitary'):
        return unyt_quantity(quantity, 'dimensionless')
    else :
        return unyt_quantity(quantity, units)


def yaml_to_unyt_array(yaml_array):
    exp = re.match(r'(\[((-)?[0-9]*(\.[0-9]+)*(, )?)*\]) ([A-Za-z]*)$', yaml_array)
    array = [float(x) for x in exp.group(1).strip('[]').split(', ')]
    units = exp.group(6)
    if (units == 'unitary'):
        return unyt_array(array, 'dimensionless')
    else :
        return unyt_array(array, units)
    

def unyt_quantity_to_yaml(quantity):
    if (not isinstance(quantity, unyt_quantity)):
            print("Error, object is not an instance of a unyt quantity!")
            quit()
    return f"{float(quantity.value)} {quantity.units}"


def unyt_array_to_yaml(quantity):
    if (not isinstance(quantity, unyt_array)):
            print("Error, object is not an instance of a unyt array!")
            quit()
    return f"{quantity.value.tolist()} {quantity.units}"


def unyt_norm(array):
    return unyt_quantity(np.linalg.norm(array.value), array.units)


def reunit(ds, val, units):
    if isinstance(val, unyt_quantity):
        func = ds.quan
    else:
        func = ds.arr
    return func(val.to(units).d, units)


def transpose_unyt(quant_arr):
    values, units = zip(*[(quantity.value.item(), quantity.units) for quantity in quant_arr])
    unit = statistics.mode(units)
    arr = unyt_array([quantity.in_units(unit).value.item() for quantity in quant_arr], unit)
    return arr
