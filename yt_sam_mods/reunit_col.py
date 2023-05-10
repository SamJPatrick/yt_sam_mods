import sys
import re
import glob
import pandas as pd
import numpy as np
from unyt import unyt_quantity


FILE_TYPE = "dispersions"
FORMATS = {'Phi velocity dispersion': '{:.3e}', 'Radius (pc)': '{:.2f}'}
UNITS = {'Phi velocity dispersion': 'km/s'}

SCI_REGEX = r'^([\.0-9]+(e(\+|\-))?[0-9]+)\s([a-zA-Z/\*]+)$'
tab_files = glob.glob(''.join(['DD[0-9][0-9][0-9][0-9]_', FILE_TYPE, '.out']))
for filename in tab_files:
    table = pd.read_csv(filename, sep='\t')
    for column, col_unit in UNITS.items():
        col_index = table.columns.values.tolist().index(column)
        for i in range (table.shape[0]):
            val = table.iloc[i, col_index]
            quantity = float(re.search(SCI_REGEX, val).group(1))
            units = re.search(SCI_REGEX, val).group(4)
            reunit_val = unyt_quantity(quantity, units).in_units(col_unit).value
            table.iloc[i, col_index] = reunit_val.tolist()

    for column, col_format in FORMATS.items():
        try :
            table[column] = table[column].map(lambda x: col_format.format(x))
        except KeyError:
                print(f"Error, column {column} not entered in dictionary")
                continue
    table.to_csv(filename, sep='\t', index=False)
