# -*- encoding: utf-8 -*-
'''
@File    :   pin_cell_test.py
@Time    :   2021/10/19 09:58:33
@Author  :   Xubo GU 
@Email   :   guxubo@alumni.sjtu.edu.cn
'''
# here put the import lib

import openmc
from neorl import JAYA
import os
import time

start = time.time()
os.environ['OPENMC_CROSS_SECTIONS'] = '/home/adt/neorl/openmc/endfb71_hdf5/cross_sections.xml'

## OpenMC pwr pin cell case ##  
def pwr_pin_cell(U_enrich):
    """Create a PWR pin-cell model.

    This model is a single fuel pin with 2.4 w/o enriched UO2 corresponding to a
    beginning-of-cycle condition and borated water. The specifications are from
    the `BEAVRS <http://crpg.mit.edu/research/beavrs>`_ benchmark. Note that the
    number of particles/batches is initially set very low for testing purposes.

    Returns
    -------
    model : openmc.model.Model
        A PWR pin-cell model

    """
    model = openmc.model.Model()

    # Define materials.
    fuel = openmc.Material(name='UO2')
    fuel.set_density('g/cm3', 10.29769)
    fuel.add_element('U', 1.0, enrichment=U_enrich)
    fuel.add_element('O', 2.0)

    clad = openmc.Material(name='Zircaloy')
    clad.set_density('g/cm3', 6.55)
    clad.add_nuclide('Zr90', 2.1827e-2)
    clad.add_nuclide('Zr91', 4.7600e-3)
    clad.add_nuclide('Zr92', 7.2758e-3)
    clad.add_nuclide('Zr94', 7.3734e-3)
    clad.add_nuclide('Zr96', 1.1879e-3)

    hot_water = openmc.Material(name='Hot borated water')
    hot_water.set_density('g/cm3', 0.740582)
    hot_water.add_nuclide('H1', 4.9457e-2)
    hot_water.add_nuclide('O16', 2.4672e-2)
    hot_water.add_nuclide('B10', 8.0042e-6)
    hot_water.add_nuclide('B11', 3.2218e-5)
    hot_water.add_s_alpha_beta('c_H_in_H2O')

    # Define the materials file.
    model.materials = (fuel, clad, hot_water)

    # Instantiate ZCylinder surfaces
    pitch = 1.26
    fuel_or = openmc.ZCylinder(x0=0, y0=0, r=0.39218, name='Fuel OR')
    clad_or = openmc.ZCylinder(x0=0, y0=0, r=0.45720, name='Clad OR')
    left = openmc.XPlane(x0=-pitch/2, name='left', boundary_type='reflective')
    right = openmc.XPlane(x0=pitch/2, name='right', boundary_type='reflective')
    bottom = openmc.YPlane(y0=-pitch/2, name='bottom',
                           boundary_type='reflective')
    top = openmc.YPlane(y0=pitch/2, name='top', boundary_type='reflective')

    # Instantiate Cells
    fuel_pin = openmc.Cell(name='Fuel', fill=fuel)
    cladding = openmc.Cell(name='Cladding', fill=clad)
    water = openmc.Cell(name='Water', fill=hot_water)

    # Use surface half-spaces to define regions
    fuel_pin.region = -fuel_or
    cladding.region = +fuel_or & -clad_or
    water.region = +clad_or & +left & -right & +bottom & -top

    # Create root universe
    model.geometry.root_universe = openmc.Universe(0, name='root universe')
    model.geometry.root_universe.add_cells([fuel_pin, cladding, water])

    model.settings.batches = 150
    model.settings.inactive = 30
    model.settings.particles = 10000
    model.settings.source = openmc.Source(space=openmc.stats.Box(
        [-pitch/2, -pitch/2, -1], [pitch/2, pitch/2, 1], only_fissionable=True))

    return model


## call NEORL to find the optimal enrichment ## 
# Define the fitness function
def FIT(x):
    """Find the enrichment of U to achieve the wanted k-eff = 1.10
    """
    model = pwr_pin_cell(U_enrich=x[0])
    result_r = model.run(output=False)
    sp = openmc.StatePoint(result_r)
    k_combined = sp.k_combined
    k_combined_nom = k_combined.nominal_value
    k_combined_stddev = k_combined.std_dev
    return_val = abs(k_combined_nom - 1.10) # fitness = k_combined - k_target 

    return return_val

# Setup the parameter space(enrichment of U belongs to[0,4.0])
nx=1
BOUNDS={}
for i in range(1,nx+1):
    BOUNDS['x'+str(i)]=['float', 0.0, 4.0]

# use JAYA to fine the optimal U enrichment
jaya=JAYA(mode='min', bounds=BOUNDS, fit=FIT, npop=10, ncores=1, seed=100)
x_best, y_best, jaya_hist=jaya.evolute(ngen=5, verbose=1)
print('---JAYA Results---', )
print('x:', x_best)
print('y:', y_best)
print('JAYA History:\n', jaya_hist)

end = time.time()
running_time = end - start
print('running time:\n', running_time)
