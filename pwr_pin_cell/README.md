# PWR Pin-cell test

## PWR Pin-cell Model

PWR pin-cell model is a single fuel pin with UO2 corresponding to a beginning-of-cycle condition and borated water. The specifications are from  the [BEAVRS](http://crpg.mit.edu/research/beavrs) benchmark.

The pwr_pin_cell model is built-in module in Openmc, you can see here [openmc.examples.pwr_pin_cell](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.examples.pwr_pin_cell.html#openmc.examples.pwr_pin_cell)

## PWR Pin-cell Optimization 

### Optimization Goal

PWR pin-cell model is a single fuel pin with UO2. You boss told you that 'Well, I want to get a beginning k-eff of 1.1 from the single pin-cell, please tell me what is the enrichment of UO2'.

### Method

Here we call the [JAYA](https://neorl.readthedocs.io/en/latest/modules/jaya.html?highlight=JAYA) evolutionary algorithm to find the optimal UO2 enrichment.

* Note that the number of particles/batches is initially set very low for testing purposes.

```python
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

## Configure enviromental variable here ## 
# os.environ['OPENMC_CROSS_SECTIONS'] = '/PATH/cross_sections.xml'

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

    # Define materials
    fuel = openmc.Material(name='UO2') # 'UO2' fuel
    fuel.set_density('g/cm3', 10.29769) # add fuel's density infomation 
    fuel.add_element('U', 1.0, enrichment=U_enrich) # add Uranium element
    fuel.add_element('O', 2.0) # add Oxygen element

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
    pitch = 1.26 # in `cm`
    fuel_or = openmc.ZCylinder(x0=0, y0=0, r=0.39218, name='Fuel OR') # fuel outer radius
    clad_or = openmc.ZCylinder(x0=0, y0=0, r=0.45720, name='Clad OR') # clad outer radius
    left = openmc.XPlane(x0=-pitch/2, name='left', boundary_type='reflective') # left boundary
    right = openmc.XPlane(x0=pitch/2, name='right', boundary_type='reflective') # right boundary
    bottom = openmc.YPlane(y0=-pitch/2, name='bottom',
                           boundary_type='reflective') # bottom boundary
    top = openmc.YPlane(y0=pitch/2, name='top', boundary_type='reflective') # top boundary

    # Instantiate Cells
    fuel_pin = openmc.Cell(name='Fuel', fill=fuel) # fuel pin cell
    cladding = openmc.Cell(name='Cladding', fill=clad) # cladding cell
    water = openmc.Cell(name='Water', fill=hot_water) # water cell

    # Use surface half-spaces to define regions
    fuel_pin.region = -fuel_or # fuel pin cell's region
    cladding.region = +fuel_or & -clad_or # cladding cell'region 
    water.region = +clad_or & +left & -right & +bottom & -top # water cell's region

    # Create root universe
    model.geometry.root_universe = openmc.Universe(0, name='root universe') # define the root universe
    model.geometry.root_universe.add_cells([fuel_pin, cladding, water]) # add cells to the root universe

    model.settings.batches = 150 # OpenMC total calculating batches
    model.settings.inactive = 30 # the first 30 batches are inactive 
    model.settings.particles = 10000 # number of simulation particles
    model.settings.source = openmc.Source(space=openmc.stats.Box(
        [-pitch/2, -pitch/2, -1], [pitch/2, pitch/2, 1], only_fissionable=True)) # define the source

    return model


## call NEORL to find the optimal enrichment ## 
# Define the fitness function
def FIT(x):
    """Find the enrichment of U to achieve the wanted k-eff = 1.10
    """
    model = pwr_pin_cell(U_enrich=x[0])
    result_r = model.run(output=False)
    sp = openmc.StatePoint(result_r)
    k_combined = sp.k_combined # the combined k-eff
    k_combined_nom = k_combined.nominal_value # the nominal value of k-eff
    k_combined_stddev = k_combined.std_dev # the standard deviation of k-eff
    return_val = abs(k_combined_nom - 1.10) # fitness = k_combined - k_target 

    return return_val

# Setup the parameter space(enrichment of U belongs to[0,4.0])
nx=1
BOUNDS={}
for i in range(1,nx+1):
    BOUNDS['x'+str(i)]=['float', 0.0, 4.0]

# use JAYA to find the optimal U enrichment
jaya=JAYA(mode='min', bounds=BOUNDS, fit=FIT, npop=10, ncores=1, seed=100)
x_best, y_best, jaya_hist=jaya.evolute(ngen=5, verbose=1)
print('---JAYA Results---', )
print('x:', x_best)
print('y:', y_best)
print('JAYA History:\n', jaya_hist)

end = time.time()
running_time = end - start
print('running time:\n', running_time)

```

### Results

The program runs around 3 minutes(may vary depending on your CPU), and the results are:

```python
------------------------ JAYA Summary --------------------------
Best fitness (y) found: 0.0015497217274231812
Best individual (x) found: [2.01355604]
--------------------------------------------------------------
---JAYA Results---
x: [2.01355604]
y: 0.0015497217274231812
JAYA History:
 [0.018311916874464318, 0.0017114252626817539, 0.0017114252626817539, 0.0017114252626817539, 0.0015497217274231812]
running time:
 155.2281835079193
```



## Reference

JAYA Original paper: Rao, R. (2016). **Jaya**: A simple and new optimization algorithm for solving constrained and unconstrained optimization problems. International Journal of Industrial Engineering Computations, 7(1), 19-34.