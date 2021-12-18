import numpy as np
import os
import time

from neorl import DE  # JAYA,MFO,PSO,HHO
import openmc

start = time.time()

os.environ['OPENMC_CROSS_SECTIONS'] = '/home/super/nuclear_data/endfb71_hdf5/cross_sections.xml'

assm_width = 22

def pwr_assembly(void_loc_x=np.array([]), void_loc_y=np.array([])):
    """Create a PWR assembly model.

    This model is a reflected 17x17 fuel assembly from the the `BEAVRS
    <http://crpg.mit.edu/research/beavrs>`_ benchmark. The fuel is 2.4 w/o
    enriched UO2 corresponding to a beginning-of-cycle condition. Note that the
    number of particles/batches is initially set very low for testing purposes.

    Returns
    -------
    model : openmc.model.Model
        A PWR assembly model

    """

    model = openmc.model.Model()

    # Define materials
    fuel = openmc.Material(name='Fuel')  # UO2 fuel 
    fuel.set_density('g/cm3', 10.29769)
    fuel.add_nuclide('U234', 4.4843e-6)
    fuel.add_nuclide('U235', 5.5815e-4)
    fuel.add_nuclide('U238', 2.2408e-2)
    fuel.add_nuclide('O16', 4.5829e-2)

    clad = openmc.Material(name='Cladding') # Zr cladding
    clad.set_density('g/cm3', 6.55)
    clad.add_nuclide('Zr90', 2.1827e-2)
    clad.add_nuclide('Zr91', 4.7600e-3)
    clad.add_nuclide('Zr92', 7.2758e-3)
    clad.add_nuclide('Zr94', 7.3734e-3)
    clad.add_nuclide('Zr96', 1.1879e-3)

    hot_water = openmc.Material(name='Hot borated water') # Borated water as moderator
    hot_water.set_density('g/cm3', 0.740582)
    hot_water.add_nuclide('H1', 4.9457e-2)
    hot_water.add_nuclide('O16', 2.4672e-2)
    hot_water.add_nuclide('B10', 8.0042e-6)
    hot_water.add_nuclide('B11', 3.2218e-5)
    hot_water.add_s_alpha_beta('c_H_in_H2O')

    # Define the materials file.
    model.materials = (fuel, clad, hot_water)

    # Instantiate ZCylinder surfaces
    fuel_or = openmc.ZCylinder(x0=0, y0=0, r=0.75, name='Fuel OR')
    clad_or = openmc.ZCylinder(x0=0, y0=0, r=0.85, name='Clad OR')

    # Create boundary planes to surround the geometry
    pitch = assm_width
    min_x = openmc.XPlane(x0=-pitch/2, boundary_type='vacuum') # x/y boundary condition is 'void'
    max_x = openmc.XPlane(x0=+pitch/2, boundary_type='vacuum') 
    min_y = openmc.YPlane(y0=-pitch/2, boundary_type='vacuum')
    max_y = openmc.YPlane(y0=+pitch/2, boundary_type='vacuum')

    # Create a fuel pin universe
    fuel_pin_universe = openmc.Universe(name='Fuel Pin') 
    fuel_cell = openmc.Cell(name='fuel', fill=fuel, region=-fuel_or)
    clad_cell = openmc.Cell(name='clad', fill=clad, region=+fuel_or & -clad_or)
    hot_water_cell = openmc.Cell(name='hot water', fill=hot_water, region=+clad_or)
    fuel_pin_universe.add_cells([fuel_cell, clad_cell, hot_water_cell])

    # Create a void pin universe
    void_pitch = pitch/11
    min_void_x = openmc.XPlane(x0=-void_pitch/2)
    max_void_x = openmc.XPlane(x0=+void_pitch/2) 
    min_void_y = openmc.YPlane(y0=-void_pitch/2)
    max_void_y = openmc.YPlane(y0=+void_pitch/2)
    void_universe = openmc.Universe(universe_id=900001, name='void pin')
    void_cell = openmc.Cell(name='void', region=+min_void_x & -max_void_x & \
        +min_void_y & -max_void_y)
    void_universe.add_cell(void_cell)

    # Create fuel assembly Lattice
    assembly = openmc.RectLattice(name='Fuel Assembly')
    assembly.pitch = (pitch/11, pitch/11)
    assembly.lower_left = (-pitch/2, -pitch/2)

    # Create 17x17 array of universes
    assembly.universes = np.tile(fuel_pin_universe, (11, 11))
    assembly.universes[void_loc_x, void_loc_y] = void_universe

    # Create root Cell
    root_cell = openmc.Cell(name='root cell', fill=assembly)
    root_cell.region = +min_x & -max_x & +min_y & -max_y

    # Create root Universe
    model.geometry.root_universe = openmc.Universe(name='root universe')
    model.geometry.root_universe.add_cell(root_cell)

    model.settings.batches = 100 # Total simulation batches 
    model.settings.inactive = 30 # Inactive simulation batches
    model.settings.particles = 10000 # Number of simulation particles
    model.settings.source = openmc.Source(space=openmc.stats.Box(
        [-pitch/2, -pitch/2, -1], [pitch/2, pitch/2, 1], only_fissionable=True)) # define source

    plot = openmc.Plot() # 
    plot.origin = (0.0, 0.0, 0)
    plot.width = (assm_width, assm_width)
    plot.pixels = (300, 300)
    plot.color_by = 'material'
    model.plots.append(plot)

    return model


## call NEORL to find the optimal geometry config to max k-eff ## 
# Define the fitness function
def FIT(arr):

    total_pin = 121 # the assembly has 121 pins totally
    fuel_limit = 61 # limit of fuel units

    list_x, list_y = [], [] # store locations of void pin
    for idx, val in enumerate(arr):
        row, col = idx//11, idx%11 # pin location-(row, col) in the assembly
        if val == 0: # void pin
            list_x.append(row)
            list_y.append(col)
    lx = np.array(list_x)
    ly = np.array(list_y)
    model = pwr_assembly(void_loc_x=lx, void_loc_y=ly)
    
    # in case the program interrupted due to all neutrons leak, use try-except
    try: 
        result_r = model.run(output=True, threads=128)

        sp = openmc.StatePoint(result_r)
        k_combined = sp.k_combined
        k_combined_nom = k_combined.nominal_value
        k_combined_stddev = k_combined.std_dev

        # penalty of over-use fuel
        penalty = -1e5
        used_fuel = total_pin - len(lx)
        if used_fuel > fuel_limit: return_val = k_combined_nom + penalty
        else: return_val = k_combined_nom

    except:  
        print('All neutrons leak')
        return 0.0

    return np.round(return_val,5)

# Setup the parameter space(enrichment of U belongs to[0,4.0])
nx=121
BOUNDS={}
for i in range(1,nx+1):
    BOUNDS['x'+str(i)]=['int', 0, 1]


# setup and evolute DE
de=DE(mode='max', bounds=BOUNDS, fit=FIT, npop=50, F=0.5, CR=0.3,  ncores=1, seed=100)
x_best, y_best, de_hist=de.evolute(ngen=400, x0=None, verbose=1)
print('---DE Results---', )
print('x:', x_best)
print('y:', y_best)
print('DE History:\n', de_hist)
end = time.time()
running_time = end - start
print('running time:\n', running_time)




