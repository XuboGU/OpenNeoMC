# Assembly's max k-eff - a discrete problem with constraint 

## Assembly Model

The assembly model refers from the reference paper [*Gradient Informed Design Optimization of Select Nuclear Systems*](https://www.tandfonline.com/doi/abs/10.1080/00295639.2021.1987133)[1].  

**Geometry & Materials**

The assembly model is an 11 X 11 voxelated square grid  with a total width of around 1 meter. Each voxel wad filled with material 1 (void) or material 2 (a fuel-moderator mixture). The fuel-moderator mixture is the traditional UN TRISO fuel embeded in a SiC shell, and the moderator is yttrium hydride[2]. 

The fuel-moderator materials amount is limited: for example, if there only have 57 units of fuel/moderator, then one can set at most 57 voxel as material 2. 

**Boundary condition**

* In X and Y directions, the boundary condition is set as 'void' (i.e.: 'vaccum' in openMC)
* In Z direction, the boundary condition is set as 'reflective' (same in openMC)

**Optimal geometric configuration**

The reason the authors set X and Y directions of the model as void boundary conditions is that the optimal geometric configuration can be predictable - a circle, wherein the center of geometry should be fuel pins while edging area the void pin. This circle geometry can be a benchmark. The following picture is an optimal geometic configuration for 47 fuel pins.

<center>    
    <img style="border-radius: 0.3125em;    
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"     
    src=".\pics\assm_max_paper.png">    
    <br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    
    display: inline-block;    
    color: #999;    
   	padding: 2px;">Optimal geometric configuration (47 fuel pins)[1]</div> 
</center>



## Assembly Model Optimization 

**Problem Definition**

The optimization objective of this case is to find the optimal geometric configuration to maximize k-eff using limited amount of fuel-moderator materials in this 11 X 11 geometry. To be more specific, one needs to determine which material (1 or 2) to fill in every voxel of the geometry. 

In this case, we use up to 61 fuel pins to fill the assembly and need to find the optimal geometric configuration to make the k-eff max. Equation (1) is the math definition of the problem. Where $(p_i = 0)$ represents the void pin and $(p_i = 1)$ represents the fuel pin.

$$
\begin{equation} 
\begin{aligned} 
max_{\vec p}f(\vec p), \\
subject \ to,\\ 
\sum_ip_i\leq61, \\
where \ p_i = 0 \ or \ p_i = 1 
\end{aligned}
\tag 1
\end{equation}
$$



**Possible combinations**

The possible combinations of filling the assembly can be calculated by eqution (2).
$$
C_{121}^{61} \approx 7.973^{116} \tag 2
$$
This is a quite large number. We need to find a good enough combination among so many possible choices.

### Method

Here we call the [Defferential Evolution](https://neorl.readthedocs.io/en/latest/modules/de.html#differential-evolution-de) algorithm to find the optimal geometric configuration.

```python
import numpy as np
import os
import time

from neorl import DE  
import openmc

start = time.time()

# cross_sections config  
# os.environ['OPENMC_CROSS_SECTIONS'] = 'PATH/endfb71_hdf5/cross_sections.xml'

assm_width = 22 # the total width and lengh of assembly is 22 cm

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

    plot = openmc.Plot() # plot the assembly geometry
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
        result_r = model.run(output=True, threads=128) # output the result and run in 128 threads

        sp = openmc.StatePoint(result_r) 
        k_combined = sp.k_combined # combined k-eff
        k_combined_nom = k_combined.nominal_value # nominal value of k-eff
        k_combined_stddev = k_combined.std_dev # standard deviation of k-eff

        # the penalty of fuel units exceed the limit
        penalty = -1e5
        used_fuel = total_pin - len(lx) 
        if used_fuel > fuel_limit: return_val = k_combined_nom + penalty # excced the limit 
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

```

### Results

```python
DE step 20000/20000, F=0.5, CR=0.3, Ncores=1
Best fitness (y) found: 0.20112
Best individual (x) found: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
```

The following figure showes the converge curve of Differential Evolution. It showes that the algorithm find the best combination at around generation 300.  The best k-eff found is 0.20112. 

<img src=".\pics\de_curve.png" alt="converge curve" style="zoom:80%;" />


The geometry of the best solution is shown in the following figure, where a good circle is found. This result coincides with the benchmark.

<center>    
    <img style="border-radius: 0.3125em;    
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"     
    src=".\pics\assm_max.png">    
    <br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    
    display: inline-block;    
    color: #999;    
   	padding: 2px;">Optimal geometric configuration found by DE (61 fuel pins)</div> 
</center>

#### Further exploration #### 

* Physical-informed optimization: introduce the plane symmetry property of assembly design to reduce the search space.
* More evolutionary algorithms will be added

## Reference

[1] John Pevey, Briana Hiscox, Austin Williams, Ondřej Chvála, Vladimir Sobes & J. Wesley Hines (2021) Gradient-Informed Design Optimization of Select Nuclear Systems, Nuclear Science and Engineering

[2] Transformational Challenge Reactor Preconceptual Core Design Studie: https://www.osti.gov/servlets/purl/1651338

