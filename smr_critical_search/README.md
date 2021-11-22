## SMR criticality search

**Small modular reactors** (**SMR**s) are [nuclear fission reactors](https://en.wikipedia.org/wiki/Nuclear_fission_reactor) that are a fraction of the size of conventional reactors. They are designed to be manufactured at a plant and transported to a site to be installed. Modular reactors will reduce on-site construction, increase containment efficiency, and enhance safety. The greater safety will come via the use of [passive safety](https://en.wikipedia.org/wiki/Passive_nuclear_safety) features that operate without human intervention[1]. 

**Criticality**. In a nuclear reactor, the neutron population at any instant is a function of the rate of neutron production (due to fission processes) and the rate of neutron losses (due to non-fission absorption mechanisms and leakage from the system). When a reactor’s neutron population remains steady from one generation to the next (creating as many new neutrons as are lost), the fission chain reaction is self-sustaining and the reactor's condition is referred to as "critical"[2].

In this case, we aim to find how to set the control rod banks in the SMR to make it critical. 

## SMR criticality search

### SMR Benchmark Model

The SMR OpenMC model refers to the [mit-crpg/ecp-benchmarks](https://github.com/mit-crpg/ecp-benchmarks), which is a 3D Small Modular Reactor (SMR) model that roughly mimics the design of the NuScale reactor. 

Below is an axial x-z slice of the SMR core, we could see various kinds of assembly in it, such as control rod assembly, instrumentation tube assembly, and burnable absorber assembly, etc.  

<img src='.\pics\axial_xz_slice.png' width='70%' height='70%' alt='axial_xz_slice' align=middle />

### Core Structure and Assembles

In this case, we design the SMR with 4 regulating groups of control rod banks (bank A, B, C, and D), which is shown as follows. The most outside of the core are the reflector, and the inner ring is the in-core instrumentation assembles, then the regulation control rod back groups and the burnable absorber assemblies, finally, in the center is the in-core instrumentation assembly.  

<img src='.\pics\radial_xy_slice_demonstration.png' width='120%' height='120%' alt='axial_xz_slice' align=middle />



### Search Objective 

To achieve the critical condition, we have to find a combination insertion depth for the 4 groups of the control rod banks. For every group of the control rod bank, the insertion depth can vary from 0 to 359.634 centimeters, thus, we have to determine the insertion depth for each control rod bank group from this range. In the reality, the insertion depth is discrete (about 228 steps in SMR), but in this case, we regard it as a **continuous** problem, which provides with more accurate solution.

### Method

Here we call the [JAYA](https://neorl.readthedocs.io/en/latest/modules/jaya.html?highlight=JAYA) and [Deferential Evolutionary](https://neorl.readthedocs.io/en/latest/modules/de.html)  algorithms from NEORL to find an excellent solution.

Some key parameters are set as:

```python
# OpenMC
Monte-Carlo simulation total batches = 150
Monte-Carlo simulation inactive batches = 30
Monte-Carlo simulation particles = 10000

# NEORL: keep the same parameters for every algorithm to compare their performance
Population of the optimization algorithm: 10
Generation of the optimization algoritm: 30
Random seed: 100
```

Follows the complete code:

```python
import numpy as np
import os
import time

import openmc
from smr.materials import materials
from smr.plots import core_plots
from smr.surfaces import lattice_pitch, bottom_fuel_stack, top_active_core
from smr.core import core_geometry
from smr import inlet_temperature
from neorl import JAYA, DE

start = time.time()

# os.environ['OPENMC_CROSS_SECTIONS'] = '/home/adt/neorl/openmc/endfb71_hdf5/cross_sections.xml'

def smr(CRs):
    '''Build the SMR model
    param CRs: list of insertion depths for bank A,B,C and D 
    '''
    multipole = False
    rings, axial, depleted = 10, 196, False

    model = openmc.model.Model()
    
    model.geometry = core_geometry(rings, axial, depleted, CRs) # Geometry

    all_materials = model.geometry.get_all_materials() # Materials
    model.materials = openmc.Materials(all_materials.values())

    # Construct uniform initial source distribution over fissionable zones
    lower_left = [-7.*lattice_pitch/2., -7.*lattice_pitch/2., bottom_fuel_stack]
    upper_right = [+7.*lattice_pitch/2., +7.*lattice_pitch/2., top_active_core]
    source = openmc.source.Source(space=openmc.stats.Box(lower_left, upper_right))
    source.space.only_fissionable = True

    model.settings = openmc.Settings()
    model.settings.batches = 150
    model.settings.inactive = 30
    model.settings.particles = 10000
    model.settings.output = {'tallies': False, 'summary': False}
    model.settings.source = source
    model.settings.sourcepoint_write = False
    model.settings.temperature = {
        'default': inlet_temperature,
        'method': 'interpolation',
        'range': (300.0, 1500.0),
    }
    if multipole:
        model.settings.temperature['multipole'] = True
        model.settings.temperature['tolerance'] = 1000

    return model


## call NEORL for optimization ## 
# Define the fitness function
def FIT(x):

    model = smr(CRs=x)
    result_r = model.run(output=False, threads=16)
    sp = openmc.StatePoint(result_r)
    k_combined = sp.k_combined
    k_combined_nom = k_combined.nominal_value
    k_combined_stddev = k_combined.std_dev
    return_val = abs(k_combined_nom - 1.0) 

    return return_val

# Setup the problem search space
nx=4
BOUNDS={}
for i in range(1,nx+1):
    BOUNDS['x'+str(i)]=['float', 0, 359.634] 

# use JAYA to find the optimal parameters
jaya=JAYA(mode='min', bounds=BOUNDS, fit=FIT, npop=10, ncores=1, seed=100)
x_best, y_best, jaya_hist=jaya.evolute(ngen=30, verbose=1)
print('---JAYA Results---', )
print('x:', x_best)
print('y:', y_best)
print('JAYA History:\n', jaya_hist)


# use DE to find the optimal parameters
de=DE(mode='min', bounds=BOUNDS, fit=FIT, npop=10, F=0.5, CR=0.3,  ncores=1, seed=100)
x_best, y_best, de_hist=de.evolute(ngen=30, verbose=1)
print('---DE Results---', )
print('x:', x_best)
print('y:', y_best)
print('DE History:\n', de_hist)
```

### Results

```python
# JAYA
------------------------ JAYA Summary --------------------------
Best individual (x) found: [281.16951924 334.89014362  90.18330327 299.00796372]
Best individual (y) found: 6.633297668923177e-06
JAYA History:
[0.0027812217461274935, 0.0011909121198945272, 0.0011909121198945272, 0.0011909121198945272, 0.0011909121198945272, 5.672083497421099e-05, 5.672083497421099e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05,
2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 2.7849536240354134e-05, 6.633297668923177e-06, 6.633297668923177e-06, 6.633297668923177e-06, 6.633297668923177e-06, 6.633297668923177e-06]

# DE 
------------------------ DE Summary --------------------------
Best individual (x) found: [32.51209409591024, 359.634, 141.91959001818773, 355.32843192667775]
Best fitness (y) found: 1.9786120794562656e-05
DE History:
[0.0003842112016596566, 0.0003842112016612109, 0.00038421120166332035, 0.0003842112016559929, 0.00038421120166076683, 0.0003842112016619881, 0.0003842112016655408, 0.0003361088373750043, 0.0003361088373690091, 0.00033610883737167363, 0.0003361088373661225, 0.0002679615958958159, 0.0002679615958991466, 0.0002679615958931514, 0.0002679615958942616, 0.0002679615958902648, 0.0002679615958951498, 0.00026796159589603796, 0.00026796159589204116, 0.0002679615958951498, 0.00012359618140034279, 0.00012359618140567186, 0.0001235961814058939, 0.0001235961813810249, 0.0001235961813910169, 1.9786120792675277e-05, 1.9786120793452433e-05, 1.9786120789233586e-05, 1.9786120789566652e-05, 1.9786120794562656e-05]

```

<img src='.\pics\jaya and de compare.png' width='60%' height='60%' align=middle />

After searching for 30 generations: 

* JAYA found the best combination of control rod insertion depths is [281.16951924 334.89014362  90.18330327 299.00796372], and its  absolute error is 6.63e-06(with the goal k-eff=1.0)
* DE found the best combination is [32.51209409591024, 359.634, 141.91959001818773, 355.32843192667775], whose absolute error is 1.98e-05

JAYA has a high search efficiency during the search process, for the result in generation 8 is already pretty good.  

## Reference

[1] https://en.wikipedia.org/wiki/Small_modular_reactor

[2] https://en.wikipedia.org/wiki/Nuclear_reactor_physics#Criticality

[3] MIT CRPG smr benchmark: https://github.com/mit-crpg/ecp-benchmarks

