import numpy as np
import os
import time

import openmc
from smr.materials import materials
from smr.plots import core_plots
from smr.surfaces import lattice_pitch, bottom_fuel_stack, top_active_core
from smr.core import core_geometry
from smr import inlet_temperature
from neorl import JAYA, DE, MFO

start = time.time()

# os.environ['OPENMC_CROSS_SECTIONS'] = '/home/adt/neorl/openmc/endfb71_hdf5/cross_sections.xml'

def smr(CRs):
    '''Build the SMR model

    param CRs: list of insertion depths for bank A,B,C and D 
    '''
    multipole = False
    rings, axial, depleted = 10, 196, False

    model = openmc.model.Model()

    model.geometry = core_geometry(rings, axial, depleted, CRs)

    all_materials = model.geometry.get_all_materials()
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


## call NEORL to find the optimal enrichment ## 
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

# Setup the parameter space(enrichment of U belongs to[0,4.0])
nx=4
BOUNDS={}
for i in range(1,nx+1):
    BOUNDS['x'+str(i)]=['float', 0, 359.634]

# use JAYA to find the optimal U enrichment
jaya=JAYA(mode='min', bounds=BOUNDS, fit=FIT, npop=10, ncores=1, seed=100)
x_best, y_best, jaya_hist=jaya.evolute(ngen=30, verbose=1)
print('---JAYA Results---', )
print('x:', x_best)
print('y:', y_best)
print('JAYA History:\n', jaya_hist)
end = time.time()
running_time = end - start
print('running time:\n', running_time)

# use DE to find the optimal U enrichment
# de=DE(mode='min', bounds=BOUNDS, fit=FIT, npop=10, F=0.5, CR=0.3,  ncores=1, seed=100)
# x_best, y_best, de_hist=de.evolute(ngen=30, verbose=1)
# print('---DE Results---', )
# print('x:', x_best)
# print('y:', y_best)
# print('DE History:\n', de_hist)
# end = time.time()
# running_time = end - start
# print('running time:\n', running_time)

# # use MFO to find the optimal U enrichment
# mfo=MFO(mode='min', bounds=BOUNDS, fit=FIT, nmoths=5, ncores=1, seed=100)
# x_best, y_best, mfo_hist=mfo.evolute(ngen=5, verbose=1)
# print('---MFO Results---', )
# print('x:', x_best)
# print('y:', y_best)
# print('MFO History:\n', mfo_hist)
# end = time.time()
# running_time = end - start
# print('running time:\n', running_time)