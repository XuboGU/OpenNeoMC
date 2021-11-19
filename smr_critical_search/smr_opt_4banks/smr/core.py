"""Instantiate the main core lattice."""

import numpy as np

import openmc

from .materials import mats
from .surfaces import surfs, lattice_pitch
from .reflector import reflector_universes
from .assemblies import assembly_universes


def core_geometry(num_rings, num_axial, depleted, CRs):
    """Generate full core SMR geometry.

    Parameters
    ----------
    num_rings : int
        Number of annual regions in fuel
    num_axial : int
        Number of axial subdivisions in fuel
    depleted : bool
        Whether fuel should contain nuclides as though it were depleted

    Returns
    -------
    openmc.Geometry
        SMR full core geometry

    """
    assembly = assembly_universes(num_rings, num_axial, depleted, CRs)
    reflector = reflector_universes()

    # Construct main core lattice
    core = openmc.RectLattice(name='Main core')
    core.lower_left = (-9*lattice_pitch/2, -9*lattice_pitch/2)
    core.pitch = (lattice_pitch, lattice_pitch)
    universes = np.tile(reflector['solid'], (9, 9))

    universes[0, 2] = reflector['0,2']
    universes[0, 3] = reflector['0,3']
    universes[0, 4] = reflector['0,4']
    universes[0, 5] = reflector['0,5']
    universes[0, 6] = reflector['0,6']

    universes[1, 1] = reflector['1,1']
    universes[1, 2] = reflector['NW']
    universes[1, 3] = assembly['Assembly (3.1%) instr']
    universes[1, 4] = assembly['Assembly (2.4%) CR B']
    universes[1, 5] = assembly['Assembly (3.1%) instr']
    universes[1, 6] = reflector['NE']
    universes[1, 7] = reflector['1,7']

    universes[2, 0] = reflector['2,0']
    universes[2, 1] = reflector['NW']
    universes[2, 2] = assembly['Assembly (3.1%) instr']
    universes[2, 3] = assembly['Assembly (2.4%) CR D']
    universes[2, 4] = assembly['Assembly (3.1%) 16BA']
    universes[2, 5] = assembly['Assembly (2.4%) CR C']
    universes[2, 6] = assembly['Assembly (3.1%) instr']
    universes[2, 7] = reflector['NE']
    universes[2, 8] = reflector['2,8']

    universes[3, 0] = reflector['3,0']
    universes[3, 1] = assembly['Assembly (3.1%) instr']
    universes[3, 2] = assembly['Assembly (2.4%) CR C']
    universes[3, 3] = assembly['Assembly (3.1%) 16BA']
    universes[3, 4] = assembly['Assembly (2.4%) CR A']
    universes[3, 5] = assembly['Assembly (3.1%) 16BA']
    universes[3, 6] = assembly['Assembly (2.4%) CR D']
    universes[3, 7] = assembly['Assembly (3.1%) instr']
    universes[3, 8] = reflector['3,8']

    universes[4, 0] = reflector['4,0']
    universes[4, 1] = assembly['Assembly (2.4%) CR B']
    universes[4, 2] = assembly['Assembly (3.1%) 16BA']
    universes[4, 3] = assembly['Assembly (2.4%) CR A']
    universes[4, 4] = assembly['Assembly (1.6%) instr']
    universes[4, 5] = assembly['Assembly (2.4%) CR A']
    universes[4, 6] = assembly['Assembly (3.1%) 16BA']
    universes[4, 7] = assembly['Assembly (2.4%) CR B']
    universes[4, 8] = reflector['4,8']

    universes[5, 0] = reflector['5,0']
    universes[5, 1] = assembly['Assembly (3.1%) instr']
    universes[5, 2] = assembly['Assembly (2.4%) CR D']
    universes[5, 3] = assembly['Assembly (3.1%) 16BA']
    universes[5, 4] = assembly['Assembly (2.4%) CR A']
    universes[5, 5] = assembly['Assembly (3.1%) 16BA']
    universes[5, 6] = assembly['Assembly (2.4%) CR C']
    universes[5, 7] = assembly['Assembly (3.1%) instr']
    universes[5, 8] = reflector['5,8']

    universes[6, 0] = reflector['6,0']
    universes[6, 1] = reflector['SW']
    universes[6, 2] = assembly['Assembly (3.1%) instr']
    universes[6, 3] = assembly['Assembly (2.4%) CR C']
    universes[6, 4] = assembly['Assembly (3.1%) 16BA']
    universes[6, 5] = assembly['Assembly (2.4%) CR D']
    universes[6, 6] = assembly['Assembly (3.1%) instr']
    universes[6, 7] = reflector['SE']
    universes[6, 8] = reflector['6,8']

    universes[7, 1] = reflector['7,1']
    universes[7, 2] = reflector['SW']
    universes[7, 3] = assembly['Assembly (3.1%) instr']
    universes[7, 4] = assembly['Assembly (2.4%) CR B']
    universes[7, 5] = assembly['Assembly (3.1%) instr']
    universes[7, 6] = reflector['SE']
    universes[7, 7] = reflector['7,7']

    universes[8, 2] = reflector['8,2']
    universes[8, 3] = reflector['8,3']
    universes[8, 4] = reflector['8,4']
    universes[8, 5] = reflector['8,5']
    universes[8, 6] = reflector['8,6']

    core.universes = universes

    root_univ = openmc.Universe(universe_id=0, name='root universe')

    # Cylinder filled with core lattice
    cell = openmc.Cell(name='Main core')
    cell.fill = core
    cell.region = \
        -surfs['core barrel IR'] & +surfs['lower bound'] & -surfs['upper bound']
    root_univ.add_cell(cell)

    # Core barrel
    cell = openmc.Cell(name='core barrel')
    cell.fill = mats['SS']      # stainless steel material
    cell.region = (+surfs['core barrel IR'] & -surfs['core barrel OR'] &
                   +surfs['lower bound'] & -surfs['upper bound'])
    root_univ.add_cell(cell)

    # Downcomer
    cell = openmc.Cell(name='downcomer')
    cell.fill = mats['H2O']  
    cell.region = (+surfs['core barrel OR'] & -surfs['RPV IR'] &
                   +surfs['lower bound'] & -surfs['upper bound'])
    root_univ.add_cell(cell)

    # Reactor pressure vessel
    cell = openmc.Cell(name='reactor pressure vessel')
    cell.fill = mats['CS'] # carbon steel material
    cell.region = (+surfs['RPV IR'] & -surfs['RPV OR'] &
                   +surfs['lower bound'] & -surfs['upper bound'])
    root_univ.add_cell(cell)

    # Return geometry
    return openmc.Geometry(root_univ)
