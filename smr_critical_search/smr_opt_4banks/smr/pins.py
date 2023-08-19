"""Instantiate pin cell Cells and Universes for core model."""

from math import sqrt
import numpy as np
import openmc
from openmc.model import subdivide

from .materials import mats
from .surfaces import surfs, pellet_OR, bottom_fuel_stack, top_active_core


def make_pin(name, surfaces, materials, grid=None):
    """Construct a pin cell Universes with radially layered Cells.

    Parameters
    ----------
    name: str
        The string name to assign to the Universe and each of its Cells
    surfaces: Iterable of openmc.ZCylinder
        The surfaces between which Cells are constructed
        to comprise the radially layered pin cell Universe.
    materials: Iterable of openmc.Material
        The Materials used within each radial layer. This collection
        must be one unit longer than the collection of surfaces.
    grid: str, optional
        The type of grid spacer to wrap around the pin cell universe.
        Accepted types include 'bottom' and 'intermediate'.

    Returns
    -------
    universe: openmc.Universe
        The pin cell Universe
    """

    universe = openmc.Universe(name=name)

    # Create cell for interior of innermost ZCylinder
    cell_name = name + ' (0)'
    cell = openmc.Cell(name=cell_name, fill=materials[0], region=-surfaces[0])
    universe.add_cell(cell)

    # Create cells between two ZCylinders
    for i, (mat, surf) in enumerate(zip(materials[1:-1], surfaces[:-1])): # idx, (material, surface)
        cell_name = name + ' ({})'.format(i+1)
        cell = openmc.Cell(name=cell_name)
        cell.fill = mat
        cell.region = +surf & -surfaces[i+1]
        universe.add_cell(cell)

    # Create cell for exterior of outermost ZCylinder
    cell_name = name + ' (last)'
    cell = openmc.Cell(name=cell_name, fill=materials[-1], region=+surfaces[-1])
    universe.add_cell(cell)

    # Add spacer grid cells if specified
    if grid:
        cell.region &= surfs['rod grid box']

        cell_name = name + ' (grid)'
        cell = openmc.Cell(name=cell_name, region=~surfs['rod grid box'])

        if grid == 'bottom':
            cell.fill = mats['In']
        else:
            cell.fill = mats['Zr']

        universe.add_cell(cell)

    return universe


def make_stack(name, surfaces, universes):
    """Construct a Universe of axially stacked pin cell Universes.

    Parameters
    ----------
    name: str
        The string name to assign to the Universe and each of its Cells
    surfaces: Iterable of openmc.ZPlane
        A collection of axial surfaces between which pin cells are
        filled to comprise an axially stacked pin cell
    universes: Iterable of openmc.Universe
        The Universes used within each axial layer. This collection
        must be one unit longer than the collection of surfaces.

    Returns
    -------
    universe: openmc.Universe
        The pin cell Universe
    """

    universe = openmc.Universe(name=name)

    # Create cells for each axial segment
    for i, (univ, region) in enumerate(zip(universes, subdivide(surfaces))):
        cell_name = '{} ({})'.format(name, i)
        cell = openmc.Cell(name=cell_name, fill=univ, region=region) # fill materials
        universe.add_cell(cell)

    return universe


def make_pin_stack(name, zsurfaces, universes, boundary, fuel_fill):
    """Construct a Universe of axially stacked universes with a single inner fuel
    pin universe.

    Parameters
    ----------
    name: str
        The string name to assign to the Universe and each of its Cells
    zsurfaces: Iterable of openmc.ZPlane
        A collection of axial surfaces between which pin cells are
        filled to comprise an axially stacked pin cell
    universes: Iterable of openmc.Universe
        The Universes used within each axial layer. This collection
        must be one unit longer than the collection of surfaces.
    boundary : openmc.Surface
        Boundary between the fuel pin itself and everything outside (gap, clad,
        moderator)
    fuel_fill : openmc.Universe or openmc.Material
        Universe or material for (possibly subdivided) fuel

    Returns
    -------
    universe: openmc.Universe
        The pin cell Universe

    """

    universe = openmc.Universe(name=name)

    for i, (univ, r) in enumerate(zip(universes, subdivide(zsurfaces))):
        cell_name = '{} (o{})'.format(name, i)
        cell = openmc.Cell(name=cell_name, fill=univ, region=r & +boundary)
        universe.add_cell(cell)

    cell_name = '{} (i)'.format(name)
    cell = openmc.Cell(name=cell_name, fill=fuel_fill, region=-boundary)
    universe.add_cell(cell)

    return universe


def pin_universes(num_rings=10, num_axial=196, depleted=False, CRs=[0,0,0,0]):
    """Generate universes for SMR fuel pins.

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
    dict
        Dictionary mapping a universe name to a openmc.Universe object

    """
    fuel = 'depleted' if depleted else 'fresh'

    bank_top = 786.348
    bank_bot = 405.713



    surfs['bankD top'] = openmc.ZPlane(
    z0= bank_top - CRs[3], name='CR bank D top')
    surfs['bankD bot'] = openmc.ZPlane(
    z0= bank_bot - CRs[3], name='CR bank D bottom')

    surfs['bankC top'] = openmc.ZPlane(
    z0= bank_top - CRs[2], name='CR bank C top')
    surfs['bankC bot'] = openmc.ZPlane(
    z0= bank_bot - CRs[2], name='CR bank C bottom')

    surfs['bankB top'] = openmc.ZPlane(
    z0= bank_top - CRs[1], name='CR bank B top')
    surfs['bankB bot'] = openmc.ZPlane(
    z0= bank_bot - CRs[1], name='CR bank B bottom')

    surfs['bankA top'] = openmc.ZPlane(
    z0= bank_top - CRs[0], name='CR bank A top')
    surfs['bankA bot'] = openmc.ZPlane(
    z0= bank_bot - CRs[0], name='CR bank A bottom')

    # Create dictionary to store pin universes
    univs = {}

    # Dummy water cell
    cell = openmc.Cell(name='water pin', fill=mats['H2O'])
    univs['water pin'] = openmc.Universe(name='Empty water pin cell universe')
    univs['water pin'].add_cell(cell)

    # Guide tube pin cells
    univs['GT empty'] = make_pin(
        'GT empty',
        [surfs['GT IR'], surfs['GT OR']],
        [mats['H2O'], mats['Zr'], mats['H2O']])
    univs['GT empty grid (bottom)'] = make_pin(
        'GT empty grid (bottom)',
        [surfs['GT IR'], surfs['GT OR']],
        [mats['H2O'], mats['Zr'], mats['H2O']], grid='bottom')
    univs['GT empty grid (intermediate)'] = make_pin(
        'GT empty grid (intermediate)',
        [surfs['GT IR'], surfs['GT OR']],
        [mats['H2O'], mats['Zr'], mats['H2O']], grid='intermediate')
    univs['GT empty nozzle'] = make_pin(
        'GT empty nozzle',
        [surfs['GT IR'], surfs['GT OR']],
        [mats['H2O'], mats['Zr'], mats['H2O']])

    univs['GTd empty'] = make_pin(
        'GT empty at dashpot',
        [surfs['GT dashpot IR'], surfs['GT dashpot OR']],
        [mats['H2O'], mats['Zr'], mats['H2O']])
    univs['GTd empty grid (bottom)'] = make_pin(
        'GT empty at dashpot grid (bottom)',
        [surfs['GT dashpot IR'], surfs['GT dashpot OR']],
        [mats['H2O'], mats['Zr'], mats['H2O']], grid='bottom')
    univs['GTd empty grid (intermediate)'] = make_pin(
        'GT empty at dashpot grid (intermediate)',
        [surfs['GT dashpot IR'], surfs['GT dashpot OR']],
        [mats['H2O'], mats['Zr'], mats['H2O']], grid='intermediate')
    univs['GTd empty nozzle'] = make_pin(
        'GT empty nozzle',
        [surfs['GT dashpot IR'], surfs['GT dashpot OR']],
        [mats['H2O'], mats['Zr'], mats['H2O']])

    # Stack all axial pieces of guide tube together 
    # number: 20
    stack_surfs = [                    
        surfs['bot support plate'],
        surfs['top support plate'],
        surfs['top lower nozzle'],
        surfs['top lower thimble'],
        surfs['grid1bot'],
        surfs['grid1top'],
        surfs['dashpot top'],
        surfs['grid2bot'],
        surfs['grid2top'],
        surfs['grid3bot'],
        surfs['grid3top'],
        surfs['grid4bot'],
        surfs['grid4top'],
        surfs['top active core'],
        surfs['grid5bot'],
        surfs['grid5top'],
        surfs['top pin plenum'],
        surfs['top FR'],
        surfs['bot upper nozzle'],
        surfs['top upper nozzle']
    ]

    univs['GT empty'] = make_stack(
        'GT empty', surfaces=stack_surfs,
        universes=[univs['water pin'],
                   univs['water pin'],
                   univs['water pin'],
                   univs['GTd empty'],
                   univs['GTd empty'],
                   univs['GTd empty grid (bottom)'],
                   univs['GTd empty'],
                   univs['GT empty'],
                   univs['GT empty grid (intermediate)'],
                   univs['GT empty'],
                   univs['GT empty grid (intermediate)'],
                   univs['GT empty'],
                   univs['GT empty grid (intermediate)'],
                   univs['GT empty'],
                   univs['GT empty'],
                   univs['GT empty grid (intermediate)'],
                   univs['GT empty'],
                   univs['GT empty'],
                   univs['GT empty'],
                   univs['water pin'],
                   univs['water pin']])

    univs['GT empty instr'] = make_stack(
        'GT empty instr', surfaces=stack_surfs,
        universes=[univs['water pin'],
                   univs['water pin'],
                   univs['water pin'],
                   univs['GT empty'],
                   univs['GT empty'],
                   univs['GT empty grid (bottom)'],
                   univs['GT empty'],
                   univs['GT empty'],
                   univs['GT empty grid (intermediate)'],
                   univs['GT empty'],
                   univs['GT empty grid (intermediate)'],
                   univs['GT empty'],
                   univs['GT empty grid (intermediate)'],
                   univs['GT empty'],
                   univs['GT empty'],
                   univs['GT empty grid (intermediate)'],
                   univs['GT empty'],
                   univs['GT empty'],
                   univs['GT empty'],
                   univs['water pin'],
                   univs['water pin']])


    # Instrument tube pin cell
    univs['IT'] = make_pin(
        'IT',
        [surfs['IT IR'], surfs['IT OR'], surfs['GT IR'], surfs['GT OR']],
        [mats['Air'], mats['Zr'], mats['H2O'], mats['Zr'], mats['H2O']])
    univs['IT grid (bottom)'] = make_pin(
        'IT grid (bottom)',
        [surfs['IT IR'], surfs['IT OR'], surfs['GT IR'], surfs['GT OR']],
        [mats['Air'], mats['Zr'], mats['H2O'], mats['Zr'], mats['H2O']],
        grid='bottom')
    univs['IT grid (intermediate)'] = make_pin(
        'IT grid (intermediate)',
        [surfs['IT IR'], surfs['IT OR'], surfs['GT IR'], surfs['GT OR']],
        [mats['Air'], mats['Zr'], mats['H2O'], mats['Zr'], mats['H2O']],
        grid='intermediate')

    univs['IT nozzle'] = make_pin(
        'IT nozzle',
        [surfs['IT IR'], surfs['IT OR'], surfs['GT IR'], surfs['GT OR']],
        [mats['Air'], mats['Zr'], mats['H2O'], mats['Zr'], mats['SS']])
    univs['IT dashpot'] = make_pin(
        'IT dashpot',
        [surfs['IT IR'], surfs['IT OR']],
        [mats['Air'], mats['Zr'], mats['H2O']])

    # Stack all axial pieces of instrument tube together
    univs['IT stack'] = make_stack(
        'GT instr',
        surfaces=stack_surfs,
        universes=[univs['IT dashpot'],
                   univs['IT dashpot'],
                   univs['IT dashpot'],
                   univs['IT'],
                   univs['IT'],
                   univs['IT grid (bottom)'],
                   univs['IT'],
                   univs['IT'],
                   univs['IT grid (intermediate)'],
                   univs['IT'],
                   univs['IT grid (intermediate)'],
                   univs['IT'],
                   univs['IT grid (intermediate)'],
                   univs['IT'],
                   univs['IT'],
                   univs['IT grid (intermediate)'],
                   univs['IT'],
                   univs['IT'],
                   univs['IT'],
                   univs['IT dashpot'],
                   univs['water pin']])

    # Control rod pin cells
    univs['CR'] = make_pin( 
        'CR',
        [surfs['CP OR'], surfs['CR IR'], surfs['GT IR'], surfs['GT OR']],  # 'CP OR' - Control Poison OR
        [mats['AIC'], mats['Air'], mats['SS'], mats['H2O'], mats['Zr'], mats['H2O']]) # AIC - Ag-In-Cd control rod material
    univs['CR grid (bottom)'] = make_pin(
        'CR grid (bottom)',
        [surfs['CP OR'], surfs['CR IR'], surfs['GT IR'], surfs['GT OR']],
        [mats['AIC'], mats['Air'], mats['SS'], mats['H2O'], mats['Zr'], mats['H2O']],
        grid='bottom')
    univs['CR grid (intermediate)'] = make_pin(
        'CR grid (intermediate)',
        [surfs['CP OR'], surfs['CR IR'], surfs['GT IR'], surfs['GT OR']],
        [mats['AIC'], mats['Air'], mats['SS'], mats['H2O'], mats['Zr'], mats['H2O']],
        grid='intermediate')
    univs['CR nozzle'] = make_pin(
        'CR nozzle',
        [surfs['CP OR'], surfs['CR IR'], surfs['CR OR']],
        [mats['AIC'], mats['Air'], mats['SS'], mats['H2O']])

    univs['CR blank'] = make_pin(
        'CR blank',
        [surfs['CP OR'], surfs['CR IR'], surfs['CR OR'], surfs['GT IR'], surfs['GT OR']],
        [mats['SS'], mats['Air'], mats['SS'], mats['H2O'], mats['Zr'], mats['H2O']])
    univs['CR blank grid (bottom)'] = make_pin(
        'CR blank grid (bottom)',
        [surfs['CP OR'], surfs['CR IR'], surfs['CR OR'], surfs['GT IR'], surfs['GT OR']],
        [mats['SS'], mats['Air'], mats['SS'], mats['H2O'], mats['Zr'], mats['H2O']],
        grid='bottom')
    univs['CR blank grid (intermediate)'] = make_pin(
        'CR blank grid (intermediate)',
        [surfs['CP OR'], surfs['CR IR'], surfs['CR OR'], surfs['GT IR'], surfs['GT OR']],
        [mats['SS'], mats['Air'], mats['SS'], mats['H2O'], mats['Zr'], mats['H2O']],
        grid='intermediate')
    univs['CR blank nozzle'] = make_pin(
        'CR blank nozzle',
        [surfs['CP OR'], surfs['CR IR'], surfs['CR OR']],
        [mats['SS'], mats['Air'], mats['SS'], mats['H2O']])
    univs['CR blank bare'] = make_pin(
        'CR blank bare',
        [surfs['CP OR'], surfs['CR IR'], surfs['CR OR']],
        [mats['SS'], mats['Air'], mats['SS'], mats['H2O']])
    univs['CR bare'] = make_pin(
        'CR bare',
        [surfs['CP OR'], surfs['CR IR'], surfs['CR OR']],
        [mats['AIC'], mats['Air'], mats['SS'], mats['H2O']])

    # Stack all axial pieces of control rod tubes together for each bank

    for b in ['A', 'B', 'C', 'D', 'SA', 'SB', 'SC', 'SD', 'SE']:
        # no grid, no nozzle
        univs['GT CR bank {} dummy'.format(b)] = make_stack(
            'GT CR bank {} dummy'.format(b),
            surfaces=[surfs['bottom FR'],
                      surfs['dashpot top'],
                      surfs['bank{} bot'.format(b)],
                      surfs['bank{} top'.format(b)]],
            universes=[univs['water pin'],
                       univs['GTd empty'],
                       univs['GT empty'],
                       univs['CR'],
                       univs['CR blank']])

        # bottom grid
        univs['GT CR bank {} dummy grid (bottom)'.format(b)] = make_stack(
            'GT CR bank {} dummy grid (bottom)'.format(b),
            surfaces=[surfs['bottom FR'],
                      surfs['dashpot top'],
                      surfs['bank{} bot'.format(b)],
                      surfs['bank{} top'.format(b)]],
            universes=[univs['water pin'],
                       univs['GTd empty grid (bottom)'],
                       univs['GT empty grid (bottom)'],
                       univs['CR grid (bottom)'],
                       univs['CR blank grid (bottom)']])

        # intermediate grid
        univs['GT CR bank {} dummy grid (intermediate)'.format(b)] = make_stack(
            'GT CR bank {} dummy grid (intermediate)'.format(b),
            surfaces=[surfs['bottom FR'],
                      surfs['dashpot top'],
                      surfs['bank{} bot'.format(b)],
                      surfs['bank{} top'.format(b)]],
            universes=[univs['water pin'],
                       univs['GTd empty grid (intermediate)'],
                       univs['GT empty grid (intermediate)'],
                       univs['CR grid (intermediate)'],
                       univs['CR blank grid (intermediate)']])

        # nozzle
        univs['GT CR bank {} dummy nozzle'.format(b)] = make_stack(
            'GT CR bank {} dummy nozzle'.format(b),
            surfaces=[surfs['bottom FR'],
                      surfs['dashpot top'],
                      surfs['bank{} bot'.format(b)],
                      surfs['bank{} top'.format(b)]],
            universes=[univs['water pin'],
                       univs['GTd empty nozzle'],
                       univs['GT empty nozzle'],
                       univs['CR nozzle'],
                       univs['CR blank nozzle']])

        # bare
        univs['GT CR bank {} dummy bare'.format(b)] = make_stack(
            'GT CR bank {} dummy bare'.format(b),
            surfaces=[surfs['bottom FR'],
                      surfs['dashpot top'],
                      surfs['bank{} bot'.format(b)],
                      surfs['bank{} top'.format(b)]],
            universes=[univs['water pin'],
                       univs['GTd empty nozzle'],
                       univs['GT empty nozzle'],
                       univs['CR bare'],
                       univs['CR blank bare']])

        # final combination of all axial pieces for control rod bank "b"
        univs['GT CR bank {}'.format(b)] = make_stack(
            'GT CR bank {}'.format(b), stack_surfs,
            universes=[univs['water pin'],   # number: 21
                       univs['GT CR bank {} dummy nozzle'.format(b)],
                       univs['GT CR bank {} dummy nozzle'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy grid (bottom)'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy grid (intermediate)'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy grid (intermediate)'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy grid (intermediate)'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy grid (bottom)'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy'.format(b)],
                       univs['GT CR bank {} dummy bare'.format(b)],
                       univs['GT CR bank {} dummy bare'.format(b)]])


    #### BURNABLE ABSORBER PIN CELLS
    univs['BA'] = make_pin(
        'BA',
        surfaces=[surfs['BA IR 1'],
                  surfs['BA IR 2'],
                  surfs['BA IR 3'],
                  surfs['BA IR 4'],
                  surfs['BA IR 5'],
                  surfs['BA IR 6'],
                  surfs['BA IR 7'],
                  surfs['BA IR 8']],
        materials=[mats['Air'],
                   mats['SS'],
                   mats['Air'],
                   mats['BSG'],
                   mats['Air'],
                   mats['SS'],
                   mats['H2O'],
                   mats['Zr'],
                   mats['H2O']])

    univs['BA grid (bottom)'] = make_pin(
        'BA grid (bottom)',
        surfaces=[surfs['BA IR 1'],
                  surfs['BA IR 2'],
                  surfs['BA IR 3'],
                  surfs['BA IR 4'],
                  surfs['BA IR 5'],
                  surfs['BA IR 6'],
                  surfs['BA IR 7'],
                  surfs['BA IR 8']],
        materials=[mats['Air'],
                   mats['SS'],
                   mats['Air'],
                   mats['BSG'],
                   mats['Air'],
                   mats['SS'],
                   mats['H2O'],
                   mats['Zr'],
                   mats['H2O']],
        grid='bottom')

    univs['BA grid (intermediate)'] = make_pin(
        'BA grid (intermediate)',
        surfaces=[surfs['BA IR 1'],
                  surfs['BA IR 2'],
                  surfs['BA IR 3'],
                  surfs['BA IR 4'],
                  surfs['BA IR 5'],
                  surfs['BA IR 6'],
                  surfs['BA IR 7'],
                  surfs['BA IR 8']],
        materials=[mats['Air'],
                   mats['SS'],
                   mats['Air'],
                   mats['BSG'],
                   mats['Air'],
                   mats['SS'],
                   mats['H2O'],
                   mats['Zr'],
                   mats['H2O']],
        grid='intermediate')

    univs['BA dashpot'] = make_pin(
        'BA dashpot',
        surfaces=[surfs['BA IR 1'],
                  surfs['BA IR 2'],
                  surfs['BA IR 3'],
                  surfs['BA IR 4'],
                  surfs['BA IR 5'],
                  surfs['BA IR 6'],
                  surfs['GT dashpot IR'],
                  surfs['GT dashpot OR']],
        materials=[mats['Air'],
                   mats['SS'],
                   mats['Air'],
                   mats['BSG'],
                   mats['Air'],
                   mats['SS'],
                   mats['H2O'],
                   mats['Zr'],
                   mats['H2O']])

    univs['BA dashpot grid (bottom)'] = make_pin(
        'BA dashpot grid (bottom)',
        surfaces=[surfs['BA IR 1'],
                  surfs['BA IR 2'],
                  surfs['BA IR 3'],
                  surfs['BA IR 4'],
                  surfs['BA IR 5'],
                  surfs['BA IR 6'],
                  surfs['GT dashpot IR'],
                  surfs['GT dashpot OR']],
        materials=[mats['Air'],
                   mats['SS'],
                   mats['Air'],
                   mats['BSG'],
                   mats['Air'],
                   mats['SS'],
                   mats['H2O'],
                   mats['Zr'],
                   mats['H2O']],
        grid='bottom')

    univs['BA dashpot grid (intermediate)'] = make_pin(
        'BA dashpot grid (intermediate)',
        surfaces=[surfs['BA IR 1'],
                  surfs['BA IR 2'],
                  surfs['BA IR 3'],
                  surfs['BA IR 4'],
                  surfs['BA IR 5'],
                  surfs['BA IR 6'],
                  surfs['GT dashpot IR'],
                  surfs['GT dashpot OR']],
        materials=[mats['Air'],
                   mats['SS'],
                   mats['Air'],
                   mats['BSG'],
                   mats['Air'],
                   mats['SS'],
                   mats['H2O'],
                   mats['Zr'],
                   mats['H2O']],
        grid='intermediate')

    univs['BA blank SS'] = make_pin(
        'BA blank SS',
        surfaces=[surfs['BA IR 6'],
                  surfs['BA IR 7'],
                  surfs['BA IR 8']],
        materials=[mats['SS'],
                   mats['H2O'],
                   mats['Zr'],
                   mats['H2O']])

    univs['BA blank SS bare'] = make_pin(
        'BA blank SS bare',
        surfaces=[surfs['BA IR 6']],
        materials=[mats['SS'],
                   mats['H2O']])

    stack_surfs_BA = [
        surfs['bot support plate'],
        surfs['top support plate'],
        surfs['top lower nozzle'],
        surfs['top lower thimble'],
        surfs['grid1bot'],
        surfs['BA bot'],
        surfs['grid1top'],
        surfs['dashpot top'],
        surfs['grid2bot'],
        surfs['grid2top'],
        surfs['grid3bot'],
        surfs['grid3top'],
        surfs['grid4bot'],
        surfs['grid4top'],
        surfs['top active core'],
        surfs['grid5bot'],
        surfs['grid5top'],
        surfs['top pin plenum'],
        surfs['top FR'],
        surfs['bot upper nozzle'],
        surfs['top upper nozzle']]

    # Stack all axial pieces of control rod tubes together for each bank
    univs['BA stack'] = make_stack(
        'BA stack', stack_surfs_BA,
        universes=[univs['water pin'],
                   univs['water pin'],
                   univs['water pin'],
                   univs['GTd empty'],
                   univs['GTd empty'],
                   univs['GTd empty grid (bottom)'],
                   univs['BA dashpot grid (bottom)'],
                   univs['BA dashpot'],
                   univs['BA'],
                   univs['BA grid (intermediate)'],
                   univs['BA'],
                   univs['BA grid (intermediate)'],
                   univs['BA'],
                   univs['BA grid (intermediate)'],
                   univs['BA'],
                   univs['BA blank SS'],
                   univs['BA blank SS'],
                   univs['BA blank SS'],
                   univs['BA blank SS'],
                   univs['BA blank SS'],
                   univs['BA blank SS bare'],
                   univs['water pin']])


    # Fuel pin cells
    univs['SS pin'] = make_pin(
        'SS pin',
        [surfs['clad OR']],
        [mats['SS'], mats['H2O']])

    univs['end plug'] = make_pin(
        'end plug',
        [surfs['clad OR']],
        [mats['M5'], mats['H2O']])

    univs['pin plenum'] = make_pin(
        'pin plenum',
        surfaces=[surfs['plenum spring OR'],
                  surfs['clad IR'],
                  surfs['clad OR']],
        materials=[mats['SS302'],
                   mats['He'],
                   mats['M5'],
                   mats['H2O']])

    univs['pin plenum grid (intermediate)'] = make_pin(
        'pin plenum grid (intermediate)',
        surfaces=[surfs['plenum spring OR'],
                  surfs['clad IR'],
                  surfs['clad OR']],
        materials=[mats['In'],
                   mats['He'],
                   mats['Zr'],
                   mats['H2O']],
        grid='intermediate')


    #### 1.6% ENRICHED FUEL PIN CELL

    if num_axial > 1:
        # Determine z position between each fuel pellet, omitting the surfaces
        # corresponding to the very bottom and top of the active fuel length
        axial_splits = np.linspace(bottom_fuel_stack, top_active_core, num_axial + 1)[1:-1]
        axial_surfs = [openmc.ZPlane(z0=z) for z in axial_splits]

    if num_rings > 1:
        # Get z-cylinder surfaces for each ring
        rings = []
        for i in range(1, num_rings):
            R = sqrt(i*pellet_OR**2/num_rings)
            cyl = openmc.ZCylinder(r=R, name='fuel ring {}'.format(i))
            rings.append(cyl)

    def subdivided_fuel(fill):
        # Create universe for UO2 alone with axial/radial subdivision
        uo2_cells = []
        if num_axial > 1:
            for axial_region in subdivide(axial_surfs):
                if num_rings > 1:
                    for ring_region in subdivide(rings):
                        cell = openmc.Cell(fill=fill, region=axial_region & ring_region)
                        uo2_cells.append(cell)
                else:
                    uo2_cells.append(openmc.Cell(fill=fill, region=axial_region))
        else:
            if num_rings > 1:
                for ring_region in subdivide(rings):
                    cell = openmc.Cell(fill=fill, region=ring_region)
                    uo2_cells.append(cell)
            else:
                raise RuntimeError("Shouldn't call with 1 ring and 1 axial segment")

        return openmc.Universe(cells=uo2_cells)

    # If rings/axial segments are present, create a universe for the subdivided
    # fuel. Otherwise just use a plain material.
    if num_rings > 1 or num_axial > 1:
        fuel_fill = subdivided_fuel(mats['UO2 1.6 {}'.format(fuel)])
    else:
        fuel_fill = mats['UO2 1.6 {}'.format(fuel)]

    outside_pin_surfaces = [surfs['clad IR'], surfs['clad OR']]
    outside_pin_mats = [mats['He'], mats['M5'], mats['H2O']]

    univs['Outside pin'] = make_pin(
        'Outside pin',
        surfaces=outside_pin_surfaces,
        materials=outside_pin_mats)

    univs['Outside pin grid (bottom)'] = make_pin(
        'Outside pin grid (bottom)',
        surfaces=outside_pin_surfaces,
        materials=outside_pin_mats,
        grid='bottom')

    univs['Outside pin grid (intermediate)'] = make_pin(
        'Outside pin grid (intermediate)',
        surfaces=outside_pin_surfaces,
        materials=outside_pin_mats,
        grid='intermediate')

    # Stack all axial pieces of 1.6% enriched fuel pin cell

    within_fuel_surfs = [
        surfs['grid1bot'],
        surfs['grid1top'],
        surfs['dashpot top'],
        surfs['grid2bot'],
        surfs['grid2top'],
        surfs['grid3bot'],
        surfs['grid3top'],
        surfs['grid4bot'],
        surfs['grid4top']
    ]

    univs['Fuel pin (1.6%) stack'] = make_pin_stack(
        'Fuel pin (1.6%) stack',
        zsurfaces=within_fuel_surfs,
        universes=[
            univs['Outside pin'],
            univs['Outside pin grid (bottom)'],
            univs['Outside pin'],
            univs['Outside pin'],
            univs['Outside pin grid (intermediate)'],
            univs['Outside pin'],
            univs['Outside pin grid (intermediate)'],
            univs['Outside pin'],
            univs['Outside pin grid (intermediate)'],
            univs['Outside pin']
        ],
        boundary=surfs['pellet OR'],
        fuel_fill=fuel_fill)

    fuel_stack_surfs = [
        surfs['bot support plate'],
        surfs['top support plate'],
        surfs['top lower nozzle'],
        surfs['top lower thimble'],
        surfs['top active core'],
        surfs['grid5bot'],
        surfs['grid5top'],
        surfs['top pin plenum'],
        surfs['top FR'],
        surfs['bot upper nozzle'],
        surfs['top upper nozzle']
    ]

    univs['Fuel (1.6%) stack'] = make_stack(
        'Fuel (1.6%) stack',
        surfaces=fuel_stack_surfs,
        universes=[univs['water pin'],
                   univs['SS pin'],
                   univs['SS pin'],
                   univs['end plug'],
                   univs['Fuel pin (1.6%) stack'],
                   univs['pin plenum'],
                   univs['pin plenum grid (intermediate)'],
                   univs['pin plenum'],
                   univs['end plug'],
                   univs['water pin'],
                   univs['SS pin'],
                   univs['water pin']])


    #### 2.4% ENRICHED FUEL PIN CELL

    # If rings/axial segments are present, create a universe for the subdivided
    # fuel. Otherwise just use a plain material.
    if num_rings > 1 or num_axial > 1:
        fuel_fill = subdivided_fuel(mats['UO2 2.4 {}'.format(fuel)])
    else:
        fuel_fill = mats['UO2 2.4 {}'.format(fuel)]

    # Stack all axial pieces of 2.4% enriched fuel pin cell

    univs['Fuel pin (2.4%) stack'] = make_pin_stack(
        'Fuel pin (2.4%) stack',
        zsurfaces=within_fuel_surfs,
        universes=[
            univs['Outside pin'],
            univs['Outside pin grid (bottom)'],
            univs['Outside pin'],
            univs['Outside pin'],
            univs['Outside pin grid (intermediate)'],
            univs['Outside pin'],
            univs['Outside pin grid (intermediate)'],
            univs['Outside pin'],
            univs['Outside pin grid (intermediate)'],
            univs['Outside pin']
        ],
        boundary=surfs['pellet OR'],
        fuel_fill=fuel_fill)

    univs['Fuel (2.4%) stack'] = make_stack(
        'Fuel (2.4%) stack',
        surfaces=fuel_stack_surfs,
        universes=[univs['water pin'],
                   univs['SS pin'],
                   univs['SS pin'],
                   univs['end plug'],
                   univs['Fuel pin (2.4%) stack'],
                   univs['pin plenum'],
                   univs['pin plenum grid (intermediate)'],
                   univs['pin plenum'],
                   univs['end plug'],
                   univs['water pin'],
                   univs['SS pin'],
                   univs['water pin']])


    #### 3.1% ENRICHED FUEL PIN CELL

    # If rings/axial segments are present, create a universe for the subdivided
    # fuel. Otherwise just use a plain material.
    if num_rings > 1 or num_axial > 1:
        fuel_fill = subdivided_fuel(mats['UO2 3.1 {}'.format(fuel)])
    else:
        fuel_fill = mats['UO2 3.1 {}'.format(fuel)]

    # Stack all axial pieces of 3.1% enriched fuel pin cell

    univs['Fuel pin (3.1%) stack'] = make_pin_stack(
        'Fuel pin (3.1%) stack',
        zsurfaces=within_fuel_surfs,
        universes=[
            univs['Outside pin'],
            univs['Outside pin grid (bottom)'],
            univs['Outside pin'],
            univs['Outside pin'],
            univs['Outside pin grid (intermediate)'],
            univs['Outside pin'],
            univs['Outside pin grid (intermediate)'],
            univs['Outside pin'],
            univs['Outside pin grid (intermediate)'],
            univs['Outside pin']
        ],
        boundary=surfs['pellet OR'],
        fuel_fill=fuel_fill)

    univs['Fuel (3.1%) stack'] = make_stack(
        'Fuel (3.1%) stack',
        surfaces=fuel_stack_surfs,
        universes=[univs['water pin'],
                   univs['SS pin'],
                   univs['SS pin'],
                   univs['end plug'],
                   univs['Fuel pin (3.1%) stack'],
                   univs['pin plenum'],
                   univs['pin plenum grid (intermediate)'],
                   univs['pin plenum'],
                   univs['end plug'],
                   univs['water pin'],
                   univs['SS pin'],
                   univs['water pin']])

    return univs
