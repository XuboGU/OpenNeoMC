"""Stainless steel heavy reflector.

None of the public NuScale documents give information about the dimensions and
location of the water holes in the heavy neutron reflector. Thus, all the
dimensions of the water holes in the heavy reflectors were eyeballed from the
design certification application (ML17013A274), Figure 4.3-25. A screenshot was
used to determine the size and location of the water holes which were then
converted to actual dimensions by scaling according to the width of an assembly.

"""

import openmc

from .materials import mats
from .surfaces import lattice_pitch


def make_reflector(name, parameters):
    """Make an assembly-sized heavy neutron reflector block with cooling holes.

    Parameters
    ----------
    name : str
        Name of the universe to create
    parameters : iterable of 3-tuples
        Iterable containing tuple with the (x,y) coordinates of the center and
        the radius of a Z-cylinder and the

    Returns
    -------
    openmc.Universe
        Universe containing reflector block

    """
    water_holes = []
    for x, y, r in parameters:
        zcyl = openmc.ZCylinder(x0=x, y0=y, r=r)
        hole = openmc.Cell(fill=mats['H2O'], region=-zcyl)
        water_holes.append(hole)

    ss_region = openmc.Intersection(~c.region for c in water_holes)
    ss_cell = openmc.Cell(name='reflector {} SS'.format(name), fill=mats['SS'],
                          region=ss_region)

    univ = openmc.Universe(name='reflector {}'.format(name))
    univ.add_cells(water_holes)
    univ.add_cell(ss_cell)
    return univ


def reflector_universes():
    """Generate universes for SMR heavy neutron reflector blocks.

    Parameters
    ----------
    num_rings : int
        Number of annual regions in fuel
    num_axial : int
        Number of axial subdivisions in fuel

    Returns
    -------
    dict
        Dictionary mapping a universe name to a openmc.Universe object

    """
    # Create dictionary to store universes
    univs = {}

    # Reflector at northwest corner (fuel assemblies to the right and below)
    width = 276
    p1 = 59
    p2 = 126
    p3 = 196
    p4 = 264

    p5 = 105

    p6 = 122
    p7 = 164

    p8 = 138
    p9 = 222

    p10 = 247

    # There are 8 large water holes and all others appear to have the same, smaller
    # diameter
    d_small = 13
    d_large = 30

    # All pixel widths are scaled according to the actual width of an assembly
    # divided by the width of an assembly in pixels
    scale = lattice_pitch/width

    # Physical positions
    x1 = -lattice_pitch/2 + scale*(width - p4)
    x2 = -lattice_pitch/2 + scale*(width - p3)
    x3 = -lattice_pitch/2 + scale*(width - p2)
    x4 = -lattice_pitch/2 + scale*(width - p1)
    y1 = -lattice_pitch/2 + scale*p1
    y2 = -lattice_pitch/2 + scale*p2
    y3 = -lattice_pitch/2 + scale*p3
    y4 = -lattice_pitch/2 + scale*p4

    x5 = -lattice_pitch/2 + scale*(width - p5)
    y5 = -lattice_pitch/2 + scale*p5
    x6 = -lattice_pitch/2 + scale*(width - p7)
    y6 = -lattice_pitch/2 + scale*p6
    x7 = -lattice_pitch/2 + scale*(width - p6)
    y7 = -lattice_pitch/2 + scale*p7
    x8 = -lattice_pitch/2 + scale*(width - p9)
    y8 = -lattice_pitch/2 + scale*p8
    x9 = -lattice_pitch/2 + scale*(width - p8)
    y9 = -lattice_pitch/2 + scale*p9

    y10 = -lattice_pitch/2 + scale*p10

    # Radius of small/large water holes
    r1 = scale*d_small/2
    r2 = scale*d_large/2

    params = [
        (x1, y1, r1), (x2, y1, r1), (x3, y1, r1), (x4, y1, r2),
        (x4, y2, r1), (x4, y3, r1), (x4, y4, r1), (x5, y5, r1),
        (x6, y6, r1), (x7, y7, r1), (x8, y8, r1), (x9, y9, r1),
        (x1, y10, r1)
    ]
    univs['NW'] = make_reflector('NW', params)

    # Reflector at (1, 1)

    params = [
        (x4, y1, r1),
        (lattice_pitch/2 - scale*103, -lattice_pitch/2 + scale*156, r1),
        (lattice_pitch/2 - scale*158, -lattice_pitch/2 + scale*103, r1)
    ]
    univs['1,1'] = make_reflector('1,1', params)

    # Left reflector (4,0)

    left1 = 58
    left2 = 118
    left3 = 173
    up3 = 76

    x1 = -lattice_pitch/2 + scale*(width - left1)
    x2 = -lattice_pitch/2 + scale*(width - left2)
    d_y = scale*67
    x3 = -lattice_pitch/2 + scale*(width - left3)
    y3 = scale*up3

    params = [
        (x1, 0, r1), (x1, d_y, r1), (x1, 2*d_y, r1), (x1, -d_y, r1), (x1, -2*d_y, r1),
        (x2, d_y/2, r1), (x2, 3/2*d_y, r1), (x2, -d_y/2, r1), (x2, -3/2*d_y, r1),
        (x3, y3, r1), (x3, -y3, r1)
       ]
    univs['4,0'] = make_reflector('4,0', params)

    # Reflector at (3,0)

    params = []
    for i in range(2, 7):
        params.append((x1, i*d_y - lattice_pitch, r1))
    for i in (5, 7, 11):
        params.append((x2, i*d_y/2 - lattice_pitch, r1))

    left3 = 140
    left4 = 183
    up3 = 159
    up4 = 47

    x3 = -lattice_pitch/2 + scale*(width - left3)
    y3 = -lattice_pitch/2 + scale*up3
    x4 = -lattice_pitch/2 + scale*(width - left4)
    y4 = -lattice_pitch/2 + scale*up4
    params += [(x3, y3, r1), (x4, y4, r1)]

    univs['3,0'] = make_reflector('3,0', params)

    # Reflector at (5,0)
    params = [(x, -y, r) for x, y, r in params]
    univs['5,0'] = make_reflector('5,0', params)

    # Reflector at (2, 0)

    params = [(-lattice_pitch/2 + scale*(width - 78),
               -lattice_pitch/2 + scale*98, r1)]
    univs['2,0'] = make_reflector('2,0', params)

    ################################################################################
    # Beyond this point, all universes are just copies of the ones previously
    # created with a rotation applied

    # First define helper function to create new universe by rotating an
    # existing one
    def rotate_universe(univ, rotation, name):
        cell = openmc.Cell(name='reflector {}'.format(name), fill=univ)
        cell.rotation = rotation
        return openmc.Universe(name=name, cells=[cell])

    univs['NE'] = rotate_universe(univs['NW'], (0, 0, -90), 'NE')
    univs['SW'] = rotate_universe(univs['NW'], (0, 0, 90), 'SW')
    univs['SE'] = rotate_universe(univs['NW'], (0, 0, 180), 'SE')
    univs['0,2'] = rotate_universe(univs['2,0'], (0, 180, -90), '0,2')
    univs['0,3'] = rotate_universe(univs['5,0'], (0, 0, -90), '0,3')
    univs['0,4'] = rotate_universe(univs['4,0'], (0, 0, -90), '0,4')
    univs['0,5'] = rotate_universe(univs['3,0'], (0, 0, -90), '0,5')
    univs['0,6'] = rotate_universe(univs['2,0'], (0, 0, -90), '0,6')
    univs['1,7'] = rotate_universe(univs['1,1'], (0, 0, -90), '1,7')
    univs['2,8'] = rotate_universe(univs['2,0'], (0, 180, 0), '2,8')
    univs['3,8'] = rotate_universe(univs['3,0'], (0, 180, 0), '3,8')
    univs['4,8'] = rotate_universe(univs['4,0'], (0, 180, 0), '4,8')
    univs['5,8'] = rotate_universe(univs['3,0'], (0, 0, 180), '5,8')
    univs['6,0'] = rotate_universe(univs['2,0'], (180, 0, 0), '6,0')
    univs['6,8'] = rotate_universe(univs['2,0'], (0, 0, 180), '6,8')
    univs['7,1'] = rotate_universe(univs['1,1'], (180, 0, 0), '7,1')
    univs['7,7'] = rotate_universe(univs['1,1'], (0, 0, 180), '7,7')
    univs['8,2'] = rotate_universe(univs['2,0'], (0, 0, 90), '8,2')
    univs['8,3'] = rotate_universe(univs['3,0'], (0, 0, 90), '8,3')
    univs['8,4'] = rotate_universe(univs['4,0'], (0, 0, 90), '8,4')
    univs['8,5'] = rotate_universe(univs['5,0'], (0, 0, 90), '8,5')
    univs['8,6'] = rotate_universe(univs['2,0'], (0, 0, 180), '8,6')

    # Solid stainless steel universe
    all_ss = openmc.Cell(name='heavy reflector', fill=mats['SS'])
    univs['solid'] = openmc.Universe(name='solid', cells=[all_ss])

    return univs
