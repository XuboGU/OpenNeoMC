"""Instantiate OpenMC Plots to visualize the core model."""

from random import randint, seed

import openmc

from .surfaces import lowest_extent, highest_extent, lattice_pitch, rpv_OR, \
    pin_pitch
from .materials import mats


# color specifications
colors = {
    mats['H2O']: (198, 226, 255),  # light blue
    mats['In']: (101, 101, 101),  # dgray
    mats['CS']: (  0,   0,   0),  # carbons black
    mats['Zr']: (201, 201, 201),  # gray
    mats['SS']: (  0,   0,   0),  # black
    mats['Air']: (255, 255, 255),  # white
    mats['He']: (255, 218, 185),  # light orange
    mats['BSG']: (  0, 255,   0),  # green
    mats['AIC']: (255,   0,   0),  # bright red
    mats['UO2 1.6 fresh']: (142,  35,  35),  # light red
    mats['UO2 2.4 fresh']: (255, 215,   0),  # gold
    mats['UO2 3.1 fresh']: (  0,   0, 128)   # dark blue
}


def core_plots():
    # Create a collection of plots
    plots = openmc.Plots()

    plot = openmc.Plot(name='radial slice ')
    plot.basis = 'xy'
    plot.color_by = 'material'
    plot.origin = [0., 0., (highest_extent-lowest_extent)/2.]
    plot.width = [25*lattice_pitch/2, 25*lattice_pitch/2.]
    plot.filename = 'radial_xy_slice'
    plot.colors = colors
    plot.background = [255, 255, 255]
    plot.pixels = [1000, 1000]
    plots.append(plot)

    plot = openmc.Plot(name='axial slice')
    plot.basis = 'xz'
    plot.color_by = 'material'
    plot.origin = [0., 0., (highest_extent-lowest_extent)/2.]
    plot.width = [rpv_OR*2., (highest_extent-lowest_extent)]
    plot.filename = 'axial_xz_slice'
    plot.colors = colors
    plot.background = [255, 255, 255]
    plot.pixels = [1000, 1000]
    plots.append(plot)

    plot = openmc.Plot(name='assembly grid spacer')
    plot.basis = 'xy'
    plot.color_by = 'material'
    plot.origin = [0., 0., 95.0]
    plot.width = [lattice_pitch*1.5, lattice_pitch*1.5]
    plot.filename = 'assm_grid_spacer'
    plot.colors = colors
    plot.background = [255, 255, 255]
    plot.pixels = [2000, 2000]
    plots.append(plot)

    plot = openmc.Plot(name='assembly no spacer')
    plot.basis = 'xy'
    plot.color_by = 'material'
    plot.origin = [0., 0., 90.]
    plot.width = [lattice_pitch*1.5, lattice_pitch*1.5]
    plot.filename = 'assm_no_spacer'
    plot.colors = colors
    plot.background = [255, 255, 255]
    plot.pixels = [2000, 2000]
    plots.append(plot)

    plot = openmc.Plot(name='assembly no spacer cell')
    plot.basis = 'xy'
    plot.color_by = 'cell'
    plot.origin = [0., 0., 90.]
    plot.width = [lattice_pitch*1.5, lattice_pitch*1.5]
    plot.filename = 'assm_no_spacer_cell'
    plot.background = [255, 255, 255]
    plot.pixels = [2000, 2000]
    plots.append(plot)

    plot = openmc.Plot(name='z slice')
    plot.basis = 'xz'
    plot.color_by = 'cell'
    plot.filename = 'assm_xz'
    plot.origin = (0., 0., lowest_extent + (highest_extent - lowest_extent)/2)
    plot.width = (lattice_pitch*1.5, highest_extent - lowest_extent)
    plot.pixels = (int(plot.width[0]/plot.width[1]*4000), 4000)
    plots.append(plot)

    return plots


def assembly_plots(univ):
    colors = {
        mats['H2O']: (198, 226, 255),  # light blue
        mats['In']: (101, 101, 101),  # dgray
        mats['Zr']: (201, 201, 201),  # gray
        mats['SS']: (  0,   0,   0),  # black
        mats['Air']: (255, 255, 255),  # white
        mats['He']: (255, 218, 185),  # light orange
        mats['BSG']: (  0, 255,   0),  # green
    }

    # Set colors for fuel as shades of yellow
    seed(1)
    for m in univ.get_all_materials().values():
        if 'UO2 Fuel' in m.name:
            colors[m] = (randint(200, 255), randint(200, 255), randint(0, 150))

    pin = openmc.Plot(name='fuel pin')
    pin.basis = 'xy'
    pin.color_by = 'material'
    pin.colors = colors
    pin.filename = 'smr_pin'
    pin.origin = (pin_pitch, 0., (highest_extent-lowest_extent)/2.)
    pin.width = (pin_pitch, pin_pitch)
    pin.pixels = (400, 400)

    assm = openmc.Plot(name='fuel assembly')
    assm.basis = 'xy'
    assm.color_by = 'material'
    assm.colors = colors
    assm.filename = 'smr_assm'
    assm.origin = (0., 0., (highest_extent-lowest_extent)/2.)
    assm.width = (lattice_pitch, lattice_pitch)
    assm.pixels = (2000, 2000)

    return openmc.Plots([pin, assm])
