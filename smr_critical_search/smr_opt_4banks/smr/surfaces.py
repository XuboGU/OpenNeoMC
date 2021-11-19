"""Instantiate the OpenMC Surfaces needed by the core model.

The geometric parameters defining the core model are tabulated here.
The geometric specifications are loosely based upon NuScale's Small
Modular Pressurized Water Reactor concept as detailed here in their
NRC design certification (DC) documentation:

https://www.nrc.gov/docs/ML1618/ML16187A017.pdf
https://www.nrc.gov/docs/ML1700/ML17007A001.pdf

NuScale DC application, chapter 1: https://www.nrc.gov/docs/ML1701/ML17013A264.pdf
NuScale DC application, chapter 4: https://www.nrc.gov/docs/ML1701/ML17013A274.pdf

"""

import copy
from math import tan, pi

import numpy as np
import openmc

INCHES = 2.54

# Notation
# GT: Guide Tube
# BA: Burnable Absorber
# CP: Control Poison
# FR: Fuel Rod
# IR: Inner Radius
# OR: Outer Radius
# IT: Instrument Tube
# FA: Fuel Assembly
# RPV: Reactor Pressure Vessel

# fuel rod parameters
pellet_OR          = 0.3195*INCHES/2  # ML17013A274, Table 4.1-2
pellet_length      = 0.4*INCHES       # ML17013A274, Table 4.1-2
clad_IR            = 0.326*INCHES/2   # ML17013A274, Table 4.1-2
clad_OR            = 0.374*INCHES/2   # ML17013A274, Table 4.1-2
active_fuel_length = 78.74*INCHES     # ML17013A274, Figure 4.2-10
plenum_length      = 5.311*INCHES     # ML17013A274, Figure 4.2-10
fuel_rod_length    = 85.00*INCHES     # ML17013A274, Table 4.1-2
lower_end_cap_length = 0.575*INCHES   # ML17007A001, Table 3-2

# pin cell parameters
guide_tube_IR      = 0.450*INCHES/2  # ML17013A274, Table 4.1-2
guide_tube_OR      = 0.482*INCHES/2  # ML17013A274, Table 4.1-2
guide_tube_dash_IR = 0.397*INCHES/2  # ML17013A274, Table 4.1-2
guide_tube_dash_OR = guide_tube_OR
boron_carbide_OR   = 0.333*INCHES/2  # ML17013A274, Table 4.1-3
ag_in_cd_OR        = 0.336*INCHES/2  # ML17013A274, Table 4.1-3
control_rod_IR     = 0.344*INCHES/2  # ML17013A274, Table 4.1-3
control_rod_OR     = 0.381*INCHES/2  # ML17013A274, Table 4.1-3
burn_abs_r1        = 0.21400
burn_abs_r2        = 0.23051
burn_abs_r3        = 0.24130
burn_abs_r4        = 0.42672
burn_abs_r5        = 0.43688
burn_abs_r6        = 0.48387
burn_abs_r7        = 0.56134
burn_abs_r8        = 0.60198
instr_tube_IR      = 0.450*INCHES/2  # ML17007A001, Table 3-1
instr_tube_OR      = 0.482*INCHES/2  # ML17007A001, Table 3-1
plenum_spring_OR   = 0.06459  # Estimate, actual is ECI

# grid spacer parameters
rod_grid_side = 1.24416
spacer_height = 1.750*INCHES  # ML17013A274, Figure 4.2-7

# assembly parameters
assembly_length   = 95.89*INCHES  # ML17013A274, Table 4.1-2
pin_pitch         = 0.496*INCHES  # ML17013A274, Table 4.1-2
lattice_pitch     = 8.466*INCHES  # ML17013A274, Table 4.1-2
grid_strap_side   = 21.47270
top_nozzle_height = 3.551*INCHES  # ML17013A274, Figure 4.2-2
top_nozzle_width  = 8.406*INCHES  # ML17013A274, Figure 4.2-2

# core radial parameters
core_barrel_IR     = 74*INCHES/2  # ML17013A274, Table 4.1-2
core_barrel_OR     = 78*INCHES/2  # ML17013A274, Table 4.1-2
neutron_shield_OR  = core_barrel_OR + 2.0
rpv_IR             = 96.5*INCHES/2  # ML17013A274, Table 5.3-1
rpv_OR             = 105*INCHES/2   # ML17013A274, Table 5.3-1

# axial parameters
reference_z = 0.0
lowest_extent        = reference_z
bottom_support_plate = lowest_extent + 20.000
top_support_plate    = bottom_support_plate + 5.000
bottom_lower_nozzle  = bottom_support_plate + 5.000
top_lower_nozzle     = bottom_lower_nozzle + 4.0*INCHES  # nozzle 喷嘴
bottom_fuel_rod      = bottom_lower_nozzle + 4.0*INCHES  # 35.16
top_lower_thimble    = bottom_fuel_rod + lower_end_cap_length # thimble 顶针
bottom_fuel_stack    = bottom_fuel_rod + lower_end_cap_length
bot_burn_abs         = bottom_fuel_stack + 2.0*INCHES
top_active_core = bottom_fuel_stack + active_fuel_length
top_plenum = top_active_core + plenum_length
top_fuel_rod = bottom_fuel_rod + fuel_rod_length
bottom_upper_nozzle = top_fuel_rod + (423.049 - 419.704)     # BEAVRS, Fig. 32
top_upper_nozzle = bottom_upper_nozzle + (431.876 - 423.049) # BEAVRS, Fig. 32
highest_extent = top_upper_nozzle + 20.0

# The grid spacer locations are eyeball estimated from Figure 3-1 in NuScale's
# FA design certification doc, ML17007A001. This assumes 6cm and 2cm spacings
# between the bottom and top of the fuel rods and the bottom and top grid
# spacers.
first_grid_bot = bottom_fuel_rod + 6.0 # 41.16
last_grid_top = top_fuel_rod - 2.0
last_grid_bot = last_grid_top - spacer_height

grid_bottom = np.linspace(first_grid_bot, last_grid_bot, 5)
grid_top = grid_bottom + spacer_height

# control rod step heights - taken from BEAVRS, use with caution for NuScale
step0H                =    46.079   # temporary value for now
step102H              =   206.415
step248H              =   269.122
step_width            =     1.58173
bank_bot              =   405.713
bank_step             =   248.
bank_top              =   786.348

neutron_shield_NWbot_SEtop = tan(pi/3)
neutron_shield_NWtop_SEbot = tan(pi/6)
neutron_shield_NEbot_SWtop = tan(-pi/3)
neutron_shield_NEtop_SWbot = tan(-pi/6)


surfs = {}

surfs['pellet OR'] = openmc.ZCylinder(
    r=pellet_OR, name='Pellet OR')
surfs['plenum spring OR'] = openmc.ZCylinder(
    r=plenum_spring_OR, name='FR Plenum Spring OR')
surfs['clad IR'] = openmc.ZCylinder(
    r=clad_IR, name='Clad IR')
surfs['clad OR'] = openmc.ZCylinder(
    r=clad_OR, name='Clad OR')
surfs['GT IR'] = openmc.ZCylinder(
    r=guide_tube_IR, name='GT IR (above dashpot)')
surfs['GT OR'] = openmc.ZCylinder(
    r=guide_tube_OR, name='GT OR (above dashpot)')
surfs['GT dashpot IR'] = openmc.ZCylinder(
    r=guide_tube_dash_IR, name='GT IR (at dashpot)')
surfs['GT dashpot OR'] = openmc.ZCylinder(
    r=guide_tube_dash_OR, name='GT OR (at dashpot)')
surfs['CP OR'] = openmc.ZCylinder(
    r=boron_carbide_OR, name='Control Poison OR')
surfs['CR IR'] = openmc.ZCylinder(
    r=control_rod_IR, name='CR Clad IR')
surfs['CR OR'] = openmc.ZCylinder(
    r=control_rod_OR, name='CR Clad OR')
surfs['BA IR 1'] = openmc.ZCylinder(
    r=burn_abs_r1, name='BA IR 1')
surfs['BA IR 2'] = openmc.ZCylinder(
    r=burn_abs_r2, name='BA IR 2')
surfs['BA IR 3'] = openmc.ZCylinder(
    r=burn_abs_r3, name='BA IR 3')
surfs['BA IR 4'] = openmc.ZCylinder(
    r=burn_abs_r4, name='BA IR 4')
surfs['BA IR 5'] = openmc.ZCylinder(
    r=burn_abs_r5, name='BA IR 5')
surfs['BA IR 6'] = openmc.ZCylinder(
    r=burn_abs_r6, name='BA IR 6')
surfs['BA IR 7'] = openmc.ZCylinder(
    r=burn_abs_r7, name='BA IR 7')
surfs['BA IR 8'] = openmc.ZCylinder(
    r=burn_abs_r8, name='BA IR 8')
# surfs['IT IR'] = copy.deepcopy(surfs['BA IR 5'])
# surfs['IT OR'] = copy.deepcopy(surfs['BA IR 6'])
surfs['IT IR'] = surfs['BA IR 5'].clone()
surfs['IT OR'] = surfs['BA IR 6'].clone()

# Rectangular prisms for grid spacers
surfs['rod grid box'] = \
    openmc.rectangular_prism(rod_grid_side, rod_grid_side)

# Rectangular prisms for lattice grid sleeves
surfs['lat grid box inner'] = \
    openmc.rectangular_prism(17.*pin_pitch, 17.*pin_pitch)
surfs['lat grid box outer'] = \
    openmc.rectangular_prism(grid_strap_side, grid_strap_side)

surfs['bot support plate'] = openmc.ZPlane(
    z0=bottom_support_plate, name='bot support plate')
surfs['top support plate'] = openmc.ZPlane(
    z0=top_support_plate, name='top support plate')
surfs['bottom FR'] = openmc.ZPlane(z0=bottom_fuel_rod, name='bottom FR')
# surfs['top lower nozzle'] = copy.deepcopy(surfs['bottom FR'])
# surfs['bot lower nozzle'] = copy.deepcopy(surfs['top support plate'])
surfs['top lower nozzle'] = surfs['bottom FR'].clone()
surfs['bot lower nozzle'] = surfs['top support plate'].clone()

# axial surfaces
surfs['bot active core'] = openmc.ZPlane(
    z0=bottom_fuel_stack, name='bot active core')
surfs['top active core'] = openmc.ZPlane(
    z0=top_active_core, name='top active core')

# surfs['top lower thimble'] = copy.deepcopy(surfs['bot active core'])
surfs['top lower thimble'] = surfs['bot active core'].clone()

surfs['BA bot'] = openmc.ZPlane(
    z0=bot_burn_abs, name='bottom of BA')

for i, (bottom, top) in enumerate(zip(grid_bottom, grid_top)):
    # Create plane for bottom of spacer grid
    key = 'grid{}bot'.format(i + 1)
    name = 'bottom grid {}'.format(i + 1)
    surfs[key] = openmc.ZPlane(z0=bottom, name=name)

    # Create plane for top of spacer grid
    key = 'grid{}top'.format(i + 1)
    name = 'top of grid {}'.format(i + 1)
    surfs[key] = openmc.ZPlane(z0=top, name=name)

surfs['dashpot top'] = openmc.ZPlane(
    z0=step0H, name='top dashpot')

surfs['top pin plenum'] = openmc.ZPlane(
    z0=top_plenum, name='top pin plenum')
surfs['top FR'] = openmc.ZPlane(
    z0=top_fuel_rod, name='top FR')
surfs['bot upper nozzle'] = openmc.ZPlane(
    z0=bottom_upper_nozzle, name='bottom upper nozzle')
surfs['top upper nozzle'] = openmc.ZPlane(
    z0=top_upper_nozzle, name='top upper nozzle')

# Control rod bank surfaces for ARO configuration
for bank in ['A','B','C','D','E',]:
    surfs['bankS{} top'.format(bank)] = openmc.ZPlane(
        z0=step248H+step_width*228, name='CR bankS{} top'.format(bank))
    surfs['bankS{} bot'.format(bank)] = openmc.ZPlane(
        z0=step248H, name='CR bankS{} bottom'.format(bank))

# surfs['bankA top'] = openmc.ZPlane(
#     z0=bank_top, name='CR bank A top')
# surfs['bankA bot'] = openmc.ZPlane(
#     z0=bank_bot, name='CR bank A bottom')
# surfs['bankB top'] = openmc.ZPlane(
#     z0=bank_top, name='CR bank B top')
# surfs['bankB bot'] = openmc.ZPlane(
#     z0=bank_bot, name='CR bank B bottom')
# surfs['bankC top'] = openmc.ZPlane(
#     z0=bank_top, name='CR bank C top')
# surfs['bankC bot'] = openmc.ZPlane(
#     z0=bank_bot, name='CR bank C bottom')
# surfs['bankD top'] = openmc.ZPlane(s
#     z0=bank_top, name='CR bank D top')
# surfs['bankD bot'] = openmc.ZPlane(
#     z0=bank_bot, name='CR bank D bottom')

# outer radial surfaces
surfs['core barrel IR'] = openmc.ZCylinder(
    r=core_barrel_IR, name='core barrel IR')
surfs['core barrel OR'] = openmc.ZCylinder(
    r=core_barrel_OR, name='core barrel OR')
surfs['neutron shield OR'] = openmc.ZCylinder(
    r=neutron_shield_OR, name='neutron shield OR')

# neutron shield planes
surfs['neutron shield NWbot SEtop'] = openmc.Plane(
    a=1., b=neutron_shield_NWbot_SEtop, c=0., d=0.,
    name='neutron shield NWbot SEtop')
surfs['neutron shield NWtop SEbot'] = openmc.Plane(
    a=1., b=neutron_shield_NWtop_SEbot, c=0., d=0.,
    name='neutron shield NWtop SEbot')
surfs['neutron shield NEbot SWtop'] = openmc.Plane(
    a=1., b=neutron_shield_NEbot_SWtop, c=0., d=0.,
    name='neutron shield NEbot SWtop')
surfs['neutron shield NEtop SWbot'] = openmc.Plane(
    a=1., b=neutron_shield_NEtop_SWbot, c=0., d=0.,
    name='neutron shield NEtop SWbot')

# outer radial surfaces
surfs['RPV IR'] = openmc.ZCylinder(
    r=rpv_IR, name='RPV IR')
surfs['RPV OR'] = openmc.ZCylinder(
    r=rpv_OR, name='RPV OR', boundary_type='vacuum')

# outer axial surfaces
surfs['upper bound'] = openmc.ZPlane(
    z0=highest_extent, name='upper problem boundary', boundary_type='vacuum')
surfs['lower bound'] = openmc.ZPlane(
    z0=lowest_extent, name='lower problem boundary', boundary_type='vacuum')
