"""Instantiate the OpenMC Materials needed by the core model."""

import openmc
from openmc.data import atomic_weight, atomic_mass, water_density

from . import system_pressure, core_average_temperature


_DEPLETION_NUCLIDES = [
    "U239", "U240", "Np234", "Np235", "Np236", "Np237", "Np238", "Np239",
    "Pu236", "Pu237", "Pu238", "Pu239", "Pu240", "Pu241", "Pu242", "B11",
    "N14", "N15", "Fe57", "Fe58", "Co59", "Ni60", "Ni61", "Ni62",
    "Cu63", "Ni64", "Zn64", "Cu65", "Zn65", "Zn66", "Zn67", "Zn68",
    "Ga69", "Zn70", "Ge70", "Ga71", "Ge72", "Ge73", "Ge74", "As74",
    "Se74", "As75", "Ge76", "Se76", "Se77", "Se78", "Se79", "Br79",
    "Se80", "Kr80", "Br81", "Se82", "Kr82", "Kr83", "Kr84", "Sr84",
    "Kr85", "Rb85", "Kr86", "Rb86", "Sr86", "Rb87", "Sr87", "Sr88",
    "Sr89", "Y89", "Sr90", "Y90", "Zr90", "Y91", "Zr91", "Zr92",
    "Zr93", "Nb93", "Zr94", "Nb94", "Mo94", "Zr95", "Nb95", "Mo95",
    "Zr96", "Mo96", "Mo97", "Mo98", "Ru98", "Mo99", "Tc99", "Ru99",
    "Mo100", "Ru100", "Ru101", "Ru102", "Pd102", "Ru103", "Rh103", "Ru104",
    "Pd104", "Ru105", "Rh105", "Pd105", "Ru106", "Pd106", "Pd107", "Ag107",
    "Pd108", "Cd108", "Ag109", "Pd110", "Cd110", "Ag111", "Cd111", "Cd112",
    "Sn112", "Cd113", "In113", "Sn113", "Cd114", "Sn114", "In115", "Sn115",
    "Cd116", "Sn116", "Sn117", "Sn118", "Sn119", "Sn120", "Te120", "Sb121",
    "Sn122", "Te122", "Sn123", "Sb123", "Te123", "Sn124", "Sb124", "Te124",
    "Sn125", "Sb125", "Te125", "Sn126", "Sb126", "Te126", "Xe126", "I127",
    "Te128", "Xe128", "I129", "Xe129", "Te130", "I130", "Xe130", "I131",
    "Xe131", "Te132", "Xe132", "Ba132", "Xe133", "Cs133", "Ba133", "Xe134",
    "Cs134", "Ba134", "I135", "Xe135", "Cs135", "Ba135", "Xe136", "Cs136",
    "Ba136", "Cs137", "Ba137", "Ba138", "La138", "Ce138", "La139", "Ce139",
    "Ba140", "La140", "Ce140", "Ce141", "Pr141", "Ce142", "Pr142", "Nd142",
    "Ce143", "Pr143", "Nd143", "Ce144", "Nd144", "Nd145", "Nd146", "Nd147",
    "Pm147", "Sm147", "Nd148", "Pm148", "Sm148", "Pm149", "Sm149", "Nd150",
    "Sm150", "Pm151", "Sm151", "Eu151", "Sm152", "Eu152", "Gd152", "Sm153",
    "Eu153", "Gd153", "Sm154", "Eu154", "Gd154", "Eu155", "Gd155", "Eu156",
    "Gd156", "Eu157", "Gd157", "Gd158", "Dy158", "Tb159", "Gd160", "Tb160",
    "Dy160", "Dy161", "Dy162", "Dy163", "Dy164", "Er164", "Ho165", "Er166",
    "Er167", "Er168", "Tm168", "Tm169", "Er170", "Tm170"]


mats = {}

# Create He gas material for fuel pin gap
mats['He'] = openmc.Material(name='Helium')
mats['He'].set_density('g/cc', 0.0015981)
mats['He'].add_element('He', 1.0, 'ao')

# Create air material for instrument tubes
mats['Air'] = openmc.Material(name='Air')
mats['Air'].set_density('g/cc', 0.00616)
mats['Air'].add_element('O', 0.2095, 'ao')
mats['Air'].add_element('N', 0.7809, 'ao')
mats['Air'].add_element('Ar', 0.00933, 'ao')
mats['Air'].add_element('C', 0.00027, 'ao')

# Create inconel 718 material
mats['In'] = openmc.Material(name='Inconel')
mats['In'].set_density('g/cc', 8.2)
mats['In'].add_element('Si', 0.0035, 'wo')
mats['In'].add_element('Cr', 0.1896, 'wo')
mats['In'].add_element('Mn', 0.0087, 'wo')
mats['In'].add_element('Fe', 0.2863, 'wo')
mats['In'].add_element('Ni', 0.5119, 'wo')

# Create stainless steel 302
mats['SS302'] = openmc.Material(name='SS302')
mats['SS302'].set_density('g/cm3', 7.86)
mats['SS302'].add_element('Si', 0.01, 'wo')
mats['SS302'].add_element('Cr', 0.18, 'wo')
mats['SS302'].add_element('Mn', 0.02, 'wo')
mats['SS302'].add_element('Fe', 0.70, 'wo')
mats['SS302'].add_element('Ni', 0.09, 'wo')

# Create stainless steel material
mats['SS'] = openmc.Material(name='SS304')
mats['SS'].set_density('g/cc', 8.03)
mats['SS'].add_element('Si', 0.0060, 'wo')
mats['SS'].add_element('Cr', 0.1900, 'wo')
mats['SS'].add_element('Mn', 0.0200, 'wo')
mats['SS'].add_element('Fe', 0.6840, 'wo')
mats['SS'].add_element('Ni', 0.1000, 'wo')

# Create carbon steel material
mats['CS'] = openmc.Material(name='Carbon Steel')
mats['CS'].set_density('g/cc', 7.8)
mats['CS'].add_element('C', 0.00270, 'wo')
mats['CS'].add_element('Mn', 0.00750, 'wo')
mats['CS'].add_element('P', 0.00025, 'wo')
mats['CS'].add_element('S', 0.00025, 'wo')
mats['CS'].add_element('Si', 0.00400, 'wo')
mats['CS'].add_element('Ni', 0.00750, 'wo')
mats['CS'].add_element('Cr', 0.00350, 'wo')
mats['CS'].add_element('Mo', 0.00625, 'wo')
mats['CS'].add_element('V', 0.00050, 'wo')
mats['CS'].add_element('Nb', 0.00010, 'wo')
mats['CS'].add_element('Cu', 0.00200, 'wo')
mats['CS'].add_element('Ca', 0.00015, 'wo')
mats['CS'].add_element('B', 0.00003, 'wo')
mats['CS'].add_element('Ti', 0.00015, 'wo')
mats['CS'].add_element('Al', 0.00025, 'wo')
mats['CS'].add_element('Fe', 0.96487, 'wo')

# Create zircaloy 4 material
mats['Zr'] = openmc.Material(name='Zircaloy-4')
mats['Zr'].set_density('g/cc', 6.55)
mats['Zr'].add_element('O', 0.00125, 'wo')
mats['Zr'].add_element('Cr', 0.0010, 'wo')
mats['Zr'].add_element('Fe', 0.0021, 'wo')
mats['Zr'].add_element('Zr', 0.98115, 'wo')
mats['Zr'].add_element('Sn', 0.0145, 'wo')

# Create M5 alloy material
m5_niobium = 0.01    # http://publications.jrc.ec.europa.eu/repository/bitstream/JRC100644/lcna28366enn.pdf
m5_oxygen = 0.00135  # http://publications.jrc.ec.europa.eu/repository/bitstream/JRC100644/lcna28366enn.pdf
m5_density = 6.494   # 10.1039/C5DT03403E
mats['M5'] = openmc.Material(name='M5')
mats['M5'].add_element('Zr', 1.0 - m5_niobium - m5_oxygen)
mats['M5'].add_element('Nb', m5_niobium)
mats['M5'].add_element('O', m5_oxygen)
mats['M5'].set_density('g/cm3', m5_density)

# Create Ag-In-Cd control rod material
mats['AIC'] = openmc.Material(name='Ag-In-Cd')
mats['AIC'].set_density('g/cc', 10.16)
mats['AIC'].add_element('Ag', 0.80, 'wo')
mats['AIC'].add_element('In', 0.15, 'wo')
mats['AIC'].add_element('Cd', 0.05, 'wo')


#### Borated Water

# Concentration of boron at beginning of equilibrium cycle
boron_ppm = 1240  # ML17013A274, Figure 4.3-17

# Density of water
h2o_dens = water_density(core_average_temperature, system_pressure)

# Weight percent of natural boron in borated water
wB_Bh2o = boron_ppm * 1.0e-6

# Borated water density
rho_Bh2o = h2o_dens / (1 - wB_Bh2o)

# Compute weight percent of clean water in borated water
wh2o_Bh2o = 1.0 - wB_Bh2o

# Compute molecular mass of clean water
M_h2o = 2. * atomic_weight('H') + atomic_weight('O')

# Compute molecular mass of borated water
M_Bh2o = 1. / (wB_Bh2o / atomic_weight('B') + wh2o_Bh2o / M_h2o)

# Compute atom fractions of boron and water
aB_Bh2o = wB_Bh2o * M_Bh2o / atomic_weight('B')
ah2o_Bh2o = wh2o_Bh2o * M_Bh2o / M_h2o

# Compute atom fractions of hydrogen, oxygen
ah_Bh2o = 2.0 * ah2o_Bh2o
aho_Bh2o = ah2o_Bh2o

# Create borated water for coolant / moderator
mats['H2O'] = openmc.Material(name='Borated Water')
mats['H2O'].set_density('g/cc', rho_Bh2o)
mats['H2O'].add_element('B', aB_Bh2o, 'ao')
mats['H2O'].add_element('H', ah_Bh2o, 'ao')
mats['H2O'].add_element('O', aho_Bh2o, 'ao')
mats['H2O'].add_s_alpha_beta(name='c_H_in_H2O')


#### Borosilicate Glass

# CASMO weight fractions
wO_bsg = 0.5481
wAl_bsg = 0.0344
wSi_bsg = 0.3787
wB10_bsg = 0.0071
wB11_bsg = 0.0317

# Molar mass of borosilicate glass
M_bsg = 1.0 / (wO_bsg / atomic_weight('O') + wAl_bsg / atomic_weight('Al') +
               wSi_bsg /atomic_weight('Si') + wB10_bsg / atomic_mass('B10') +
               wB11_bsg / atomic_mass('B11'))

# Compute atom fractions for borosilicate glass
aO_bsg = wO_bsg * M_bsg / atomic_weight('O')
aAl_bsg = wAl_bsg * M_bsg / atomic_weight('Al')
aSi_bsg = wSi_bsg * M_bsg / atomic_weight('Si')
aB10_bsg = wB10_bsg * M_bsg / atomic_mass('B11')
aB11_bsg = wB11_bsg * M_bsg / atomic_mass('B10')
aB_bsg = aB10_bsg + aB11_bsg

# Create borosilicate glass material
mats['BSG'] = openmc.Material(name='Borosilicate Glass')
mats['BSG'].temperature = 300
mats['BSG'].set_density('g/cc', 2.26)
mats['BSG'].add_element('O', aO_bsg, 'ao')
mats['BSG'].add_element('Si', aSi_bsg, 'ao')
mats['BSG'].add_element('Al', aAl_bsg, 'ao')
mats['BSG'].add_nuclide('B10', aB10_bsg, 'ao')
mats['BSG'].add_nuclide('B11', aB11_bsg, 'ao')


#### Enriched UO2 Fuel

# Create 1.6% enriched UO2 fuel material
mat = openmc.Material(name='1.6% Enr. UO2 Fuel')
mat.temperature = 300
mat.set_density('g/cc', 10.31341)
mat.add_element('O', 2., 'ao')
mat.add_element('U', 1., 'ao', enrichment=1.61006)
mats['UO2 1.6 fresh'] = mat

# Create 2.4% enriched UO2 fuel material
mat = openmc.Material(name='2.4% Enr. UO2 Fuel')
mat.temperature = 300
mat.set_density('g/cc', 10.29748)
mat.add_element('O', 2., 'ao')
mat.add_element('U', 1., 'ao', enrichment=2.39993)
mats['UO2 2.4 fresh'] = mat

# Create 3.1% enriched UO2 fuel material
mat = openmc.Material(name='3.1% Enr. UO2 Fuel')
mat.temperature = 300
mat.set_density('g/cc', 10.30166)
mat.add_element('O', 2., 'ao')
mat.add_element('U', 1., 'ao', enrichment=3.10221)
mats['UO2 3.1 fresh'] = mat

# Depleted versions of 1.6%, 2.4%, 3.1% fuel
mat = openmc.Material(name='2.4% Enr. UO2 Fuel')
mat.temperature = 300
mat.set_density('g/cc', 10.29748)
mat.add_element('O', 2., 'ao')
mat.add_element('U', 1., 'ao', enrichment=2.39993)
for nuc in _DEPLETION_NUCLIDES:
    mat.add_nuclide(nuc, 1.0e-11)
mats['UO2 2.4 depleted'] = mat

mat = openmc.Material(name='1.6% Enr. UO2 Fuel')
mat.temperature = 300
mat.set_density('g/cc', 10.31341)
mat.add_element('O', 2., 'ao')
mat.add_element('U', 1., 'ao', enrichment=1.61006)
for nuc in _DEPLETION_NUCLIDES:
    mat.add_nuclide(nuc, 1.0e-11)
mats['UO2 1.6 depleted'] = mat

mat = openmc.Material(name='3.1% Enr. UO2 Fuel')
mat.temperature = 300
mat.set_density('g/cc', 10.30166)
mat.add_element('O', 2., 'ao')
mat.add_element('U', 1., 'ao', enrichment=3.10221)
for nuc in _DEPLETION_NUCLIDES:
    mat.add_nuclide(nuc, 1.0e-11)
mats['UO2 3.1 depleted'] = mat


# Construct a collection of Materials to export to XML

materials = openmc.Materials(mats.values())
