# For compatibility with Python2 #
from __future__ import print_function, division, absolute_import
##################################
import numpy as np

# Physical constants
electron_mass = 511.0  # [keV]
fine_structure = 1.0 / 137.0  # The fine structure constant
electron_radius = 2.81794e-13  # Classical electron radius [cm]
avogadros_number = 6.02e23  # [atoms/mole]
electron_charge = 1.602e-19  # [C]
detour_factor = 2.0  # Reciprocal of detour factor for an electron in diffusion
geometry_factor = 1.0 / (4.0 * np.pi)  # Geometry factor (converts to per unit
# ... solid angle)

# Material constants
atomic_number_tungsten = 74  # Atomic number Tungsten
atomic_weight_tungsten = 183.8  # Atomic weight Tungsten [g]
density_tungsten = 19.3  # Density Tungsten [g/cm^3]
phi0 = atomic_number_tungsten ** 2 * electron_radius ** 2 * fine_structure  #
density_air = 1.205E-03  # Density air [g/cm^3]

# Edge energies for Characteristic Emission 
energy_K_edge = 69.5 
energy_L_1_edge = 12.1
energy_L_2_edge = 11.5
energy_L_3_edge = 10.2


# nK - ratio of K-photons escaping to num of brems emitted above K-edge
# ... nK is related to SpekCalc's P parameter (see Ref. [2] in SpekModel.py, 
# ... equations 17 to 20): 0.5*(1+P)*fk=nK
# nL1/L2/L3 - ratio of L-photons escaping to num of brems emitted above L-edge
# nbr - bremsstrahlung normalization factor
# scale - if a number is specified it is intepreted as keV and electron data is
# ... extrapolated from this this value
# brxs - source of bremsstrahlung cross-section: 'nist' or 'mebh' (modified 
# ... Elwet-Bethe-Heitler)
# depth - maximum penetration depth based on Thomson-Whiddington range ('tw') 
# ... or the CSDA range ('csda')
model_param_default = {'nbr': 0.9279, 
                       'nL1': 0.5832, 
                       'nL2': 1.1107,
                       'nL3': 2.8486,
                       'nK': 0.3734,
                       'scale': None,
                       'brxs': 'nist',
                       'depth': 'csda'}

model_param_legacy = {'nbr':  0.7046, 
                      'nL1': 0., 
                      'nL2': 0., 
                      'nL3': 0., 
                      'nK': 0.4921,
                      'scale': 100,
                      'brxs': 'mebh',
                      'depth': 'tw'}
# The default SpekCalc value of nbr (i.e. nf) was 0.68. The default SpekCalc 
# ... value of nK was 0.5*(1+0.33)*(4.4-1.0)/4.4 = 0.514
# The "legacy" values are different due to the different numerical 
# ... implementation here.
# The precise "legacy" values were selected to agree closely with the 
# ... predictions of SpekCalc
# The precise "default" values were selected to agree closely with the Monte
# ...  Carlo simulations (BEAMnrc) 

# I/O constants
dir_data = 'data'
dir_matl_def = 'matl_def'
dir_matl_usr = 'matl_usr'
dir_state_def = 'state_def'
dir_state_usr = 'state_usr'
dir_tables = 'tables'
extension_matl_composition = '.comp'
extension_state_file = '.state'
extension_export_spectrum_file = '.spk'
extension_data_file = '.dat'
file_brem_data = 'nist_brem.dat'
file_pe_data = 'pe.dat'
file_ne_data = 'ne.dat'
file_range_data = 'erange.dat'
file_atomic_weight_data = 'atwts.dat'
file_line_data = 'lines.dat'

# Conversion factors
conversion_MeV2keV = 1.0e3
conversion_A2mA = 1.0e-3
conversion_mm2cm = 1.0e-1
conversion_keV2J = 1.6022e-16
conversion_per_g2per_kg = 1.0e3
conversion_keV2ev = 1.0e3
conversion_Gy2uGy = 1.0e6
conversion_g2mg = 1.0e3

# Units in calculation
units_hvl = 'mm'
units_kerma = 'uGy'
units_norm_kerma = 'uGy mAs^-1'
units_energy = 'keV'
units_boolean = 'bool'
units_flu = 'Photons cm^-2'
units_flu_norm_mas = 'Photons cm^-2 mAs^-1'
units_flu_norm_kerma = \
    lambda ref_kerma: 'Photons cm^-2 (%.2f uGy)^-1' % ref_kerma
units_flu_norm_flu = \
    lambda ref_flu: 'Photons cm^-2 (%.2f Photons cm^-2)^-1' % ref_flu
