# Scripto to demonstrate the functionality of USpek class

import numpy as np

from uspeckpy.uspekpy import USpek

# Using USpek without input files
# ----------------------------------------------------------------------------------------------------------------------

# Define x-ray beam parameters
my_beam = {
    'kVp': (60, 0.01),
    'th': (20, 0.01),
    'Al': (4, 0.01),
    'Cu': (0.6, 0.01),
    'Sn': (0, 0),
    'Pb': (0, 0),
    'Be': (0, 0),
    'Air': (1000, 0.01)
}

# Define mass transmission coefficients
my_mu = (
    np.array(
        [1.0, 1.1726, 1.25, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0, 3.2063, 3.206301, 3.22391, 3.25051, 3.5, 3.61881, 4.0,
         5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 14.0, 15.0, 17.5, 20.0, 25.0, 28.6633, 30.0, 35.0, 40.0, 50.0, 60.0,
         70.0, 80.0, 90.0, 100.0, 125.0, 140.0, 150.0, 175.0, 187.083, 200.0, 250.0, 300.0, 324.037, 350.0, 386.867,
         400.0, 474.342, 500.0, 574.456, 600.0, 673.537, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1558.93,
         1750.0, 1870.83, 2000.0, 2345.21, 2500.0, 3000.0, 3240.37, 3500.0, 4000.0, 4500.0, 5000.0, 6000.0, 6480.74,
         7000.0, 8000.0, 9000.0]),
    np.array(
        [3487.7, 2271.66, 1907.85, 1396.25, 1152.41, 746.85, 510.495, 267.712, 156.677, 128.597, 139.322, 139.22,
         136.749, 110.244, 99.9467, 74.3863, 38.3165, 22.1387, 13.8638, 9.20954762864289, 6.40746485225397,
         4.62692014996195, 2.30743153037744, 1.61505075253763, 1.30130979037554, 0.801145676760343,
         0.525990443101652, 0.262079759969996, 0.17213944380102, 0.150643910752649, 0.096166283094939,
         0.067006659158694, 0.040350724533314, 0.030060380204255, 0.025736225862943, 0.023919178948042,
         0.023273564441493, 0.023169137685326, 0.023905175732908, 0.024501233602389, 0.024926643691358,
         0.025877217107796, 0.026274062287156, 0.026655131026205, 0.027882238669658, 0.028683812435718,
         0.028972447929017, 0.029148771809045, 0.029415131704691, 0.029490249846129, 0.029698178702551,
         0.029685041261757, 0.029603764611649, 0.029549718949812, 0.029405101930449, 0.029197657857256,
         0.028890614121642, 0.028390536220701, 0.027936798166738, 0.026719405747041, 0.02562121786982,
         0.025359438553807, 0.02463108715039, 0.024240842040503, 0.023753744772526, 0.022671649894621,
         0.022309768351513, 0.021132688079239, 0.020639133403572, 0.020121558259145, 0.019347638398287,
         0.018700741656237, 0.018185053120065, 0.017326066208925, 0.016966042150857, 0.016690444038446,
         0.016199645451199, 0.015809252632524])
)

# Define conversion coefficients
my_hk = (
    np.array([7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 500,
              600, 800, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000]),
    np.array([1.2e-05, 9.5e-05, 0.00145, 0.008, 0.0331, 0.0737, 0.127, 0.19, 0.26, 0.326, 0.395, 0.466, 0.538, 0.61,
              1.1, 1.47, 1.67, 1.74, 1.72, 1.65, 1.49, 1.4, 1.31, 1.26, 1.23, 1.21, 1.19, 1.17, 1.15, 1.13, 1.12,
              1.11, 1.19, 1.09, 1.08, 1.06])
)

# Define mass transmission coefficients relative uncertainty (k=1)
my_mu_std = 0.01

# Create USpekPy object with given beam parameters, mass transmission coefficients and conversion coefficients
s = USpek(beam_parameters=my_beam, mass_transmission_coefficients=my_mu,
          mass_transmission_coefficients_uncertainty=my_mu_std, conversion_coefficients=my_hk)

# Run simulation with a given number of iterations
df1 = s.simulate(simulations_number=3)

# Using USpek with CSV input files
# ----------------------------------------------------------------------------------------------------------------------

# Define x-ray beam parameters
my_beam = {
    'kVp': (60, 0.01),
    'th': (20, 0.01),
    'Al': (4, 0.01),
    'Cu': (0.6, 0.01),
    'Sn': (0, 0),
    'Pb': (0, 0),
    'Be': (0, 0),
    'Air': (1000, 0.01)
}

# Define mass transmission coefficients
my_mu_csv = 'data/input/mu_tr_rho.csv'

# Define conversion coefficients
my_hk_csv = 'data/input/h_k_h_amb_10.csv'

# Define mass transmission coefficients relative uncertainty (k=1)
my_mu_std = 0.01

# Create USpekPy object with given beam parameters, mass transmission coefficients and conversion coefficients
s = USpek(beam_parameters=my_beam, mass_transmission_coefficients=my_mu_csv,
          mass_transmission_coefficients_uncertainty=my_mu_std, conversion_coefficients=my_hk_csv)

# Run simulation with a given number of iterations
df2 = s.simulate(simulations_number=3)
