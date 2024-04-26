# Scrip to demonstrate the usage of USpek class with CSV input files

from uspeckpy.uspekpy import USpek

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
df = s.simulate(simulations_number=3)
