# Scrip to demonstrate the usage of USpek class with CSV input files
from pathlib import Path

from uspeckpy.uspek import USpek

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

# Defines mass energy transfer coefficients of air
my_mu_csv = 'data/mu_tr_rho.csv'

# Defines conversion coefficients
my_hk_csv = 'data/h_k_h_amb_10.csv'

# Defines relative uncertainty (k=1) on mass energy transfer coefficients of air 
my_mu_std = 0.01

# Create USpekPy object with given beam parameters, mass energy transfer coefficients of air and conversion coefficients
s = USpek(beam_parameters=my_beam, mass_transmission_coefficients=my_mu_csv,
          mass_transmission_coefficients_uncertainty=my_mu_std, conversion_coefficients=my_hk_csv)

# Runs simulation with a given number of iterations
df = s.simulate(simulations_number=3)

# Defines the output folder path where the simulation results will be saved
my_folder = 'output'

# Defines the name of the output file where the simulation results will be saved
my_output_csv = 'output.csv'

# Saves results to a CSV file
df.to_csv(Path(my_folder) / my_output_csv, index=True)
