# Scrip to demonstrate the usage of USpek class with CSV input files
from uspekpy import USpek

# Define values and relative uncertainty (k=1) of x-ray beam parameters for radiation quality N-60
# (filter thickness, peak kilovoltage and anode angle)
my_beam = {
    'kVp': (60, 0.01),  # mm, fraction of one
    'th': (20, 0.01),  # mm, fraction of one
    'Al': (4, 0.01),  # mm, fraction of one
    'Cu': (0.6, 0.01),  # mm, fraction of one
    'Sn': (0, 0),  # mm, fraction of one
    'Pb': (0, 0),  # mm, fraction of one
    'Be': (0, 0),  # mm, fraction of one
    'Air': (1000, 0.01)  # mm, fraction of one
}

# Define path to CSV file containing mass energy transfer coefficients of air in terms of the energy
my_mu_csv = 'data/mu_tr_rho.csv'

# Define path to CSV file containing the monoenergetic air kerma-to-dose-equivalent conversion coefficients for H*(10)
my_hk_csv = 'data/h_k_h_amb_10.csv'

# Define relative uncertainty (k=1) of mass energy transfer coefficients of air
my_mu_std = 0.01  # fraction of one

# Initialize a USpekPy object with the defined beam parameters, mass energy transfer coefficients of air
# and monoenergetic air kerma-to-dose-equivalent conversion coefficients for H*(10)
s = USpek(beam_parameters=my_beam, mass_transfer_coefficients=my_mu_csv,
          mass_transfer_coefficients_uncertainty=my_mu_std, conversion_coefficients=my_hk_csv)

# Run simulation with a given number of iterations
df = s.simulate(simulations_number=3)

# Define the path of the output file where the simulation results will be saved
my_output_csv = 'output/output.csv'

# Save results to a CSV file
df.to_csv(my_output_csv, index=True)

# Print results
print(df.to_string())
