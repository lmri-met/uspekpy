# Scripto to demonstrate the usage of SpekWrapper class with CSV input files
from uspekpy import SpekWrapper

# Define x-ray beam parameters for radiation quality N-60 (filter thickness, peak kilovoltage and anode angle)
my_filters = [
  ('Al', 4),  # mm
  ('Cu', 0.6),  # mm
  ('Sn', 0),  # mm
  ('Pb', 0),  # mm
  ('Be', 0),  # mm
  ('Air', 1000)  # mm
]
my_kvp = 60  # kV
my_th = 10  # deg

# Define path to CSV file containing mass energy transfer coefficients of air in terms of the energy
my_mu_csv = 'data/mu_tr_rho.csv'

# Define path to CSV file containing the monoenergetic air kerma-to-dose-equivalent conversion coefficients for H*(10)
my_hk_csv = 'data/h_k_h_amb_10.csv'

# Initialize an SpeckWrapper object and add filters
spectrum = SpekWrapper(kvp=my_kvp, th=my_th)
spectrum.multi_filter(my_filters)

# Calculate half-value layers for aluminum and copper
hvl1_al = spectrum.get_hvl1()
hvl2_al = spectrum.get_hvl2()
hvl1_cu = spectrum.get_hvl1(matl='Cu')
hvl2_cu = spectrum.get_hvl2(matl='Cu')

# Calculate mean energy
mean_energy = spectrum.get_mean_energy()

# Calculate mean air kerma
mean_kerma = spectrum.get_air_kerma(mass_transfer_coefficients=my_mu_csv)

# Calculate mean air kerma-to-dose-equivalent conversion coefficient for H*(10)
mean_hk = spectrum.get_mean_conversion_coefficient(
  mass_transfer_coefficients=my_mu_csv, conversion_coefficients=my_hk_csv)

# Print results
print(f'First HVL for Al: {hvl1_al} mm')
print(f'Second HVL for Al: {hvl2_al} mm')
print(f'First HVL for Cu: {hvl1_cu} mm')
print(f'Second HVL for Cu: {hvl2_cu} mm')
print(f'Mean energy: {mean_energy} keV')
print(f'Air kerma: {mean_kerma} keV/g')
print(f'Mean conversion coefficient for H*(10): {mean_hk} Sv/Gy')