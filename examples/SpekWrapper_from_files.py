# Scripto to demonstrate the usage of SpekWrapper class with CSV input files
from uspeckpy.wrapper import SpekWrapper

# Defines x-ray beam parameters
my_filters = [
    ('Al', 4),
    ('Cu', 0.6),
    ('Sn', 0),
    ('Pb', 0),
    ('Be', 0),
    ('Air', 1000)
]

# Defines mass energy transfer coefficients of air
my_mu_csv = 'data/mu_tr_rho.csv'

# Defines conversion coefficients
my_hk_csv = 'data//h_k_h_amb_10.csv'

# Initialises an SpeckWrapper object and add filters
spectrum = SpekWrapper(kvp=60, th=20)
spectrum.multi_filter(my_filters)

# Calculate half-value layers for aluminum and copper using SpekPy method
hvl1_al = spectrum.get_hvl1()
hvl2_al = spectrum.get_hvl2()
hvl1_cu = spectrum.get_hvl1(matl='Cu')
hvl2_cu = spectrum.get_hvl2(matl='Cu')

# Gets mean energy
mean_energy = spectrum.get_mean_energy()
mean_energy_SpekPy = spectrum.get_emean()

# Gets mean air kerma
mean_kerma_SpekPy = spectrum.get_kerma()
mean_kerma = spectrum.get_mean_kerma(mass_transmission_coefficients=my_mu_csv)

# Gets mean conversion coefficient
mean_hk = spectrum.get_mean_conversion_coefficient(mass_transmission_coefficients=my_mu_csv,
                                                   conversion_coefficients=my_hk_csv)
