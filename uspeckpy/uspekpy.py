import csv

import numpy as np
import pandas as pd
from scipy.interpolate import Akima1DInterpolator
from spekpy import Spek


# TODO: include uncertainty of mass transmission coeficients
class SpecWrapper(Spek):
    def __init__(self, kvp, th):
        # Initialize SpecWrapper as a subclass of Spek with kvp and th
        Spek.__init__(self, kvp, th)

    def get_mean_energy(self):
        # Method to compute the mean energy of the spectrum
        # Get spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)
        # Compute mean energy
        return sum(fluence * energy) / fluence.sum()

    def get_mean_kerma(self, mass_transmission_coefficients):
        # Method to compute the mean kerma using mass transmission coefficients
        # Get spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)
        # Unpack energy_mu and mu from mass_transmission_coefficients
        energy_mu, mu = mass_transmission_coefficients
        # Interpolate mass attenuation coefficients for the given energy spectrum in logarithmic scale
        interpolated_mu = interpolate(x=energy_mu, y=mu, new_x=energy)
        # Compute mean kerma
        return sum(fluence * energy * interpolated_mu) / fluence.sum()

    def get_mean_conversion_coefficient(self, mass_transmission_coefficients, conversion_coefficients):
        # Method to compute the mean conversion coefficient using mass transmission coefficients
        # Get spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)
        # Unpack energy_mu and mu from mass_transmission_coefficients
        energy_mu, mu = mass_transmission_coefficients
        # Unpack energy_hk and hk from mono_energetic_conversion_coefficients
        energy_hk, hk = conversion_coefficients
        # Interpolate mass attenuation coefficients for the given energy spectrum in logarithmic scale
        interpolated_mu = interpolate(x=energy_mu, y=mu, new_x=energy)
        # Interpolate mono energetic conversion coefficients for the given energy spectrum in logarithmic scale
        interpolated_hk = interpolate(x=energy_hk, y=hk, new_x=energy)
        # Compute mean conversion coefficient
        return sum(fluence * energy * interpolated_mu * interpolated_hk) / fluence.sum()


class USpek:
    def __init__(self, beam_parameters, mass_transmission_coefficients, conversion_coefficients, angle=None):
        # Initialize USpekPy instance with beam parameters
        self.beam = beam_parameters
        self._set_mass_transmission_coefficients(mass_transmission_coefficients)
        self._set_conversion_coefficients(conversion_coefficients, angle)

    def _set_mass_transmission_coefficients(self, coefficients):
        if is_tuple_of_two_arrays(coefficients):
            self.mass_transmission_coefficients = coefficients
        elif is_csv_with_two_columns(coefficients):
            self.mass_transmission_coefficients = parse_mass_transmission_coefficients(coefficients)
        else:
            raise ValueError(f"Unsupported mass transmission coefficients format. Only a tuple of two numpy arrays and "
                             f"a CSV file with two columns are supported.")

    def _set_conversion_coefficients(self, coefficients, angle):
        if is_tuple_of_two_arrays(coefficients):
            self.conversion_coefficients = coefficients
        elif is_valid_csv(coefficients):
            self.conversion_coefficients = parse_conversion_coefficients(coefficients, angle)
        else:
            raise ValueError(f"Unsupported conversion coefficients format. Only a tuple of two numpy arrays and "
                             f"a CSV file with two or more columns are supported.")

    def _get_beam_parameters(self):
        # Method to generate random beam parameters based on the values and uncertainties of the beam parameters
        kvp = random_normal(loc=self.beam['kVp'][0], scale=self.beam['kVp'][1])
        th = random_normal(loc=self.beam['th'][0], scale=self.beam['th'][1])
        filters = [
            ('Air', random_normal(loc=self.beam['Air'][0], scale=self.beam['Air'][1])),
            ('Al', random_uniform(loc=self.beam['Al'][0], scale=self.beam['Al'][1])),
            ('Cu', random_uniform(loc=self.beam['Cu'][0], scale=self.beam['Cu'][1])),
            ('Sn', random_uniform(loc=self.beam['Sn'][0], scale=self.beam['Sn'][1])),
            ('Pb', random_uniform(loc=self.beam['Pb'][0], scale=self.beam['Pb'][1])),
            ('Be', random_uniform(loc=self.beam['Be'][0], scale=self.beam['Be'][1])),
        ]
        return kvp, th, filters

    def _iteration(self):
        # Generate random beam parameters
        kvp, th, filters = self._get_beam_parameters()

        # Initialize an SpeckWrapper object and add filters
        spectrum = SpecWrapper(kvp=kvp, th=th)
        spectrum.multi_filter(filters)

        # Calculate half-value layers for aluminum and copper
        hvl1_al = spectrum.get_hvl1()
        hvl2_al = spectrum.get_hvl2()
        hvl1_cu = spectrum.get_hvl1(matl='Cu')
        hvl2_cu = spectrum.get_hvl2(matl='Cu')

        # Get mean energy
        mean_energy = spectrum.get_mean_energy()

        # Get mean kerma
        mean_kerma = spectrum.get_mean_kerma(mass_transmission_coefficients=self.mass_transmission_coefficients)

        # Get mean conversion coefficient
        mean_hk = spectrum.get_mean_conversion_coefficient(
            mass_transmission_coefficients=self.mass_transmission_coefficients,
            conversion_coefficients=self.conversion_coefficients
        )

        return (kvp, th, filters[0][1], filters[1][1], filters[2][1], filters[3][1], filters[4][1], filters[5][1],
                hvl1_al, hvl2_al, hvl1_cu, hvl2_cu, mean_energy, mean_kerma, mean_hk)

    @staticmethod
    def _get_mean_quantities(rows):
        columns = ['#', 'kVp', 'th', 'Air', 'Al', 'Cu', 'Sn', 'Pb', 'Be', 'HVL1 Al', 'HVL2 Al', 'HVL1 Cu', 'HVL2 Cu',
                   'Mean energy', 'Mean kerma', 'Mean conv. coefficient.']
        df = pd.DataFrame(data=rows, columns=columns)

        # Calculate means, standard deviations, and relative uncertainties for the simulation results
        means = df.iloc[:, 1:].mean()
        standard_deviations = df.iloc[:, 1:].std(ddof=0)
        relative_uncertainties = standard_deviations / means

        # Append means, standard deviations, and relative uncertainties to the DataFrame
        means = ['Mean'] + list(means)
        standard_deviations = ['Standard deviation'] + list(standard_deviations)
        relative_uncertainties = ['Relative uncertainty'] + list(relative_uncertainties)
        data = [means, standard_deviations, relative_uncertainties]
        return pd.concat(objs=[df, pd.DataFrame(data=data, columns=columns)], ignore_index=True)

    def simulate(self, simulations_number):
        rows = []
        # For each iteration
        for iteration in range(simulations_number):
            # Print a message indicating the number of the current iteration
            print(f'Iteration number: {iteration}')
            rows.append([f'Iteration {iteration + 1}'] + list(self._iteration()))

        # Build results DataFrame
        return self._get_mean_quantities(rows=rows)


def interpolate(x, y, new_x):
    # Create an Akima1DInterpolator object with logarithmic x and y values
    interpolator = Akima1DInterpolator(x=np.log(x), y=np.log(y))

    # Interpolate new y values for given new_x using the Akima1DInterpolator
    new_y = np.exp(interpolator(x=np.log(new_x)))

    # Replace any NaN values with zeros in the interpolated y values
    return np.nan_to_num(new_y, nan=0)


def random_uniform(loc, scale):
    low = loc * (1 - scale * np.sqrt(3))
    high = loc * (1 - scale * np.sqrt(3))
    return float(np.random.uniform(low=low, high=high, size=1)[0])


def random_normal(loc, scale):
    # TODO: check uncertainty type for np.random.normal()
    return float(np.random.normal(loc=loc, scale=scale, size=1)[0])


def parse_mass_transmission_coefficients(file_path):
    # Get mass transmission coefficients in the format required by SpekWrapper (tuple of two numpy arrays) from CSV file

    # The mass transmission coefficients CSV file should have only 2 columns, the first one for the energy and the
    # second one for the mass transmission coefficients.

    # Load CSV file into numpy array
    array2d = np.genfromtxt(file_path, delimiter=',', skip_header=1, unpack=True)

    # Build tuple of mass transmission coefficients in the format required by SpekWrapper (tuple of two numpy arrays)
    return array2d[0], array2d[1]


def parse_conversion_coefficients(file_path, irradiation_angle):
    # Get conversion coefficients in the format required by SpekWrapper (tuple of two numpy arrays) from CSV file

    # Read CSV file into DataFrame
    df = pd.read_csv(file_path)

    # Get the energies of the mono-energetic conversion coefficients (always the first column of the Dataframe)
    energies = df.iloc[:, 0]

    # If the conversion coefficients Dataframe has only 2 columns, the first is the energy and the second is the
    # conversion coefficients.
    if df.shape[1] == 2:
        # Get the values of the mono-energetic conversion coefficients (the second column of the Dataframe)
        values = df.iloc[:, 1]
    # If the conversion coefficients Dataframe has more than 2 columns, then it contains conversion coefficients for
    # more than one irradiation angle.
    else:  # TODO: check with multicolumn input file
        # Get the label of the column which contains the irradiation angle
        column_label = next((label for label in df.columns if irradiation_angle in label), None)

        # Get the values of the mono-energetic conversion coefficients (the second column of the Dataframe)
        values = df.loc[:, column_label]

    # Build tuple of conversion coefficients in the format required by SpekWrapper (tuple of two numpy arrays)
    return energies, values


def is_tuple_of_two_arrays(arg):
    if not isinstance(arg, tuple):
        return False
    if len(arg) != 2:
        return False
    if not isinstance(arg[0], np.ndarray) or not isinstance(arg[1], np.ndarray):
        return False
    return True


def is_csv_with_two_columns(filename):
    try:
        with open(filename, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                if len(row) != 2:
                    return False
            return True
    except Exception as e:
        return False


def is_valid_csv(filename):
    try:
        with open(filename, 'r') as file:
            csv.reader(file)
        return True
    except Exception as e:
        return False
