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
    def __init__(self, beam_parameters, mass_transmission_coefficients, conversion_coefficients):
        # Initialize USpekPy instance with beam parameters
        self.beam = beam_parameters
        self.mass_transmission_coefficients = mass_transmission_coefficients
        self.conversion_coefficients = conversion_coefficients

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
