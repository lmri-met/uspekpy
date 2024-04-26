import uspeckpy.digest as dg
import numpy as np
import pandas as pd
from scipy.interpolate import Akima1DInterpolator
from spekpy import Spek


class SpecWrapper(Spek):
    """
    A subclass of Spek representing a spectrum with additional methods for computation of mean quantities.

    Parameters:
    - kvp (float): The peak kilovoltage (kVp) of the x-ray tube.
    - th (float): The anode angle of the x-ray tube.

    Attributes:
    Inherits kvp and th attributes from Spek class.

    Methods:
    - get_mean_energy(): Compute the mean energy of the spectrum.
    - get_mean_kerma(mass_transmission_coefficients): Compute the mean kerma of the spectrum using mass transmission
      coefficients.
    - get_mean_conversion_coefficient(mass_transmission_coefficients, conversion_coefficients, angle=None):
      Compute the mean conversion coefficient of the spectrum using mass transmission coefficients and conversion
      coefficients.
    """
    def __init__(self, kvp, th):
        """
        Initialize SpecWrapper instance with peak kilovoltage anode angle of the x-ray tube.
        """
        Spek.__init__(self, kvp, th)

    def get_mean_energy(self):
        """
        Compute the mean energy of the spectrum.

        Returns:
        float: The mean energy of the spectrum.
        """
        # Get spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)

        # Compute mean energy
        return sum(fluence * energy) / fluence.sum()

    def get_mean_kerma(self, mass_transmission_coefficients):
        """
        Compute the mean kerma using mass transmission coefficients.

        Parameters:
        mass_transmission_coefficients (tuple): Tuple containing the energies and values of the mass transmission
        coefficients.

        Returns:
        float: The mean kerma computed.
        """
        # Get spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)

        # Unpack the energies and values of the mass transmission coefficients
        energy_mu, mu = dg.parse_mass_transmission_coefficients(mass_transmission_coefficients)

        # Interpolate mass attenuation coefficients for the spectrum energies in logarithmic scale
        interpolated_mu = interpolate(x=energy_mu, y=mu, new_x=energy)

        # Compute mean kerma
        return sum(fluence * energy * interpolated_mu) / fluence.sum()

    def get_mean_conversion_coefficient(self, mass_transmission_coefficients, conversion_coefficients, angle=None):
        """
        Compute the mean conversion coefficient using mass transmission coefficients and conversion coefficients.

        Parameters:
        mass_transmission_coefficients (tuple): Tuple containing the energies and values of the mass transmission
        coefficients.
        conversion_coefficients (tuple): Tuple containing the energies and values of the conversion coefficients.
        angle (float, optional): The irradiation angle at which conversion coefficients are calculated.

        Returns:
        float: The mean conversion coefficient computed.
        """
        # Get spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)

        # Unpack the energies and values of the mass transmission coefficients
        energy_mu, mu = dg.parse_mass_transmission_coefficients(mass_transmission_coefficients)

        # Unpack the energies and values of the conversion coefficients
        energy_hk, hk = dg.parse_conversion_coefficients(conversion_coefficients, angle)

        # Interpolate mass attenuation coefficients for the spectrum energies in logarithmic scale
        interpolated_mu = interpolate(x=energy_mu, y=mu, new_x=energy)

        # Interpolate conversion coefficients for the spectrum energies in logarithmic scale
        interpolated_hk = interpolate(x=energy_hk, y=hk, new_x=energy)

        # Compute mean conversion coefficient
        return sum(fluence * energy * interpolated_mu * interpolated_hk) / fluence.sum()


class USpek:
    def __init__(self, beam_parameters, mass_transmission_coefficients, mass_transmission_coefficients_uncertainty,
                 conversion_coefficients, angle=None):
        # Initialize USpekPy instance with beam parameters
        self.beam = beam_parameters
        self.mass_transmission_coefficients = dg.parse_mass_transmission_coefficients(mass_transmission_coefficients)
        self.conversion_coefficients = dg.parse_conversion_coefficients(conversion_coefficients, angle)
        self.mass_transmission_coefficients_uncertainty = mass_transmission_coefficients_uncertainty

    def _get_random_values(self):
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
        energy_mu = self.mass_transmission_coefficients[0]
        nominal_mu = self.mass_transmission_coefficients[1]
        mu_std = self.mass_transmission_coefficients_uncertainty
        random_mu = np.random.normal(loc=nominal_mu, scale=nominal_mu * mu_std)
        return kvp, th, filters, (energy_mu, random_mu)

    def _iteration(self):
        # Generate random beam parameters
        kvp, th, filters, mu_tr_rho = self._get_random_values()

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
        mean_kerma = spectrum.get_mean_kerma(mass_transmission_coefficients=mu_tr_rho)

        # Get mean conversion coefficient
        mean_hk = spectrum.get_mean_conversion_coefficient(mass_transmission_coefficients=mu_tr_rho,
                                                           conversion_coefficients=self.conversion_coefficients)

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
        # Print a message indicating the start of simulation
        print('Simulation')
        rows = []
        # For each iteration
        for iteration in range(simulations_number):
            # Print a message indicating the number of the current iteration
            print(f'Iteration number: {iteration + 1}')
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
    return float(np.random.normal(loc=loc, scale=loc * scale, size=1)[0])
