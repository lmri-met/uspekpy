import numpy as np
import pandas as pd

from uspeckpy.wrapper import SpecWrapper, parse_mass_transmission_coefficients, parse_conversion_coefficients


class USpek:
    """
    Class for computing mean quantities of an x-ray spectrum with uncertainties using Monte Carlo techniques.

    This class facilitates Monte Carlo simulations to compute mean quantities, standard deviations, and relative
    uncertainties of an x-ray spectrum based on provided beam parameters and coefficients.

    Parameters:
    - beam_parameters (dict): Dictionary containing beam parameters and their uncertainties.
        The 'beam_parameters' dictionary should contain the following key-value pairs:
        - 'kVp': Tuple containing the peak kilovoltage (kVp) value and its uncertainty.
        - 'th': Tuple containing the anode angle value and its uncertainty.
        - 'Al': Tuple containing the thickness of the Aluminum filter and its uncertainty.
        - 'Cu': Tuple containing the thickness of the Copper filter and its uncertainty.
        - 'Sn': Tuple containing the thickness of the Tin filter and its uncertainty.
        - 'Pb': Tuple containing the thickness of the Lead filter and its uncertainty.
        - 'Be': Tuple containing the thickness of the Beryllium filter and its uncertainty.
        - 'Air': Tuple containing the thickness of the Air filter and its uncertainty.
    - mass_transmission_coefficients (tuple): Tuple containing the energies and values of the mass transmission
    coefficients.
    - mass_transmission_coefficients_uncertainty (float): The uncertainty associated with the mass transmission
    coefficients.
    - conversion_coefficients (tuple): Tuple containing the energies and values of the conversion coefficients.
    - angle (float, optional): The irradiation angle at which conversion coefficients are calculated.

    Methods:
    - simulate(simulations_number): Perform simulations and return a DataFrame containing mean quantities.

    Attributes:
    - beam (dict): Beam parameters and their uncertainties.
    - mass_transmission_coefficients (tuple): Energies and values of mass transmission coefficients.
    - mass_transmission_coefficients_uncertainty (float): Uncertainty associated with mass transmission coefficients.
    - conversion_coefficients (tuple): Energies and values of conversion coefficients.
    - angle (float): Irradiation angle for calculating conversion coefficients.
    """

    def __init__(self, beam_parameters, mass_transmission_coefficients, mass_transmission_coefficients_uncertainty,
                 conversion_coefficients, angle=None):
        """
        Initialize USpek instance with beam parameters and coefficients.

        Args:
        - beam_parameters (dict): Dictionary containing beam parameters and their uncertainties.
            The 'beam_parameters' dictionary should contain the following key-value pairs:
            - 'kVp': Tuple containing the peak kilovoltage (kVp) value and its uncertainty.
            - 'th': Tuple containing the anode angle value and its uncertainty.
            - 'Al': Tuple containing the thickness of the Aluminum filter and its uncertainty.
            - 'Cu': Tuple containing the thickness of the Copper filter and its uncertainty.
            - 'Sn': Tuple containing the thickness of the Tin filter and its uncertainty.
            - 'Pb': Tuple containing the thickness of the Lead filter and its uncertainty.
            - 'Be': Tuple containing the thickness of the Beryllium filter and its uncertainty.
            - 'Air': Tuple containing the thickness of the Air filter and its uncertainty.
        - mass_transmission_coefficients (tuple): Tuple containing the energies and values of the mass transmission
        coefficients.
        - mass_transmission_coefficients_uncertainty (float): The uncertainty associated with the mass transmission
        coefficients.
        - conversion_coefficients (tuple): Tuple containing the energies and values of the conversion coefficients.
        - angle (float, optional): The irradiation angle at which conversion coefficients are calculated.
        """
        self.beam = beam_parameters
        self.mass_transmission_coefficients = parse_mass_transmission_coefficients(mass_transmission_coefficients)
        self.conversion_coefficients = parse_conversion_coefficients(conversion_coefficients, angle)
        self.mass_transmission_coefficients_uncertainty = mass_transmission_coefficients_uncertainty

    def _get_random_values(self):
        """
        Generate random beam parameters and mass transmission coefficients based on their values and uncertainties.

        This method generates random beam parameters and mass transmission coefficients necessary for Monte Carlo
        simulation based on their provided values and uncertainties. It generates random values for peak kilovoltage
        and anode angle (th) beam parameters, as well as for different beam filters such as Air, Al, Cu, Sn, Pb, and Be.
        Additionally, it generates a random value for the mass transmission coefficient with the associated uncertainty.

        Returns:
        tuple: A tuple containing random beam parameters and mass transmission coefficients.
            - kvp (float): Random peak kilovoltage beam parameter.
            - th (float): Random anode angle beam parameter.
            - filters (list): List of tuples containing random values for different beam filters.
                Each tuple contains the name of the filter and its associated random value.
            - mu_tr_rho (tuple): Tuple containing the energy and random mass transmission coefficient value.
                - energy_mu (array_like): Energies for mass transmission coefficients.
                - random_mu (float): Random mass transmission coefficient value.
        """
        # Generate random peak kilovoltage beam parameter
        kvp = random_normal(loc=self.beam['kVp'][0], scale=self.beam['kVp'][1])

        # Generate random anode angle beam parameter
        th = random_normal(loc=self.beam['th'][0], scale=self.beam['th'][1])

        # Generate random values for different beam filters
        filters = [
            ('Air', random_normal(loc=self.beam['Air'][0], scale=self.beam['Air'][1])),
            ('Al', random_uniform(loc=self.beam['Al'][0], scale=self.beam['Al'][1])),
            ('Cu', random_uniform(loc=self.beam['Cu'][0], scale=self.beam['Cu'][1])),
            ('Sn', random_uniform(loc=self.beam['Sn'][0], scale=self.beam['Sn'][1])),
            ('Pb', random_uniform(loc=self.beam['Pb'][0], scale=self.beam['Pb'][1])),
            ('Be', random_uniform(loc=self.beam['Be'][0], scale=self.beam['Be'][1])),
        ]

        # Get energy and nominal value of mass transmission coefficients
        energy_mu = self.mass_transmission_coefficients[0]
        nominal_mu = self.mass_transmission_coefficients[1]

        # Get uncertainty associated with mass transmission coefficients
        mu_std = self.mass_transmission_coefficients_uncertainty

        # Generate random mass transmission coefficient
        random_mu = np.random.normal(loc=nominal_mu, scale=nominal_mu * mu_std)

        # Return tuple containing generated random values for beam parameters and mass transmission coefficients
        return kvp, th, filters, (energy_mu, random_mu)

    def _iteration(self):
        """
        Perform a single iteration of the Monte Carlo simulation.

        This method performs a single iteration of the Monte Carlo simulation by generating random beam parameters
        and simulating a spectrum with those parameters. It then calculates various quantities such as half-value
        layers for aluminum and copper, mean energy, mean kerma, and mean conversion coefficient of the spectrum.

        Returns:
        tuple: A tuple containing results of a single iteration.
            - kvp (float): Random peak kilovoltage beam parameter.
            - th (float): Random anode angle beam parameter.
            - Air (float): Random value for the Air filter thickness.
            - Al (float): Random value for the Aluminum filter thickness.
            - Cu (float): Random value for the Copper filter thickness.
            - Sn (float): Random value for the Tin filter thickness.
            - Pb (float): Random value for the Lead filter thickness.
            - Be (float): Random value for the Beryllium filter thickness.
            - hvl1_al (float): Half-value layer (HVL) for aluminum.
            - hvl2_al (float): Second HVL for aluminum.
            - hvl1_cu (float): HVL for copper.
            - hvl2_cu (float): Second HVL for copper.
            - mean_energy (float): Mean energy of the spectrum.
            - mean_kerma (float): Mean kerma calculated using mass transmission coefficients.
            - mean_hk (float): Mean conversion coefficient calculated using mass transmission and conversion
              coefficients.
        """
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
        """
        Calculate means, standard deviations, and relative uncertainties for the simulation results.

        This method takes a list of simulation results and calculates the mean values, standard deviations,
        and relative uncertainties for each parameter based on the simulations.

        Args:
        rows (list): List of simulation results, where each element represents a single simulation and contains
            values for parameters such as peak kilovoltage (kVp), anode angle (th), filter thicknesses, half-value
            layers (HVLs) for aluminum and copper, mean energy, mean kerma, and mean conversion coefficient.

        Returns:
        pandas.DataFrame: DataFrame containing mean quantities, standard deviations, and relative uncertainties
            for the simulation results. The DataFrame has the following columns:
            - # (int): Index column representing simulation number.
            - kVp (float): Peak kilovoltage.
            - th (float): Anode angle.
            - Air (float): Thickness of the Air filter.
            - Al (float): Thickness of the Aluminum filter.
            - Cu (float): Thickness of the Copper filter.
            - Sn (float): Thickness of the Tin filter.
            - Pb (float): Thickness of the Lead filter.
            - Be (float): Thickness of the Beryllium filter.
            - HVL1 Al (float): Half-value layer (HVL) for Aluminum.
            - HVL2 Al (float): Second HVL for Aluminum.
            - HVL1 Cu (float): HVL for Copper.
            - HVL2 Cu (float): Second HVL for Copper.
            - Mean energy (float): Mean energy of the spectrum.
            - Mean kerma (float): Mean kerma calculated using mass transmission coefficients.
            - Mean conv. coefficient. (float): Mean conversion coefficient calculated using mass transmission and
                conversion coefficients.
            Additionally, it includes rows for mean values, standard deviations, and relative uncertainties
            of the simulation results.
        """
        # Define column names for the DataFrame
        columns = ['#', 'kVp', 'th', 'Air', 'Al', 'Cu', 'Sn', 'Pb', 'Be', 'HVL1 Al', 'HVL2 Al', 'HVL1 Cu', 'HVL2 Cu',
                   'Mean energy', 'Mean kerma', 'Mean conv. coefficient.']

        # Create DataFrame with simulation results
        df = pd.DataFrame(data=rows, columns=columns)

        # Calculate means, standard deviations, and relative uncertainties for the simulation results
        # Calculate mean values for each column except the first
        means = df.iloc[:, 1:].mean()

        # Calculate standard deviations for each column except the first
        standard_deviations = df.iloc[:, 1:].std(ddof=0)

        # Calculate relative uncertainties for each column except the first
        relative_uncertainties = standard_deviations / means

        # Add 'Mean' label to mean values
        means = ['Mean'] + list(means)

        # Add 'Standard deviation' label to standard deviations
        standard_deviations = ['Standard deviation'] + list(standard_deviations)

        # Add 'Relative uncertainty' label to relative uncertainties
        relative_uncertainties = ['Relative uncertainty'] + list(relative_uncertainties)

        # Create list of lists containing mean, standard deviation, and relative uncertainty data
        data = [means, standard_deviations, relative_uncertainties]

        # Concatenate original DataFrame with DataFrame containing mean, standard deviation, and relative uncertainty
        return pd.concat(objs=[df, pd.DataFrame(data=data, columns=columns)], ignore_index=True)

    def simulate(self, simulations_number):
        """
        Perform Monte Carlo simulations and return mean quantities.

        This method executes Monte Carlo simulations based on the specified number of iterations and returns
        a DataFrame containing the mean quantities, standard deviations, and relative uncertainties calculated
        from the simulation results.

        Args:
        simulations_number (int): Number of simulations to perform.

        Returns:
        pandas.DataFrame: DataFrame containing mean quantities, standard deviations, and relative uncertainties
            for the simulation results. The DataFrame includes the following columns:
            - # (int): Index column representing simulation number.
            - kVp (float): Peak kilovoltage.
            - th (float): Anode angle.
            - Air (float): Thickness of the Air filter.
            - Al (float): Thickness of the Aluminum filter.
            - Cu (float): Thickness of the Copper filter.
            - Sn (float): Thickness of the Tin filter.
            - Pb (float): Thickness of the Lead filter.
            - Be (float): Thickness of the Beryllium filter.
            - HVL1 Al (float): Half-value layer (HVL) for Aluminum.
            - HVL2 Al (float): Second HVL for Aluminum.
            - HVL1 Cu (float): HVL for Copper.
            - HVL2 Cu (float): Second HVL for Copper.
            - Mean energy (float): Mean energy of the spectrum.
            - Mean kerma (float): Mean kerma calculated using mass transmission coefficients.
            - Mean conv. coefficient. (float): Mean conversion coefficient calculated using mass transmission and
                conversion coefficients.
            Additionally, it includes rows for mean values, standard deviations, and relative uncertainties
            of the simulation results.
        """
        # Print a message indicating the start of simulation
        print('Simulation')

        # Initialize an empty list to store simulation results
        rows = []

        # For each simulation iteration
        for iteration in range(simulations_number):
            # Print a message indicating the number of the current iteration
            print(f'Iteration number: {iteration + 1}')

            # Execute a single iteration of the Monte Carlo simulation and append the results to the rows list
            rows.append([f'Iteration {iteration + 1}'] + list(self._iteration()))

        # Build results DataFrame
        return self._get_mean_quantities(rows=rows)


def random_uniform(loc, scale):
    """
    Generate a random number from a uniform distribution.

    This function generates a random number from a uniform distribution with mean `loc` and standard deviation `scale`.
    The range of the uniform distribution is calculated as `low = loc * (1 - scale * sqrt(3))`
    and `high = loc * (1 + scale * sqrt(3))`.

    Parameters:
    - loc (float): Mean of the distribution.
    - scale (float): Standard deviation of the distribution.

    Returns:
    float: A random number sampled from the uniform distribution.
    """
    # Calculate the lower bound of the uniform distribution
    low = loc * (1 - scale * np.sqrt(3))

    # Calculate the upper bound of the uniform distribution
    high = loc * (1 - scale * np.sqrt(3))

    # Generate a random number from the uniform distribution
    return float(np.random.uniform(low=low, high=high, size=1)[0])


def random_normal(loc, scale):
    """
    Generate a random number from a normal (Gaussian) distribution.

    Parameters:
    - loc (float): Mean of the distribution.
    - scale (float): Standard deviation of the distribution.

    Returns:
    float: A random number sampled from the normal distribution.

    Notes:
    This function generates a random number from a normal distribution
    with mean `loc` and standard deviation `scale * loc`.
    """
    # Generate a random number from the normal distribution
    return float(np.random.normal(loc=loc, scale=loc * scale, size=1)[0])
