from csv import reader

import numpy as np
import pandas as pd
from scipy.interpolate import Akima1DInterpolator
from spekpy import Spek


class SpekWrapper(Spek):
    """
    A subclass of Spek representing a spectrum with additional methods for computation of mean quantities.

    Parameters:
    - kvp (float): The peak kilovoltage (kVp) of the x-ray tube.
    - th (float): The anode angle of the x-ray tube.

    Attributes:
    Inherits kvp and th attributes from Spek class.

    Methods:
    - get_mean_energy(): Computes the mean energy of the spectrum.
    - get_mean_kerma(mass_transmission_coefficients): Computes the air kerma of the spectrum using mass energy transfer
      coefficients of air and the photon energy fluence.
    - get_mean_conversion_coefficient(mass_transmission_coefficients, conversion_coefficients, angle=None):
      Computes the mean conversion coefficient of the spectrum using the photon energy fluence, mass energy transfer  
      coefficients of air and air kerma-to-dose-equivalent monoenergetic conversion coefficients (referred later as 
      "monoenergetic conversion coefficients" see ISO 4037-1).
    """

    def __init__(self, kvp, th):
        """
        Initialises SpecWrapper instance with peak kilovoltage and the anode angle of the x-ray tube.
        """
        Spek.__init__(self, kvp, th)

    def get_mean_energy(self):
        """
        Computes the mean energy of the spectrum.

        This method calculates the mean energy of the spectrum (see equation of section 3.8 at ISO 4037-1:2019) 
        using the photon fluence and the energy value in each bin.

        Returns:
        float: The mean energy of the spectrum.
        """
        # Gets the spectrum energy and photon energy fluence
        energy, fluence = self.get_spectrum(edges=False)

        # Computes the mean energy
        return sum(fluence * energy) / fluence.sum()

    def get_mean_kerma(self, mass_transmission_coefficients):
        """
        Computes the air kerma using the distribution of the photon fluence with respect to energy and the 
        mass energy transfer coefficients for air.

        This method calculates the kerma in air of the photon spectrum by:
        i)   first obtaining the spectrum's energy and fluence using the `get_spectrum` method,
        ii)  then unpacking the energies and values of the mass energy transfer coefficients of air,
        iii) interpolating the mass energy transfer coefficients of air for the spectrum energies 
             in logarithmic scale using the `interpolate` function and 
        iv)  finally, computing the air kerma (as defined in equation 2.16 of International Commission 
             on Radiation Units and Measurements 2016 Key data for ionizing-radiation dosimetry: measurement 
             standards and applications ICRU Report 90 vol 14 Oxford University Press).

        Parameters:
        mass_transmission_coefficients (tuple): Tuple containing the energies and values of the mass energy transfer
        coefficients of air.

        Returns:
        float: The air kerma computed.
        """
        # Gets spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)

        # Unpacks the energies and values of the mass energy transfer coefficients of air
        energy_mu, mu = parse_mass_transmission_coefficients(mass_transmission_coefficients)

        # Interpolates mass energy transfer coefficients of air for the spectrum energies in logarithmic scale
        interpolated_mu = interpolate(x=energy_mu, y=mu, new_x=energy)

        # Computes air kerma
        return sum(fluence * energy * interpolated_mu) / fluence.sum()

    def get_mean_conversion_coefficient(self, mass_transmission_coefficients, conversion_coefficients, angle=None):
        """
        Computes the mean conversion coefficient using mass energy transfer coefficients of air and monoenergetic 
        conversion coefficients.

        This method calculates the mean conversion coefficient of the spectrum by:
        i)   first obtaining the spectrum energy and fluence using the `get_spectrum` method,
        ii)  then, unpacking the energies and values of the mass energy transfer coefficients 
             and the monoenergetic conversion coefficients,
        iii) interpolating the mass enery transfer coefficients of air and the monoenergetic 
             conversion coefficients for the spectrum energies in logarithmic scale using the `interpolate` function and
        iv)  finally, computing the mean conversion coefficient.

        Parameters:
        mass_transmission_coefficients (tuple): Tuple containing the energies and values of the mass energy
        transfer coefficients of air.
        conversion_coefficients (tuple): Tuple containing the energies and values of the monoenergetic conversion 
        coefficients. 
        angle (float, optional): The angle of radiation incidence at which the mean conversion coefficient is calculated.

        Returns:
        float: The mean conversion coefficient computed.
        """
        # Gets spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)

        # Unpacks the energies and values of the mass energy transfer coefficients of air
        energy_mu, mu = parse_mass_transmission_coefficients(mass_transmission_coefficients)

        # Unpacks the energies and values of the monoenergetic conversion coefficients
        energy_hk, hk = parse_conversion_coefficients(conversion_coefficients, angle)

        # Interpolates mass energy transfer coefficients of air for the spectrum energies in logarithmic scale
        interpolated_mu = interpolate(x=energy_mu, y=mu, new_x=energy)

        # Interpolates monoenergetic conversion coefficients for the spectrum energies in logarithmic scale
        interpolated_hk = interpolate(x=energy_hk, y=hk, new_x=energy)

        # Computes the mean conversion coefficient
        return sum(fluence * energy * interpolated_mu * interpolated_hk) / sum(fluence * energy * interpolated_mu)


def interpolate(x, y, new_x):
    """
    Interpolates y values for given new_x using Akima interpolation.

    This function performs Akima interpolation on the given x and y values
    (assumed to be on a logarithmic scale) to interpolate new y values for
    the given new_x. Any resulting NaN values are replaced with zeros.

    Parameters:
    - x (array-like): The original x values.
    - y (array-like): The original y values.
    - new_x (array-like): The new x values for interpolation.

    Returns:
    array-like: The interpolated y values for the new_x.
    """
    # Creates an Akima1DInterpolator object with logarithmic x and y values
    interpolator = Akima1DInterpolator(x=np.log(x), y=np.log(y))

    # Interpolates new y values for given new_x using the Akima1DInterpolator
    new_y = np.exp(interpolator(x=np.log(new_x)))

    # Replaces any NaN values with zeros in the interpolated y values
    return np.nan_to_num(new_y, nan=0)


def parse_mass_transmission_coefficients(coefficients):
    """
    Parse mass energy transfer coefficients of air into the required format by SpekWrapper.

    This function takes mass energy transfer coefficients of air in various formats and converts them into 
    the format required by SpekWrapper, which is a tuple containing two numpy arrays representing the 
    energies and mass transmission coefficients, respectively.

    Args:
    coefficients (tuple or str): Mass energy transfer coefficients of air. This can be either a tuple of two numpy arrays
        representing energies and mass energy transfer coefficients of air, or a string representing the file path of a CSV
        file containing two columns: energy and mass energy transfer coefficients of air.

    Returns:
    tuple: A tuple containing two numpy arrays representing energies and mass energy transfer coefficients of air, respectively.

    Raises:
    ValueError: If the input format of mass energy transfer coefficients of air is not supported.
    """
    # Checks if the input is already in the required format (tuple of two arrays)
    if is_tuple_of_two_arrays(coefficients):
        return coefficients

    # If the input is a CSV file with two columns
    elif is_csv_with_two_columns(coefficients):
        # Load CSV file into a numpy array, skipping the header
        array2d = np.genfromtxt(coefficients, delimiter=',', skip_header=1, unpack=True)
        # Build tuple of mass transmission coefficients
        return array2d[0], array2d[1]
    else:
        # If the input format is not supported, raise a ValueError
        raise ValueError(f"Unsupported mass energy transfer coefficients of air format. Only a tuple of two numpy arrays and "
                         f"a CSV file with two columns are supported.")


def parse_conversion_coefficients(coefficients, irradiation_angle=None):
    """
    Parse conversion coefficients into the required format by SpekWrapper.

    This function takes conversion coefficients in various formats and converts them into the format required
    by SpekWrapper, which is a tuple containing two numpy arrays representing the energies and conversion
    coefficients, respectively.

    Args:
    coefficients (tuple or str): Conversion coefficients. This can be either a tuple of two numpy arrays
        representing energies and conversion coefficients, or a string representing the file path of a CSV
        file containing at least two columns: energy and conversion coefficients.
    irradiation_angle (int, optional): Irradiation angle for which the conversion coefficients are calculated.
        Defaults to None.

    Returns:
    tuple: A tuple containing two numpy arrays representing energies and conversion coefficients, respectively.

    Raises:
    ValueError: If the input format of conversion coefficients is not supported.

    """
    # Checks if the input is already in the required format (tuple of two arrays)
    if is_tuple_of_two_arrays(coefficients):
        return coefficients
    # If the input is a valid CSV file
    elif is_valid_csv(coefficients):
        # Reads CSV file into a DataFrame
        df = pd.read_csv(coefficients)

        # Gets the energies from the first column of the DataFrame
        energies = df.iloc[:, 0].values

        # If the DataFrame has only 2 columns, the second column contains the conversion coefficients
        if df.shape[1] == 2:
            values = df.iloc[:, 1].values
        else:
            # Finds the column containing the conversion coefficients corresponding to the specified irradiation angle
            column_label = next((label for label in df.columns if str(irradiation_angle) in label), None)

            # Gets the conversion coefficients from the identified column
            values = df.loc[:, column_label].values

        # Builds a tuple of conversion coefficients in the required format (tuple of two numpy arrays)
        return energies, values
    else:
        # If the input format is not supported, raise a ValueError
        raise ValueError("Unsupported conversion coefficients format. Only a tuple of two numpy arrays and "
                         "a CSV file with two or more columns are supported.")


def is_tuple_of_two_arrays(arg):
    """
    Checks if the input argument is a tuple containing two numpy arrays.

    This function verifies if the input argument is a tuple containing exactly two numpy arrays.
    If the input argument meets the criteria, the function returns True; otherwise, it returns False.

    Args:
    arg (any): The input argument to be validated.

    Returns:
    bool: True if the input argument is a tuple containing two numpy arrays, False otherwise.
    """
    # Checks if the input is a tuple
    if not isinstance(arg, tuple):
        return False

    # Checks if the tuple contains exactly two elements
    if len(arg) != 2:
        return False

    # Checks if both elements of the tuple are numpy arrays
    if not isinstance(arg[0], np.ndarray) or not isinstance(arg[1], np.ndarray):
        return False

    # If all conditions are met, return True
    return True


def is_csv_with_two_columns(filename):
    """
    Checks if a CSV file has exactly two columns.

    This function reads the provided CSV file and checks if each row contains exactly two columns.
    If any row has a different number of columns, the function returns False.
    If all rows have exactly two columns, the function returns True.

    Args:
    filename (str): The path to the CSV file to be validated.

    Returns:
    bool: True if the CSV file has exactly two columns in each row, False otherwise.
    """
    # Checks if the file has a CSV extension
    if not filename.lower().endswith('.csv'):
        return False

    # Opens the CSV file
    with open(filename, 'r') as file:
        rows = reader(file)

        # Iterates over each row in the CSV file
        for row in rows:
            # If any row does not contain exactly two columns, return False
            if len(row) != 2:
                return False
        # If all rows contain exactly two columns, return True
        return True


def is_valid_csv(filename):
    """
    Checks if a CSV file is valid by ensuring that all rows have the same length.

    This function reads the provided CSV file and checks if all rows have the same length.
    If any row has a different length from the first row, the function returns False.
    If all rows have the same length, the function returns True, indicating that the CSV file is valid.

    Args:
    filename (str): The path to the CSV file to be validated.

    Returns:
    bool: True if the CSV file is valid (all rows have the same length), False otherwise.
    """
    # Checks if the file has a CSV extension
    if not filename.lower().endswith('.csv'):
        return False

    # Opens the CSV file
    with open(filename, 'r') as file:
        rows = reader(file)

        # Gets the length of the first row
        first_row_length = len(next(rows))

        # Checks the length of each subsequent row
        for row in rows:
            # If any row has a different length, return False
            if len(row) != first_row_length:
                return False
        # If all rows have the same length, return True
        return True
