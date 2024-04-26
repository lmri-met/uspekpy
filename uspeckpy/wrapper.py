from csv import reader

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

        This method calculates the mean energy of the spectrum by multiplying
        the energy values of each bin by their corresponding fluence values
        and then dividing the sum of these products by the total fluence
        of the spectrum.

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

        This method calculates the mean kerma of the spectrum by first obtaining the spectrum's energy
        and fluence using the `get_spectrum` method. Then, it unpacks the energies and values of the mass transmission
        coefficients. Next, it interpolates the mass attenuation coefficients for the spectrum energies in logarithmic
        scale using the `interpolate` function. Finally, it computes the mean kerma by multiplying the fluence, energy,
        and interpolated mass attenuation coefficients, and then dividing the sum of these products by the total fluence
        of the spectrum.

        Parameters:
        mass_transmission_coefficients (tuple): Tuple containing the energies and values of the mass transmission
        coefficients.

        Returns:
        float: The mean kerma computed.
        """
        # Get spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)

        # Unpack the energies and values of the mass transmission coefficients
        energy_mu, mu = parse_mass_transmission_coefficients(mass_transmission_coefficients)

        # Interpolate mass attenuation coefficients for the spectrum energies in logarithmic scale
        interpolated_mu = interpolate(x=energy_mu, y=mu, new_x=energy)

        # Compute mean kerma
        return sum(fluence * energy * interpolated_mu) / fluence.sum()

    def get_mean_conversion_coefficient(self, mass_transmission_coefficients, conversion_coefficients, angle=None):
        """
        Compute the mean conversion coefficient using mass transmission coefficients and conversion coefficients.

        This method calculates the mean conversion coefficient of the spectrum by first obtaining the spectrum's energy
        and fluence using the `get_spectrum` method. Then, it unpacks the energies and values of the mass transmission
        coefficients and the conversion coefficients. Next, it interpolates the mass attenuation coefficients and the
        conversion coefficients for the spectrum energies in logarithmic scale using the `interpolate` function.
        Finally, it computes the mean conversion coefficient by multiplying the fluence, energy, interpolated mass
        attenuation coefficients, and interpolated conversion coefficients, and then dividing the sum of these products
        by the total fluence of the spectrum.

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
        energy_mu, mu = parse_mass_transmission_coefficients(mass_transmission_coefficients)

        # Unpack the energies and values of the conversion coefficients
        energy_hk, hk = parse_conversion_coefficients(conversion_coefficients, angle)

        # Interpolate mass attenuation coefficients for the spectrum energies in logarithmic scale
        interpolated_mu = interpolate(x=energy_mu, y=mu, new_x=energy)

        # Interpolate conversion coefficients for the spectrum energies in logarithmic scale
        interpolated_hk = interpolate(x=energy_hk, y=hk, new_x=energy)

        # Compute mean conversion coefficient
        return sum(fluence * energy * interpolated_mu * interpolated_hk) / fluence.sum()


def interpolate(x, y, new_x):
    """
    Interpolate y values for given new_x using Akima interpolation.

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
    # Create an Akima1DInterpolator object with logarithmic x and y values
    interpolator = Akima1DInterpolator(x=np.log(x), y=np.log(y))

    # Interpolate new y values for given new_x using the Akima1DInterpolator
    new_y = np.exp(interpolator(x=np.log(new_x)))

    # Replace any NaN values with zeros in the interpolated y values
    return np.nan_to_num(new_y, nan=0)


def parse_mass_transmission_coefficients(coefficients):
    """
    Parse mass transmission coefficients into the required format by SpekWrapper.

    This function takes mass transmission coefficients in various formats and converts them into the format required
    by SpekWrapper, which is a tuple containing two numpy arrays representing the energies and mass transmission
    coefficients, respectively.

    Args:
    coefficients (tuple or str): Mass transmission coefficients. This can be either a tuple of two numpy arrays
        representing energies and mass transmission coefficients, or a string representing the file path of a CSV
        file containing two columns: energy and mass transmission coefficients.

    Returns:
    tuple: A tuple containing two numpy arrays representing energies and mass transmission coefficients, respectively.

    Raises:
    ValueError: If the input format of mass transmission coefficients is not supported.
    """
    # Check if the input is already in the required format (tuple of two arrays)
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
        raise ValueError(f"Unsupported mass transmission coefficients format. Only a tuple of two numpy arrays and "
                         f"a CSV file with two columns are supported.")


def parse_conversion_coefficients(coefficients, irradiation_angle):
    """
    Parse conversion coefficients into the required format by SpekWrapper.

    This function takes conversion coefficients in various formats and converts them into the format required
    by SpekWrapper, which is a tuple containing two numpy arrays representing the energies and conversion
    coefficients, respectively.

    Args:
    coefficients (tuple or str): Conversion coefficients. This can be either a tuple of two numpy arrays
        representing energies and conversion coefficients, or a string representing the file path of a CSV
        file containing at least two columns: energy and conversion coefficients.
    irradiation_angle (str): Irradiation angle for which the conversion coefficients are calculated.

    Returns:
    tuple: A tuple containing two numpy arrays representing energies and conversion coefficients, respectively.

    Raises:
    ValueError: If the input format of conversion coefficients is not supported.

    """
    # Check if the input is already in the required format (tuple of two arrays)
    if is_tuple_of_two_arrays(coefficients):
        return coefficients
    # If the input is a valid CSV file
    elif is_valid_csv(coefficients):
        # Read CSV file into a DataFrame
        df = pd.read_csv(coefficients)

        # Get the energies from the first column of the DataFrame
        energies = df.iloc[:, 0].values

        # If the DataFrame has only 2 columns, the second column contains the conversion coefficients
        if df.shape[1] == 2:
            values = df.iloc[:, 1].values
        else:
            # Find the column containing the conversion coefficients corresponding to the specified irradiation angle
            column_label = next((label for label in df.columns if irradiation_angle in label), None)

            # Get the conversion coefficients from the identified column
            values = df.loc[:, column_label].values

        # Build a tuple of conversion coefficients in the required format (tuple of two numpy arrays)
        return energies, values
    else:
        # If the input format is not supported, raise a ValueError
        raise ValueError("Unsupported conversion coefficients format. Only a tuple of two numpy arrays and "
                         "a CSV file with two or more columns are supported.")


def is_tuple_of_two_arrays(arg):
    """
    Check if the input argument is a tuple containing two numpy arrays.

    This function verifies if the input argument is a tuple containing exactly two numpy arrays.
    If the input argument meets the criteria, the function returns True; otherwise, it returns False.

    Args:
    arg (any): The input argument to be validated.

    Returns:
    bool: True if the input argument is a tuple containing two numpy arrays, False otherwise.
    """
    # Check if the input is a tuple
    if not isinstance(arg, tuple):
        return False

    # Check if the tuple contains exactly two elements
    if len(arg) != 2:
        return False

    # Check if both elements of the tuple are numpy arrays
    if not isinstance(arg[0], np.ndarray) or not isinstance(arg[1], np.ndarray):
        return False

    # If all conditions are met, return True
    return True


def is_csv_with_two_columns(filename):
    """
    Check if a CSV file has exactly two columns.

    This function reads the provided CSV file and checks if each row contains exactly two columns.
    If any row has a different number of columns, the function returns False.
    If all rows have exactly two columns, the function returns True.

    Args:
    filename (str): The path to the CSV file to be validated.

    Returns:
    bool: True if the CSV file has exactly two columns in each row, False otherwise.
    """
    # Check if the file has a CSV extension
    if not filename.lower().endswith('.csv'):
        return False

    # Open the CSV file
    with open(filename, 'r') as file:
        rows = reader(file)

        # Iterate over each row in the CSV file
        for row in rows:
            # If any row does not contain exactly two columns, return False
            if len(row) != 2:
                return False
        # If all rows contain exactly two columns, return True
        return True


def is_valid_csv(filename):
    """
    Check if a CSV file is valid by ensuring that all rows have the same length.

    This function reads the provided CSV file and checks if all rows have the same length.
    If any row has a different length from the first row, the function returns False.
    If all rows have the same length, the function returns True, indicating that the CSV file is valid.

    Args:
    filename (str): The path to the CSV file to be validated.

    Returns:
    bool: True if the CSV file is valid (all rows have the same length), False otherwise.
    """
    # Check if the file has a CSV extension
    if not filename.lower().endswith('.csv'):
        return False

    # Open the CSV file
    with open(filename, 'r') as file:
        rows = reader(file)

        # Get the length of the first row
        first_row_length = len(next(rows))

        # Check the length of each subsequent row
        for row in rows:
            # If any row has a different length, return False
            if len(row) != first_row_length:
                return False
        # If all rows have the same length, return True
        return True
