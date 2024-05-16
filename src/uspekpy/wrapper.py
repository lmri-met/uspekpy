from csv import reader

import numpy as np
import pandas as pd
from scipy.interpolate import Akima1DInterpolator
from spekpy import Spek


class SpekWrapper(Spek):
    """A subclass of Spek representing a spectrum with additional methods for computation of mean quantities.

    Args:
        kvp (float): The peak kilovoltage (kVp) of the x-ray tube.
        th (float): The anode angle (deg) of the x-ray tube.

    Attributes:
        Inherits kvp and th attributes from Spek class.
    """

    def __init__(self, kvp, th):
        Spek.__init__(self, kvp, th)

    def get_mean_energy(self):
        """Compute the mean energy of an x-ray spectrum.

        This method calculates the mean energy of the spectrum using the photon fluence energy distribution and the
        energy as defined in equation of section 3.8 at ISO 4037-1:2019. It multiplies the energy values of each bin by
        their corresponding fluence values and then divides the sum of these products by the total fluence of the
        spectrum.

        Returns:
            float: The mean energy of the spectrum.
        """
        # Get the spectrum energy and photon energy fluence
        energy, fluence = self.get_spectrum(edges=False)

        # Compute the mean energy
        return sum(fluence * energy) / fluence.sum()

    def get_mean_kerma(self, mass_transfer_coefficients):
        """Compute the air kerma of an x-ray spectrum.

        This method calculates the air kerma of the x-ray spectrum using the photon fluence energy distribution and
        the mass energy transfer coefficients for air. The steps are:
        - Obtain the spectrum's energy and fluence using the `get_spectrum` method.
        - Unpack the energies and values of the mass energy transfer coefficients of air.
        - Interpolate the mass energy transfer coefficients of air for the spectrum energies in logarithmic scale using
          the `interpolate` function.
        - Computing the air kerma as defined in equation 2.16 of International Commission on Radiation Units and
          Measurements 2016 Key data for ionizing-radiation dosimetry: measurement  standards and applications ICRU
          Report 90 vol 14 Oxford University Press. This means by multiplying the fluence, energy, and interpolated mass
          energy transfer coefficients, and then dividing the sum of these products by the total fluence of the
          spectrum.

        Args:
            mass_transfer_coefficients (tuple): Tuple containing the energies and values of the mass energy transfer
                coefficients of air.

        Returns:
            float: The air kerma computed.
        """
        # Get spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)

        # Unpack the energies and values of the mass energy transfer coefficients of air
        energy_mu, mu = parse_mass_transfer_coefficients(mass_transfer_coefficients)

        # Interpolate mass energy transfer coefficients of air for the spectrum energies in logarithmic scale
        interpolated_mu = interpolate(x=energy_mu, y=mu, new_x=energy)

        # Compute air kerma
        return sum(fluence * energy * interpolated_mu) / fluence.sum()

    def get_mean_conversion_coefficient(self, mass_transfer_coefficients, conversion_coefficients, angle=None):
        """Compute the mean conversion coefficient of an x-ray spectrum.

        This method calculates the mean conversion coefficient of an x-ray spectrum using the photon fluence energy
        distribution, the mass energy transfer coefficients of air and the air kerma-to-dose-equivalent
        monoenergetic conversion coefficients. The steps are:
        - Obtain the spectrum energy and fluence using the `get_spectrum` method.
        - Unpack the energies and values of the mass energy transfer coefficients and the monoenergetic conversion
          coefficients.
        - Interpolate the mass energy transfer coefficients of air and the monoenergetic conversion coefficients for the
          spectrum energies in logarithmic scale using the `interpolate` function.
        - Compute the mean conversion coefficient. It first computes the sum of the product of the fluence, energy,
          interpolated mass energy transfer coefficients and interpolated conversion coefficients in each energy bin.
          Then, it computes the sum of the product of the fluence, energy and interpolated mass energy transfer
          coefficients in each energy bin. Finally, it divides the first sum of products by the second sum of products.

        Args:
            mass_transfer_coefficients (tuple): Tuple containing the energies and values of the mass energy transfer
                coefficients of air.
            conversion_coefficients (tuple): Tuple containing the energies and values of the monoenergetic conversion
                coefficients.
            angle (float, optional): The radiation incidence angle at which the mean conversion coefficient is
                calculated.

        Returns:
            float: The mean conversion coefficient computed.
        """
        # Get spectrum energy and fluence
        energy, fluence = self.get_spectrum(edges=False)

        # Unpack the energies and values of the mass energy transfer coefficients of air
        energy_mu, mu = parse_mass_transfer_coefficients(mass_transfer_coefficients)

        # Unpack the energies and values of the monoenergetic conversion coefficients
        energy_hk, hk = parse_conversion_coefficients(conversion_coefficients, angle)

        # Interpolate mass energy transfer coefficients of air for the spectrum energies in logarithmic scale
        interpolated_mu = interpolate(x=energy_mu, y=mu, new_x=energy)

        # Interpolate monoenergetic conversion coefficients for the spectrum energies in logarithmic scale
        interpolated_hk = interpolate(x=energy_hk, y=hk, new_x=energy)

        # Compute the mean conversion coefficient
        return sum(fluence * energy * interpolated_mu * interpolated_hk) / sum(fluence * energy * interpolated_mu)


def interpolate(x, y, new_x):
    """Interpolate y values for given new_x using Akima interpolation.

    This function performs Akima interpolation on the given x and y values to interpolate new y values for the given
    new_x. Given x, y and new_x values are assumed to be on a linear scale and are transformed to a logarithmic scale
    before the interpolation. Interpolated y values are obtained in logarithmic scale and transformed to a linear scale.
    Any resulting NaN values are replaced with zeros.

    Args:
        x (array-like): The original x values in linear scale.
        y (array-like): The original y values in linear scale.
        new_x (array-like): The new x values for interpolation in linear scale.

    Returns:
        array-like: The interpolated y values for the new_x in linear scale.
    """
    # Create an Akima1DInterpolator object with logarithmic x and y values
    interpolator = Akima1DInterpolator(x=np.log(x), y=np.log(y))

    # Interpolate new y values for given new_x using the Akima1DInterpolator
    new_y = np.exp(interpolator(x=np.log(new_x)))

    # Replace any NaN values with zeros in the interpolated y values
    return np.nan_to_num(new_y, nan=0)


def parse_mass_transfer_coefficients(coefficients):
    """Parse mass energy transfer coefficients of air into the required format by SpekWrapper.

    This function takes mass energy transfer coefficients of air in various formats and converts them into 
    the format required by SpekWrapper, which is a tuple containing two numpy arrays representing the 
    energies and mass energy transfer coefficients, respectively.

    Args:
        coefficients (tuple or str): Mass energy transfer coefficients of air. This can be either a tuple of two numpy
            arrays representing energies and mass energy transfer coefficients of air, or a string representing the file
            path of a CSV file containing two columns: energy and mass energy transfer coefficients of air.

    Returns:
        tuple: A tuple containing two numpy arrays, the energies and the mass energy transfer coefficients of air.

    Raises:
        ValueError: If the input format of mass energy transfer coefficients of air is not supported.
    """
    # Check if the input is already in the required format (tuple of two arrays)
    if is_tuple_of_two_arrays(coefficients):
        return coefficients

    # If the input is a CSV file with two columns
    elif is_csv_with_two_columns(coefficients):
        # Load CSV file into a numpy array, skipping the header
        array2d = np.genfromtxt(coefficients, delimiter=',', skip_header=1, unpack=True)
        # Build tuple of mass energy transfer coefficients
        return array2d[0], array2d[1]
    else:
        # If the input format is not supported, raise a ValueError
        raise ValueError(f"Unsupported format for mass energy transfer coefficients of air. Only a tuple of two numpy "
                         f"arrays or a CSV file with two columns are supported.")


def parse_conversion_coefficients(coefficients, irradiation_angle=None):
    """Parse conversion coefficients into the required format by SpekWrapper.

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
            column_label = next((label for label in df.columns if str(irradiation_angle) in label), None)

            # Gets the conversion coefficients from the identified column
            values = df.loc[:, column_label].values

        # Build a tuple of conversion coefficients in the required format (tuple of two numpy arrays)
        return energies, values
    else:
        # If the input format is not supported, raise a ValueError
        raise ValueError("Unsupported conversion coefficients format. Only a tuple of two numpy arrays and "
                         "a CSV file with two or more columns are supported.")


def is_tuple_of_two_arrays(arg):
    """Check if the input argument is a tuple containing two numpy arrays.

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
    """Check if a CSV file has exactly two columns.

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
