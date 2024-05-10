import os

import pytest
import numpy as np
from scipy.interpolate import Akima1DInterpolator

import src.uspekpy.wrapper as wrp


class TestIsTupleOfTwoArrays:
    def test_valid_input(self):
        arg = (np.array([1, 2, 3]), np.array([4, 5, 6]))
        assert wrp.is_tuple_of_two_arrays(arg) is True

    def test_invalid_input_not_tuple(self):
        arg = np.array([1, 2, 3])
        assert wrp.is_tuple_of_two_arrays(arg) is False

    def test_invalid_input_wrong_length(self):
        arg = (np.array([1, 2, 3]), np.array([4, 5, 6]), np.array([7, 8, 9]))
        assert wrp.is_tuple_of_two_arrays(arg) is False

    def test_invalid_input_not_arrays(self):
        arg = (np.array([1, 2, 3]), [4, 5, 6])
        assert wrp.is_tuple_of_two_arrays(arg) is False


class TestIsCsvWithTwoColumns:
    def test_valid_csv_with_two_columns(self, tmpdir):
        # Create a temporary CSV file with two columns
        filename = os.path.join(tmpdir, 'test.csv')
        with open(filename, 'w') as file:
            file.write('Column1,Column2\n')
            file.write('Value1,Value2\n')

        assert wrp.is_csv_with_two_columns(filename) is True

    def test_invalid_non_csv_file(self, tmpdir):
        # Create a temporary text file
        filename = os.path.join(tmpdir, 'test.txt')
        with open(filename, 'w') as file:
            file.write('This is a test file.')

        assert wrp.is_csv_with_two_columns(filename) is False

    def test_invalid_csv_with_one_column(self, tmpdir):
        # Create a temporary CSV file with one column
        filename = os.path.join(tmpdir, 'test.csv')
        with open(filename, 'w') as file:
            file.write('Column1\n')
            file.write('Value1\n')

        assert wrp.is_csv_with_two_columns(filename) is False


class TestIsValidCsv:
    def test_valid_csv(self, tmpdir):
        # Create a temporary CSV file with valid rows
        filename = os.path.join(tmpdir, 'test.csv')
        with open(filename, 'w') as file:
            file.write('Column1,Column2\n')
            file.write('Value1,Value2\n')
            file.write('Value3,Value4\n')

        assert wrp.is_valid_csv(filename) is True

    def test_invalid_non_csv_file(self, tmpdir):
        # Create a temporary text file
        filename = os.path.join(tmpdir, 'test.txt')
        with open(filename, 'w') as file:
            file.write('This is a test file.')

        assert wrp.is_valid_csv(filename) is False

    def test_invalid_csv_different_row_lengths(self, tmpdir):
        # Create a temporary CSV file with rows of different lengths
        filename = os.path.join(tmpdir, 'test.csv')
        with open(filename, 'w') as file:
            file.write('Column1,Column2\n')
            file.write('Value1,Value2\n')
            file.write('Value3\n')

        assert wrp.is_valid_csv(filename) is False


class TestParseMassTransmissionCoefficients:
    @pytest.fixture
    def coefficients_tuple(self):
        energies = np.array([1, 2, 3])
        coefficients = np.array([0.1, 0.2, 0.3])
        return energies, coefficients

    def test_parse_tuple_input(self, coefficients_tuple):
        result = wrp.parse_mass_transmission_coefficients(coefficients_tuple)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], np.ndarray)
        assert isinstance(result[1], np.ndarray)

    def test_parse_csv_input(self, tmpdir):
        # Create a temporary CSV file with two columns
        filename = tmpdir.join("test.csv")
        with open(filename, "w") as file:
            file.write("energy,mass_transmission\n")
            file.write("1,0.1\n")
            file.write("2,0.2\n")
            file.write("3,0.3\n")

        result = wrp.parse_mass_transmission_coefficients(str(filename))
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], np.ndarray)
        assert isinstance(result[1], np.ndarray)

    def test_invalid_input(self):
        with pytest.raises(ValueError):
            wrp.parse_mass_transmission_coefficients("invalid_input")


class TestParseConversionCoefficients:
    @pytest.fixture
    def coefficients_tuple(self):
        energies = np.array([1, 2, 3])
        coefficients = np.array([0.1, 0.2, 0.3])
        return energies, coefficients

    def test_parse_tuple_input(self, coefficients_tuple):
        result = wrp.parse_conversion_coefficients(coefficients_tuple, None)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], np.ndarray)
        assert isinstance(result[1], np.ndarray)

    def test_parse_csv_input_two_columns(self, tmpdir):
        # Create a temporary CSV file with two columns
        filename = tmpdir.join("test.csv")
        with open(filename, "w") as file:
            file.write("energy,angle20\n")
            file.write("1,0.2\n")
            file.write("2,0.4\n")
            file.write("3,0.6\n")

        result = wrp.parse_conversion_coefficients(str(filename), None)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], np.ndarray)
        assert isinstance(result[1], np.ndarray)

    def test_parse_csv_input_more_than_two_columns(self, tmpdir):
        # Create a temporary CSV file with more than two columns
        filename = tmpdir.join("test.csv")
        with open(filename, "w") as file:
            file.write("energy,angle10,angle20\n")
            file.write("1,0.1,0.2\n")
            file.write("2,0.3,0.4\n")
            file.write("3,0.5,0.6\n")

        result = wrp.parse_conversion_coefficients(str(filename), 20)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], np.ndarray)
        assert isinstance(result[1], np.ndarray)

    def test_invalid_input(self):
        with pytest.raises(ValueError):
            wrp.parse_conversion_coefficients("invalid_input", "angle")


class TestInterpolate:
    @staticmethod
    def test_interpolate():
        # Define original x and y values
        x = [1, 3, 5]
        y = [10, 30, 50]

        # Define new x values for interpolation
        new_x = [2, 4]

        # Call the function with the original and new x values
        new_y = wrp.interpolate(x, y, new_x)

        # Define expected interpolated y values
        expected_new_y = np.array([20.0, 40.0])

        # Assert that the result is an array-like object and contains the expected interpolated y values
        assert isinstance(new_y, np.ndarray)
        np.testing.assert_allclose(new_y, expected_new_y, rtol=1e-6)


class TestSpekWrapper:

    @pytest.fixture
    def spectrum_energy_fluence(self):
        # Define x-ray beam parameters
        my_filters = [
            ('Al', 4),
            ('Cu', 0.6),
            ('Sn', 0),
            ('Pb', 0),
            ('Be', 0),
            ('Air', 1000)
        ]

        # Initialises an SpeckWrapper object and add filters
        spectrum = wrp.SpekWrapper(kvp=60, th=20)
        spectrum.multi_filter(my_filters)

        # Get spectrum energy and fluence
        energy, fluence = spectrum.get_spectrum(edges=False)
        return spectrum, energy, fluence

    @pytest.fixture
    def mass_transmission_coefficients(self):
        # Define the energies and values of the mass energy transfer coefficients for air
        energy_mu, mu = (
            np.array(
                [1.0, 1.1726, 1.25, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0, 3.2063, 3.206301, 3.22391, 3.25051, 3.5, 3.61881,
                 4.0,
                 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 14.0, 15.0, 17.5, 20.0, 25.0, 28.6633, 30.0, 35.0, 40.0, 50.0,
                 60.0,
                 70.0, 80.0, 90.0, 100.0, 125.0, 140.0, 150.0, 175.0, 187.083, 200.0, 250.0, 300.0, 324.037, 350.0,
                 386.867,
                 400.0, 474.342, 500.0, 574.456, 600.0, 673.537, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0, 1558.93,
                 1750.0, 1870.83, 2000.0, 2345.21, 2500.0, 3000.0, 3240.37, 3500.0, 4000.0, 4500.0, 5000.0, 6000.0,
                 6480.74,
                 7000.0, 8000.0, 9000.0]),
            np.array(
                [3487.7, 2271.66, 1907.85, 1396.25, 1152.41, 746.85, 510.495, 267.712, 156.677, 128.597, 139.322,
                 139.22,
                 136.749, 110.244, 99.9467, 74.3863, 38.3165, 22.1387, 13.8638, 9.20954762864289, 6.40746485225397,
                 4.62692014996195, 2.30743153037744, 1.61505075253763, 1.30130979037554, 0.801145676760343,
                 0.525990443101652, 0.262079759969996, 0.17213944380102, 0.150643910752649, 0.096166283094939,
                 0.067006659158694, 0.040350724533314, 0.030060380204255, 0.025736225862943, 0.023919178948042,
                 0.023273564441493, 0.023169137685326, 0.023905175732908, 0.024501233602389, 0.024926643691358,
                 0.025877217107796, 0.026274062287156, 0.026655131026205, 0.027882238669658, 0.028683812435718,
                 0.028972447929017, 0.029148771809045, 0.029415131704691, 0.029490249846129, 0.029698178702551,
                 0.029685041261757, 0.029603764611649, 0.029549718949812, 0.029405101930449, 0.029197657857256,
                 0.028890614121642, 0.028390536220701, 0.027936798166738, 0.026719405747041, 0.02562121786982,
                 0.025359438553807, 0.02463108715039, 0.024240842040503, 0.023753744772526, 0.022671649894621,
                 0.022309768351513, 0.021132688079239, 0.020639133403572, 0.020121558259145, 0.019347638398287,
                 0.018700741656237, 0.018185053120065, 0.017326066208925, 0.016966042150857, 0.016690444038446,
                 0.016199645451199, 0.015809252632524])
        )
        return energy_mu, mu

    @pytest.fixture
    def conversion_coefficients(self):
        # Define the energies and values of the monoenergetic conversion coefficients
        energy_hk, hk = (
            np.array(
                [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 500,
                 600, 800, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000]),
            np.array(
                [1.2e-05, 9.5e-05, 0.00145, 0.008, 0.0331, 0.0737, 0.127, 0.19, 0.26, 0.326, 0.395, 0.466, 0.538, 0.61,
                 1.1, 1.47, 1.67, 1.74, 1.72, 1.65, 1.49, 1.4, 1.31, 1.26, 1.23, 1.21, 1.19, 1.17, 1.15, 1.13, 1.12,
                 1.11, 1.19, 1.09, 1.08, 1.06])
        )
        return energy_hk, hk

    def test_get_mean_energy(self, spectrum_energy_fluence):
        # Get spectrum and spectrum energy and fluence
        spectrum, energy, fluence = spectrum_energy_fluence

        # Compute expected mean energy
        expected_mean_energy = sum(fluence * energy) / fluence.sum()

        # Compute mean energy with SpekWrapper.get_mean_energy()
        mean_energy = spectrum.get_mean_energy()

        assert mean_energy == expected_mean_energy

    def test_get_mean_kerma(self, spectrum_energy_fluence, mass_transmission_coefficients):
        # Get spectrum and spectrum energy and fluence
        spectrum, energy, fluence = spectrum_energy_fluence

        # Define the energies and values of the mass energy transfer coefficients for air
        energy_mu, mu = mass_transmission_coefficients

        # Create an Akima1DInterpolator with logarithmic energies and mass energy transfer coefficients for air
        interpolator = Akima1DInterpolator(x=np.log(energy_mu), y=np.log(mu))

        # Interpolate mass energy transfer coefficients for air for the spectrum energies in logarithmic scale
        interpolated_mu = np.exp(interpolator(x=np.log(energy)))

        # Replace any NaN values with zeros in the interpolated mass energy transfer coefficients for air
        interpolated_mu = np.nan_to_num(interpolated_mu, nan=0)

        # Compute mean air kerma
        expected_mean_kerma = sum(fluence * energy * interpolated_mu) / fluence.sum()

        # Compute mean air kerma with SpekWrapper.get_mean_kerma()
        mean_kerma = spectrum.get_mean_kerma(mass_transmission_coefficients=(energy_mu, mu))

        assert mean_kerma == expected_mean_kerma

    def test_get_mean_conversion_coefficient(self, spectrum_energy_fluence, mass_transmission_coefficients,
                                             conversion_coefficients):
        # Get spectrum and spectrum energy and fluence
        spectrum, energy, fluence = spectrum_energy_fluence

        # Define the energies and values of the mass energy transfer coefficients for air
        energy_mu, mu = mass_transmission_coefficients

        # Define the energies and values of the monoenergetic conversion coefficients
        energy_hk, hk = conversion_coefficients

        # Create an Akima1DInterpolator with logarithmic energies and mass energy transfer coefficients for air
        interpolator = Akima1DInterpolator(x=np.log(energy_mu), y=np.log(mu))

        # Interpolate mass energy transfer coefficients for air for the spectrum energies in logarithmic scale
        interpolated_mu = np.exp(interpolator(x=np.log(energy)))

        # Replace any NaN values with zeros in the interpolated mass energy transfer coefficients for air
        interpolated_mu = np.nan_to_num(interpolated_mu, nan=0)

        # Create an Akima1DInterpolator with logarithmic energies and monoenergetic conversion coefficients
        interpolator = Akima1DInterpolator(x=np.log(energy_hk), y=np.log(hk))

        # Interpolate monoenergetic conversion coefficients for the spectrum energies in logarithmic scale
        interpolated_hk = np.exp(interpolator(x=np.log(energy)))

        # Replace any NaN values with zeros in the interpolated monoenergetic conversion coefficients
        interpolated_hk = np.nan_to_num(interpolated_hk, nan=0)

        # Compute mean conversion coefficient
        expected_hk = sum(fluence * energy * interpolated_mu * interpolated_hk) / sum(fluence * energy * interpolated_mu)

        # Compute mean conversion coefficient with SpekWrapper.get_mean_conversion_coefficient()
        mean_hk = spectrum.get_mean_conversion_coefficient(mass_transmission_coefficients=(energy_mu, mu),
                                                           conversion_coefficients=(energy_hk, hk))

        assert mean_hk == expected_hk
