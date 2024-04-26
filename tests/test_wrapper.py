import os

import pytest
import numpy as np

import uspeckpy.wrapper as wrp


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

# TODO: interpolate(), SpekWrapper()
