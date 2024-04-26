from csv import reader
from functools import reduce

import numpy as np
import pandas as pd


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


def parse_beam_parameters(df, column):
    # Get beam parameters in the format required by SpekWrapper (dictionary of tuples)

    # Extract values and uncertainties of filters width from the input DataFrame column
    keys = ['Al', 'Cu', 'Sn', 'Pb', 'Be', 'Air']
    values = [df.at[f'{key} filter width (mm)', column] for key in keys]
    uncertainties = [df.at[f'{key} filter width uncertainty', column] for key in keys]

    # Extract values and uncertainties of peak kilovoltage and anode angle filters width from the input DataFrame column
    # and append them to the previous lists
    keys += ['kVp', 'th']
    values += [df.at['Peak kilovoltage (kV)', column], df.at['Anode angle (deg)', column]]
    uncertainties += [df.at['Peak kilovoltage uncertainty', column],
                      df.at['Anode angle uncertainty', column]]

    # Build dictionary of beam parameters in the format required by SpekWrapper (dictionary of tuples)
    return dict(zip(keys, zip(values, uncertainties)))


def output_digest(input_df, output_dfs):
    """
    Generate a DataFrame combining input and simulation results.

    This function generates a DataFrame by combining the input DataFrame with the simulation
    results. It transforms the simulation results into a format suitable for merging with the input DataFrame
    and concatenates them accordingly.

    Args:
    input_df (pandas.DataFrame): Input DataFrame containing simulation parameters.
    output_dfs (list of pandas.DataFrame): List of DataFrames containing simulation results.

    Returns:
    pandas.DataFrame: DataFrame combining input and simulation results.

    """
    # Define columns to extract from the output DataFrames containing simulation results
    result_columns = ['HVL1 Al', 'HVL2 Al', 'HVL1 Cu', 'HVL2 Cu', 'Mean energy', 'Mean kerma',
                      'Mean conv. coefficient.']

    # Define rows to extract from the output DataFrames containing simulation results
    result_rows = ['Mean', 'Standard deviation', 'Relative uncertainty']

    # Initialize an empty list to store transformed simulation results
    results = []

    # Iterate over each output DataFrame containing simulation results
    for output_df in output_dfs:
        # Set the index of the DataFrame to '#' column
        output_df.set_index(keys='#', inplace=True)

        # Extract a subset of data from the DataFrame based on result_rows and result_columns
        df = output_df.loc[result_rows, result_columns]

        # Transpose the resulting DataFrame to swap rows and columns
        df = df.transpose()

        # Stack the DataFrame from wide to long format, creating a MultiIndex Series
        df = df.stack()

        # Convert the stacked Series back to a DataFrame
        df = df.to_frame()

        # Create a new index by combining the levels of the MultiIndex
        combined_index = df.index.map(lambda x: '{} {}'.format(x[0], x[1]))

        # Set the combined index to be the new index of the DataFrame
        df = df.set_index(combined_index)

        # Append the transformed DataFrame to the results list
        results.append(df)

    # Merge all transformed DataFrames using reduce and merge function
    merged_df = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), results)

    # Set the columns of the merged DataFrame to match the input DataFrame columns
    merged_df.columns = input_df.columns

    # Create a row to indicate the results section of the merged DataFrame
    data = [['Results'] + [None] * input_df.shape[1]]

    # Define columns for the DataFrame
    columns = ['Name'] + list(input_df.columns)

    # Create a DataFrame from data with columns specified
    df = pd.DataFrame(data, columns=columns)

    # Set the index of the DataFrame to 'Name'
    df.set_index(keys='Name', inplace=True)

    # Concatenate the input DataFrame, the results section DataFrame and the merged results DataFrame
    return pd.concat([input_df, df, merged_df])


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
