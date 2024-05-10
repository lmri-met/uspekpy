from functools import reduce

import pandas as pd

from src.uspekpy.uspek import USpek
from src.uspekpy.wrapper import parse_mass_transmission_coefficients, parse_conversion_coefficients


def batch_simulation(input_file_path, sheet_name=None):
    """
    Perform batch simulations and generate output data.

    This function performs batch simulations based on input data provided in an Excel or CSV file.
    It reads the input data into a DataFrame, processes each simulation, and generates output data.
    The output data is stored in a list of DataFrames, one for each simulation. Finally, it generates
    a DataFrame by combining the input DataFrame with the simulation results DataFrames.

    Args:
    input_file_path (str): The path to the input file containing simulation parameters.
    sheet_name (str, optional): The name of the sheet to read if the input file is in Excel format.

    Returns:
    list: DataFrame combining input and simulation results.

    Raises:
    ValueError: If the input file format is not supported or if there are issues with the input data.
    """
    # Print a message indicating the start of input digestion
    print('Batch simulation')

    # Print a message indicating the start of initial input digestion
    print('\nInitial input digest')

    # Read input Excel or CSV file into a DataFrame
    input_df = read_file_to_dataframe(input_file_path, sheet_name=sheet_name)

    # Read Excel file into a DataFrame and set 'Name' column as index
    input_df.set_index(keys='Name', inplace=True)

    # Initialize an empty list to store the simulations results
    output_dfs = []

    # Initialize an iterator
    i = 0

    # Iterate over each column in the input DataFrame
    for column_name in input_df.columns:
        # Print a message indicating simulation number
        print(f'\nSimulation {i + 1}')

        # Print a message indicating the start of input digestion
        print('Input digest')

        # Parse beam parameters from the input DataFrame
        beam_parameters = parse_beam_parameters(df=input_df, column=column_name)

        # Extract mass energy transfer coefficients file path from the input DataFrame
        file_path = input_df.at['mass energy transfer coefficients file', column_name]

        # Parse mass energy transfer coefficients from the file
        mass_transmission_coefficients = parse_mass_transmission_coefficients(coefficients=file_path)

        # Extract mono-energetic conversion coefficients file path from the input DataFrame
        file_path = input_df.at['Mono-energetic conversion coefficients file', column_name]

        # Extract irradiation angle from the input DataFrame
        irradiation_angle = input_df.at['Irradiation angle (deg)', column_name]

        # Parse conversion coefficients from the file
        conversion_coefficients = parse_conversion_coefficients(coefficients=file_path,
                                                                irradiation_angle=irradiation_angle)

        # Extract number of simulations from the input DataFrame
        simulations_number = input_df.at['Number of simulations', column_name]

        # Extract mass energy transfer coefficients uncertainty from the input DataFrame
        mass_transmission_coefficients_uncertainty = input_df.at[
            'mass energy transfer coefficients uncertainty', column_name]

        # Print a message indicating the start of input digestion
        print('Simulation')

        # Create a USpek object with the specified parameters
        s = USpek(beam_parameters=beam_parameters, mass_transmission_coefficients=mass_transmission_coefficients,
                  mass_transmission_coefficients_uncertainty=mass_transmission_coefficients_uncertainty,
                  conversion_coefficients=conversion_coefficients)

        # Run simulation with the specified number of iterations
        output_df = s.simulate(simulations_number=simulations_number)

        # Print a message indicating the start of input digestion
        print('Output digest')

        # Append output DataFrame to output DataFrames list
        output_dfs.append(output_df)

        # Increment the simulation iterator
        i += 1

    # Print a message indicating the start of input digestion
    print('\nFinal output digest')

    # Return DataFrame with the simulation results
    return output_digest(input_df=input_df, output_dfs=output_dfs)


def read_file_to_dataframe(file_path, sheet_name=None):
    """
    Reads data from a file into a pandas DataFrame.

    This function reads data from a file specified by the file_path parameter into a pandas DataFrame.
    The file can be either in Excel (.xlsx, .xls) or CSV (.csv) format. If the file is in Excel format,
    the optional sheet_name parameter can be used to specify the sheet name to read. If the file is in CSV
    format, the function automatically detects the delimiter.

    Args:
    file_path (str): The path to the file to be read.
    sheet_name (str, optional): The name of the sheet to read if the file is in Excel format.

    Returns:
    pandas.DataFrame: A DataFrame containing the data read from the file.

    Raises:
    ValueError: If the file format is not supported. Only Excel (.xlsx, .xls) and CSV (.csv) files are supported.
    """
    # Check the file extension to determine the file type
    if file_path.endswith('.xlsx') or file_path.endswith('.xls'):
        # Read Excel file
        df = pd.read_excel(file_path, sheet_name=sheet_name)
    elif file_path.endswith('.csv'):
        # Read CSV file
        df = pd.read_csv(file_path)
        # Iterate over all elements in the DataFrame
        for column in df.columns:
            for idx, value in enumerate(df[column]):
                try:
                    # Try to convert the value to an integer
                    df.at[idx, column] = int(value)
                except ValueError:
                    try:
                        # If conversion to int fails, try converting to float
                        df.at[idx, column] = float(value)
                    except ValueError:
                        # If conversion to float also fails, leave it unchanged
                        pass
    else:
        # Raise a ValueError for unsupported file formats
        raise ValueError("Unsupported file format. Only Excel (.xlsx, .xls) and CSV (.csv) files are supported.")

    return df


def parse_beam_parameters(df, column):
    """
    Gets beam parameters in the format required by SpekWrapper.

    The DataFrame should have rows in the specified column named:
    - '<filter name> filter width (mm)' for filter widths.
    - '<filter name> filter width uncertainty' for filter width uncertainties.
    - 'Peak kilovoltage (kV)' for peak kilovoltage.
    - 'Anode angle (deg)' for anode angle.
    - 'Peak kilovoltage uncertainty' for peak kilovoltage uncertainty.
    - 'Anode angle uncertainty' for anode angle uncertainty.

    Parameters:
        df (pd.DataFrame): Input DataFrame containing beam parameters.
        column (str): Column name in the DataFrame containing the beam parameters.

    Returns:
        Dict[str, Tuple[float, float]]: Dictionary of beam parameters in the format required by SpekWrapper,
        where keys are parameter names and values are tuples of parameter values and uncertainties.
    """
    # Keys for the filters
    keys = ['Al', 'Cu', 'Sn', 'Pb', 'Be', 'Air']

    # Extracting values for each filter from the specified DataFrame column
    values = [df.at[f'{key} filter width (mm)', column] for key in keys]

    # Extracting uncertainties for each filter from the specified DataFrame column
    uncertainties = [df.at[f'{key} filter width uncertainty', column] for key in keys]

    # Appends additional keys for peak kilovoltage and anode angle
    keys += ['kVp', 'th']

    # Extracts values for peak kilovoltage and anode angle from the specified DataFrame column and appends them
    values += [df.at['Peak kilovoltage (kV)', column], df.at['Anode angle (deg)', column]]

    # Extracts uncertainties for peak kilovoltage and anode angle from the specified DataFrame column and appends them
    uncertainties += [df.at['Peak kilovoltage uncertainty', column], df.at['Anode angle uncertainty', column]]

    # Builds dictionary of beam parameters in the format required by SpekWrapper
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

    # Initialises an empty list to store transformed simulation results
    results = []

    # Iterates over each output DataFrame containing simulation results
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

    # Merges all transformed DataFrames using reduce and merge function
    merged_df = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), results)

    # Sets the columns of the merged DataFrame to match the input DataFrame columns
    merged_df.columns = input_df.columns

    # Creates a row to indicate the results section of the merged DataFrame
    data = [['Results'] + [None] * input_df.shape[1]]

    # Defines columns for the DataFrame
    columns = ['Name'] + list(input_df.columns)

    # Creates a DataFrame from data with columns specified
    df = pd.DataFrame(data, columns=columns)

    # Sets the index of the DataFrame to 'Name'
    df.set_index(keys='Name', inplace=True)

    # Concatenates the input DataFrame, the results section DataFrame and the merged results DataFrame
    return pd.concat([input_df, df, merged_df])