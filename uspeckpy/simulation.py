from pathlib import Path

import pandas as pd

import uspeckpy.digest as dg
from uspeckpy.uspekpy import USpek


def batch_simulation(input_file_path, sheet_name=None):
    """
    Perform batch simulations and generate output data.

    This function performs batch simulations based on input data provided in an Excel or CSV file.
    It reads the input data into a DataFrame, processes each simulation, and generates output data.
    The output data is stored in a list of DataFrames, one for each simulation, and returned. TODO

    Args:
    input_file_path (str): The path to the input file containing simulation parameters.
    sheet_name (str, optional): The name of the sheet to read if the input file is in Excel format.

    Returns:
    list: A list of DataFrames containing the output data for each simulation. TODO

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
        beam_parameters = dg.parse_beam_parameters(df=input_df, column=column_name)

        # Extract mass transmission coefficients file path from the input DataFrame
        file_path = input_df.at['Mass transmission coefficients file', column_name]

        # Parse mass transmission coefficients from the file
        mass_transmission_coefficients = dg.parse_mass_transmission_coefficients(coefficients=file_path)

        # Extract mono-energetic conversion coefficients file path from the input DataFrame
        file_path = input_df.at['Mono-energetic conversion coefficients file', column_name]

        # Extract irradiation angle from the input DataFrame
        irradiation_angle = input_df.at['Irradiation angle (deg)', column_name]

        # Parse conversion coefficients from the file
        conversion_coefficients = dg.parse_conversion_coefficients(coefficients=file_path,
                                                                   irradiation_angle=irradiation_angle)

        # Extract number of simulations from the input DataFrame
        simulations_number = input_df.at['Number of simulations', column_name]

        # Extract mass transmission coefficients uncertainty from the input DataFrame
        mass_transmission_coefficients_uncertainty = input_df.at['Mass transmission coefficients uncertainty', column_name]

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
    return dg.output_digest(input_df=input_df, output_dfs=output_dfs)


def read_file_to_dataframe(file_path, sheet_name=None):
    """
    Read data from a file into a pandas DataFrame.

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
