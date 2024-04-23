from pathlib import Path
from time import time

import pandas as pd

import uspeckpy.digest as dg
from uspeckpy.uspekpy import USpek


def batch_simulation(excel_file_path, sheet_name, output_folder):
    initial_time = time()

    # Print a message indicating the start of input digestion
    print('Batch simulation')

    # Print a message indicating the start of input digestion
    print('\nInitial input digest')

    # Read Excel file into a DataFrame and set 'Name' column as index
    input_df = pd.read_excel(excel_file_path, sheet_name=sheet_name)
    input_df.set_index(keys='Name', inplace=True)

    # Initialize an empty list to store the simulations results
    output_dfs = []

    # Initialize an iterator
    i = 0

    for column_name in input_df.columns:
        # Print a message indicating simulation number
        print(f'\nSimulation {i + 1}')

        # Print a message indicating the start of input digestion
        print('Input digest')

        # Get beam parameters in the format required by SpekWrapper (dictionary of tuples)
        beam_parameters = dg.parse_beam_parameters(df=input_df, column=column_name)

        # Get mass transmission coefficients in the format required by SpekWrapper (tuple of two numpy arrays)
        mass_transmission_coefficients = dg.parse_mass_transmission_coefficients(df=input_df, column=column_name)

        # Get conversion coefficients in the format required by SpekWrapper (tuple of two numpy arrays)
        conversion_coefficients = dg.parse_conversion_coefficients(df=input_df, column=column_name)

        # Extract number of simulations from the input DataFrame column
        simulations_number = input_df.at['Number of simulations', column_name]

        # Print a message indicating the start of input digestion
        print('Simulation')

        # Create USpekPy object with given beam parameters, mass transmission coefficients and conversion coefficients
        s = USpek(beam_parameters=beam_parameters, mass_transmission_coefficients=mass_transmission_coefficients,
                  conversion_coefficients=conversion_coefficients)

        # Run simulation with a given number of iterations
        output_df = s.simulate(simulations_number=simulations_number)

        output_df.to_csv(Path(output_folder) / f'output_case{i + 1}.csv', index=True)

        # Print a message indicating the start of input digestion
        print('Output digest')

        # Append output DataFrame to output DataFrames list
        output_dfs.append(output_df)

        # Increment iterator
        i += 1

    # Print a message indicating the start of input digestion
    print('\nFinal output digest')

    results = dg.output_digest(input_df=input_df, output_dfs=output_dfs)

    results.to_csv(Path(output_folder) / 'output.csv', index=True)

    print(f'\nExecution time: {time() - initial_time} s')

    return results
