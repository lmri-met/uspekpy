from functools import reduce

import numpy as np
import pandas as pd


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


def parse_mass_transmission_coefficients(df, column):
    # Get mass transmission coefficients in the format required by SpekWrapper (tuple of two numpy arrays)
    # The mass transmission coefficients CSV file should have only 2 columns, the first one for the energy and the
    # second one for the mass transmission coefficients.

    # Extract mass transmission coefficients CSV file path from input DataFrame column
    file_path = df.at['Mass transmission coefficients file', column]

    # Load CSV file into numpy array
    array2d = np.genfromtxt(file_path, delimiter=',', skip_header=1, unpack=True)

    # Build tuple of mass transmission coefficients in the format required by SpekWrapper (tuple of two numpy arrays)
    return array2d[0], array2d[1]


def parse_conversion_coefficients(df, column):
    # Get conversion coefficients in the format required by SpekWrapper (tuple of two numpy arrays)

    # Extract conversion coefficients CSV file path from input DataFrame column
    file_path = df.at['Mono-energetic conversion coefficients file', column]

    # Read CSV files into DataFrames for mono-energetic conversion coefficients
    df = pd.read_csv(file_path)

    # Get the energies of the mono-energetic conversion coefficients (always the first column of the Dataframe)
    energies = df.iloc[:, 0]

    # If the conversion coefficients Dataframe has only 2 columns, the first is the energy and the second is the
    # conversion coefficients.
    if df.shape[1] == 2:
        # Get the values of the mono-energetic conversion coefficients (the second column of the Dataframe)
        values = df.iloc[:, 1]
    # If the conversion coefficients Dataframe has more than 2 columns, then it contains conversion coefficients for
    # more than one irradiation angle.
    else:  # TODO: check with multicolumn input file
        # Extract irradiation angle from the input DataFrame column
        irradiation_angle = df.at['Irradiation angle (deg)', column]

        # Get the label of the column which contains the irradiation angle
        column_label = next((label for label in df.columns if irradiation_angle in label), None)

        # Get the values of the mono-energetic conversion coefficients (the second column of the Dataframe)
        values = df.loc[:, column_label]

    # Build tuple of conversion coefficients in the format required by SpekWrapper (tuple of two numpy arrays)
    return energies, values


def output_digest(input_df, output_dfs):
    # Build a DataFrame with the simulations result in the format of the input DataFrame
    result_columns = ['HVL1 Al', 'HVL2 Al', 'HVL1 Cu', 'HVL2 Cu', 'Mean energy', 'Mean kerma',
                      'Mean conv. coefficient.']
    result_rows = ['Mean', 'Standard deviation', 'Relative uncertainty']
    results = []
    for output_df in output_dfs:
        # Set the index of the DataFrame to a column named '#', modifying the DataFrame in place
        output_df.set_index(keys='#', inplace=True)

        # Extract a subset of data from the DataFrame using specific rows and columns
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

        results.append(df)

    # Merge the results DataFrames on their indices using reduce
    merged_df = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), results)
    merged_df.columns = input_df.columns
    # df = pd.DataFrame(data=rows2)

    # Create DataFrame with one row using a list
    data = [['Results'] + [None] * input_df.shape[1]]
    columns = ['Name'] + list(input_df.columns)
    df = pd.DataFrame(data, columns=columns)
    df.set_index(keys='Name', inplace=True)

    # Concatenate along rows
    return pd.concat([input_df, df, merged_df])
