# Script to calculate the mean conversion coefficient and related quantities from x-ray spectrum for
# an operational quantity and an x-ray quality
from time import time

import pandas as pd
from matplotlib import pyplot as plt


def input_digest(excel_file_path, sheet_name):
    # Print a message indicating the start of input digestion
    print('Input digest')

    # Read Excel file into a DataFrame and set 'Name' column as index
    df = pd.read_excel(excel_file_path, sheet_name=sheet_name)
    df.set_index(keys='Name', inplace=True)

    # Get the column label from the first column
    # TODO: Later on this will be a loop over all columns
    column = df.columns[0]

    # Extract numeric parameters from the DataFrame column
    quality = df.at['Quality', column]
    operational_quantity = df.at['Operational quantity', column]
    irradiation_angle = df.at['Angle (deg)', column]
    simulations_number = df.at['Number of simulations', column]
    mass_transmission_coefficients_uncertainty = df.at['Mass transmission coefficients uncertainty', column]

    # Read CSV files into DataFrames for mass transmission coefficients and mono-energetic conversion coefficients
    mass_transmission_coefficients = pd.read_csv(df.at['Mass transmission coefficients file', column])
    mono_energetic_conversion_coefficients = pd.read_csv(
        df.at['Mono-energetic conversion coefficients file', column])

    # Extract conversion coefficients corresponding to the specified irradiation angle
    column_angle = [col for col in mono_energetic_conversion_coefficients.columns if str(irradiation_angle) in col][0]
    mono_energetic_conversion_coefficients_angle = mono_energetic_conversion_coefficients[['E (keV)', column_angle]]

    # Extract values and uncertainties for various beam parameters in SpekPy format
    names = ['Al', 'Cu', 'Sn', 'Pb', 'Be', 'Air']
    values = [df.at[f'{name} filter width (mm)', column] for name in names]
    uncertainties = [df.at[f'{name} filter width uncertainty', column] for name in names]
    names += ['kVp', 'angle']
    values += [df.at['Peak kilovoltage (kV)', column], df.at['Anode angle (deg)', column]]
    uncertainties += [df.at['Peak kilovoltage uncertainty', column], df.at['Anode angle uncertainty', column]]
    beam_parameters = [(name, value, uncertainty) for name, value, uncertainty in zip(names, values, uncertainties)]

    # Return all extracted parameters as a tuple
    # TODO: function returns individual parameters (strings and numbers) and sets of parameters (list of tuples and
    #  DataFrames). Would it be better to be all dataframes?
    return (quality, operational_quantity, irradiation_angle, simulations_number, beam_parameters,
            mass_transmission_coefficients_uncertainty, mass_transmission_coefficients,
            mono_energetic_conversion_coefficients_angle)


# INPUT DATA
# ----------------------------------------------------------------------------------------------------------------------
# Define the path to the input Excel file and specify the sheet name
input_excel_file_path = './input.xlsx'
input_excel_file_sheet = 'input'

# INPUT DIGEST
# ----------------------------------------------------------------------------------------------------------------------
# Call the input_digest function to process the input Excel file and extract relevant parameters
results = input_digest(excel_file_path=input_excel_file_path, sheet_name=input_excel_file_sheet)

# OUTPUT DIGEST
# ----------------------------------------------------------------------------------------------------------------------
# Print a message indicating the beginning of the output digest
print('Output digest')

# Create a figure to plot each column of the DataFrame in subplots
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(30, 30))
df.iloc[:-3, :].plot(subplots=True, ax=axes)

# Save the figure
fig.savefig('subplots.png')

# Display the plot (only for Pycharm)
plt.show()  # Optionally display the plot

def output_digest(operational_magnitude, quality, mean_magnitudes, conversion_coefficients, execution_time,
                  output_folder):
    """
    Generate a summary of mean magnitudes and conversion coefficients and save it as both a text file and a CSV file.

    Parameters:
    - operational_magnitude (str): The operational magnitude.
    - quality (str): The x-ray radiation quality.
    - mean_magnitudes (list): A list of tuples containing mean magnitudes (value, standard deviation and percentage uncertainty).
    - conversion_coefficients (list): A list of tuples containing conversion coefficients (angle, value, standard deviation and percentage uncertainty).
    - execution_time (float): The time taken for execution in seconds.
    - output_folder (str): The directory where output files will be saved.

    Returns:
    - df (pandas.DataFrame): The DataFrame containing the summary data.

    Example:
    >>> mean_magnitudes = [(10, 0.5, 1.5), (20, 0.8, 2.0)]
    >>> conversion_coefficients = [(0, 0.1, 0.01, 0.05), (15, 0.2, 0.02, 0.08)]
    >>> output_digest('h_amb_10', 'N-60', mean_magnitudes, conversion_coefficients, 30.5, 'results')
    """

    # Create names for the data columns
    names = ['Mean energy', 'Mean kerma', 'Mean HVL1 Al', 'Mean HVL2 Al', 'Mean HVL1 Cu', 'Mean HVL2 Cu']
    names += [f'Conv. Coeff. {angle[0]} degrees' for angle in conversion_coefficients]

    # Gather values for the data columns
    values = [mean_magnitude[0] for mean_magnitude in mean_magnitudes]
    values += [conversion_coefficient[1] for conversion_coefficient in conversion_coefficients]

    # Gather standard deviations for the data columns
    standard_deviations = [mean_magnitude[1] for mean_magnitude in mean_magnitudes]
    standard_deviations += [conversion_coefficient[2] for conversion_coefficient in conversion_coefficients]

    # Gather percentage uncertainties for the data columns
    percentage_uncertainty = [mean_magnitude[2] for mean_magnitude in mean_magnitudes]
    percentage_uncertainty += [conversion_coefficient[3] for conversion_coefficient in conversion_coefficients]

    # Create a dictionary to hold the data
    data = {'Magnitude': names, 'Value': values, 'Std. Dev.': standard_deviations, '% Unc.': percentage_uncertainty}

    # Create a pandas DataFrame from the data
    df = pd.DataFrame(data)

    # Build path of output file for text file and CSV file
    file_path_txt = os.path.join(output_folder, f"{operational_magnitude}_{quality}_results.txt")
    file_path_csv = os.path.join(output_folder, f"{operational_magnitude}_{quality}_results.csv")

    # Create text file
    with open(file_path_txt, 'w') as file:
        file.write(f'Conversion coefficients for magnitude {operational_magnitude} and quality {quality}\n\n')
        file.write(df.to_string())
        file.write(f'\n\nExecution time: {execution_time} s.')

    # Save DataFrame to a CSV file
    df.to_csv(file_path_csv, index=False)

    return df