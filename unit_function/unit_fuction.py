# Script to calculate the mean conversion coefficient and related quantities from x-ray spectrum for
# an operational quantity and an x-ray quality
from time import time
import spekpy as sp
import pandas as pd
import numpy as np
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
input_excel_file_path = 'C:/Users/admin/PycharmProjects/uspekpy/unit_function/input.xlsx'
input_excel_file_sheet = 'input'

# INPUT DIGEST
# ----------------------------------------------------------------------------------------------------------------------
# Measure the initial time before executing the input digestion process
initial_time = time()

# Call the input_digest function to process the input Excel file and extract relevant parameters
results = input_digest(excel_file_path=input_excel_file_path, sheet_name=input_excel_file_sheet)

# COMPUTE MEAN CONVERSION COEFFICIENTS FOR A GIVEN X-RAY QUALITY, OPERATIONAL QUANTITY AND IRRADIATION ANGLE
# ----------------------------------------------------------------------------------------------------------------------
# Print a message indicating the beginning of the calculation of mean conversion coefficient
print('Calculate mean conversion coefficient')

# Extract relevant data from the results obtained from input_digest function
random_beam_parameters = results[4]
mass_transmission_coefficients = results[6]
mono_energetic_conversion_coefficients = results[7]
simulations_number = results[3]
# simulations_number = 100  # TODO: remove

# Extract beam parameters for further processing
names = [parameter[0] for parameter in random_beam_parameters]
values = [parameter[1] for parameter in random_beam_parameters]
uncertainties = [parameter[2] for parameter in random_beam_parameters]

# Create a DataFrame for beam parameters
beam_df = pd.DataFrame({'Names': names, 'Value': values, 'Uncertainty': uncertainties})
beam_df.set_index(keys='Names', inplace=True)

# Calculate minimum and maximum values for beam parameters (filters width) for later Monte Carlo calculations
beam_df['Minimum'] = beam_df.loc[:'Be', 'Value'] * (1 - beam_df.loc[:'Be', 'Uncertainty'] * np.sqrt(3))
beam_df['Maximum'] = beam_df.loc[:'Be', 'Value'] * (1 + beam_df.loc[:'Be', 'Uncertainty'] * np.sqrt(3))

# Take the natural logarithm of mass transmission coefficients and mono-energetic conversion coefficients for later
# interpolation
for column in mass_transmission_coefficients.columns:
    mass_transmission_coefficients[f'ln({column})'] = np.log(mass_transmission_coefficients[column])
for column in mono_energetic_conversion_coefficients.columns:
    mono_energetic_conversion_coefficients[f'ln({column})'] = np.log(mono_energetic_conversion_coefficients[column])

# Initialize an empty DataFrame for storing simulation results
columns = ['#', 'Al (mm)', 'Cu (mm)', 'Sn (mm)', 'Pb (mm)', 'Be (mm)', 'Air (mm)', 'kVp (kV)', 'angle (deg)',
           'HVL1 Al (mm)', 'HVL2 Al (mm)', 'HVL1 Cu (mm)', 'HVL2 Cu (mm)']
df = pd.DataFrame(columns=columns)

# Perform simulations for the specified number of iterations
for i in range(simulations_number):
    # Print a message indicating the number of the current iteration
    print(f'Iteration number: {i}')

    # Initialize an empty list for storing random values of beam parameters in SpekPy format
    random_beam_parameters = []

    # Generate random values for beam parameters: filters width, uniform distribution
    for index in ['Al', 'Cu', 'Sn', 'Pb', 'Be']:
        random_value = np.random.uniform(beam_df.loc[index, 'Minimum'], beam_df.loc[index, 'Maximum'], 1)
        random_beam_parameters.append((index, float(random_value[0])))

    # Generate random values for beam parameters: air width, peak kilovoltage and anode angle, normal distribution
    for index in ['Air', 'kVp', 'angle']:
        random_value = np.random.normal(beam_df.loc[index, 'Value'], beam_df.loc[index, 'Uncertainty'], 1)
        random_beam_parameters.append((index, float(random_value[0])))

    # Generate spectrum based on the random beam parameters with SpekPy
    spectrum = sp.Spek(kvp=random_beam_parameters[6][1], th=random_beam_parameters[7][1])
    spectrum.multi_filter(random_beam_parameters[:6])
    energy, fluence = spectrum.get_spectrum(edges=False)

    # Calculate half-value layers for aluminum and copper with SpekPy
    hvl1_al = spectrum.get_hvl1()
    hvl2_al = spectrum.get_hvl2()
    hvl1_cu = spectrum.get_hvl1(matl='Cu')
    hvl2_cu = spectrum.get_hvl2(matl='Cu')

    # Store spectrum in a DataFrame
    spectrum_df = pd.DataFrame({'E (keV)': energy, 'Fluence ()': fluence})

    # Calculate mean energy
    # Calculate mean air kerma
    # Take logarithms for the spectrum for later interpolation?
    # Interpolate mass transmission coefficients to the energy bins of the spectrum using Akima method
    # Interpolate conversion coefficients to the energy bins of the spectrum using Akima method

    # Store simulation results in a DataFrame
    row = [f'Iteration {i + 1}'] + [parameter[1] for parameter in random_beam_parameters] + [hvl1_al, hvl2_al, hvl1_cu, hvl2_cu]
    df = pd.concat(objs=[df, pd.DataFrame(data=[row], columns=columns)], ignore_index=True)

# Calculate means, standard deviations, and relative uncertainties for the simulation results
means = df.iloc[:, 1:].mean()
standard_deviations = df.iloc[:, 1:].std(ddof=0)
relative_uncertainties = standard_deviations / means

# Append means, standard deviations, and relative uncertainties to the DataFrame
means = ['Mean'] + list(means)
standard_deviations = ['Standard deviation'] + list(standard_deviations)
relative_uncertainties = ['Relative uncertainty'] + list(relative_uncertainties)
data = [means, standard_deviations, relative_uncertainties]
df = pd.concat(objs=[df, pd.DataFrame(data=data, columns=columns)], ignore_index=True)

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

# EXECUTION TIME
# ----------------------------------------------------------------------------------------------------------------------
# Measure the final time after completing the calculations
final_time = time()

# Calculate the execution time
execution_time = final_time - initial_time

# Print the execution time
print(f'Execution time: {execution_time} s')
