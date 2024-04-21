# Script to calculate the mean conversion coefficient and related quantities from x-ray spectrum for
# an operational quantity and an x-ray quality
from time import time
import spekpy as sp
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


def input_digest(excel_file_path, sheet_name):
    print('Input digest')
    # Read input file to a dataframe
    df = pd.read_excel(excel_file_path, sheet_name=sheet_name)
    # Set first column as dataframe index
    df.set_index(keys='Name', inplace=True)
    # Select column to read: the first one. TODO: Later on this will be a loop over all columns
    column = df.columns[0]

    # Read numerical data from selected column
    quality = df.at['Quality', column]
    operational_quantity = df.at['Operational quantity', column]
    irradiation_angle = df.at['Angle (deg)', column]
    simulations_number = df.at['Number of simulations', column]
    mass_transmission_coefficients_uncertainty = df.at['Mass transmission coefficients uncertainty', column]

    # Read CSV files into DataFrames
    mass_transmission_coefficients = pd.read_csv(df.at['Mass transmission coefficients file', column])
    mono_energetic_conversion_coefficients = pd.read_csv(
        df.at['Mono-energetic conversion coefficients file', column])

    # Select the column of the conversion coefficients CSV which contains the selected angle
    column_angle = [col for col in mono_energetic_conversion_coefficients.columns if str(irradiation_angle) in col][0]
    mono_energetic_conversion_coefficients_angle = mono_energetic_conversion_coefficients[['E (keV)', column_angle]]

    # Read beam parameters (values and uncertainties) in SpekPy format
    names = ['Al', 'Cu', 'Sn', 'Pb', 'Be', 'Air']
    values = [df.at[f'{name} filter width (mm)', column] for name in names]
    uncertainties = [df.at[f'{name} filter width uncertainty', column] for name in names]
    names += ['kVp', 'angle']
    values += [df.at['Peak kilovoltage (kV)', column], df.at['Anode angle (deg)', column]]
    uncertainties += [df.at['Peak kilovoltage uncertainty', column], df.at['Anode angle uncertainty', column]]
    beam_parameters = [(name, value, uncertainty) for name, value, uncertainty in zip(names, values, uncertainties)]

    return (quality, operational_quantity, irradiation_angle, simulations_number, beam_parameters,
            mass_transmission_coefficients_uncertainty, mass_transmission_coefficients,
            mono_energetic_conversion_coefficients_angle)


# Input data
# ----------------------------------------------------------------------------------------------------------------------
input_excel_file_path = 'C:/Users/admin/PycharmProjects/uspekpy/unit_function/input.xlsx'
input_excel_file_sheet = 'input'

# Measure initial time
initial_time = time()

# Input digest
# ----------------------------------------------------------------------------------------------------------------------
results = input_digest(excel_file_path=input_excel_file_path, sheet_name=input_excel_file_sheet)

# Compute mean conversion coefficient for the given x-ray quality, operational quantity and irradiation angle
# ----------------------------------------------------------------------------------------------------------------------
print('Calculate mean conversion coefficient')
# Unpack input digest results
random_beam_parameters = results[4]
mass_transmission_coefficients = results[6]
mono_energetic_conversion_coefficients = results[7]
simulations_number = results[3]
#simulations_number = 100  # TODO: remove

# Build dataframe from filters
names = [parameter[0] for parameter in random_beam_parameters]
values = [parameter[1] for parameter in random_beam_parameters]
uncertainties = [parameter[2] for parameter in random_beam_parameters]
beam_df = pd.DataFrame({'Names': names, 'Value': values, 'Uncertainty': uncertainties})
beam_df.set_index(keys='Names', inplace=True)

# Calculate upper and lower limits of filters width for later Monte Carlo calculations
beam_df['Minimum'] = beam_df.loc[:'Be', 'Value'] * (1 - beam_df.loc[:'Be', 'Uncertainty'] * np.sqrt(3))
beam_df['Maximum'] = beam_df.loc[:'Be', 'Value'] * (1 + beam_df.loc[:'Be', 'Uncertainty'] * np.sqrt(3))

# Take logarithms for mono energetic conversion coefficients and mass transmission coefficients for later interpolation?
for column in mass_transmission_coefficients.columns:
    mass_transmission_coefficients[f'ln({column})'] = np.log(mass_transmission_coefficients[column])
for column in mono_energetic_conversion_coefficients.columns:
    mono_energetic_conversion_coefficients[f'ln({column})'] = np.log(mono_energetic_conversion_coefficients[column])

# Define empty dataframe to store iteration results
columns = ['#', 'Al (mm)', 'Cu (mm)', 'Sn (mm)', 'Pb (mm)', 'Be (mm)', 'Air (mm)', 'kVp (kV)', 'angle (deg)',
           'HVL1 Al (mm)', 'HVL2 Al (mm)', 'HVL1 Cu (mm)', 'HVL2 Cu (mm)']
df = pd.DataFrame(columns=columns)

# Iterate over the number of simulations
for i in range(simulations_number):
    print(f'Iteration number: {i}')

    # Define empty list to store beam parameters in SpekPy format
    random_beam_parameters = []

    # Get random values for filters width (uniform distribution)
    for index in ['Al', 'Cu', 'Sn', 'Pb', 'Be']:
        random_value = np.random.uniform(beam_df.loc[index, 'Minimum'], beam_df.loc[index, 'Maximum'], 1)
        random_beam_parameters.append((index, float(random_value[0])))

    # Get random values for air width, peak kilovoltage and anode angle  (normal distribution)
    for index in ['Air', 'kVp', 'angle']:
        random_value = np.random.normal(beam_df.loc[index, 'Value'], beam_df.loc[index, 'Uncertainty'], 1)
        random_beam_parameters.append((index, float(random_value[0])))

    # Calculate the spectrum with SpekPy
    spectrum = sp.Spek(kvp=random_beam_parameters[6][1], th=random_beam_parameters[7][1])
    spectrum.multi_filter(random_beam_parameters[:6])
    energy, fluence = spectrum.get_spectrum(edges=False)

    # Calculate first and second HVL for Al and Cu  with SpekPy
    hvl1_al = spectrum.get_hvl1()
    hvl2_al = spectrum.get_hvl2()
    hvl1_cu = spectrum.get_hvl1(matl='Cu')
    hvl2_cu = spectrum.get_hvl2(matl='Cu')

    # Store spectrum in dataframe
    spectrum_df = pd.DataFrame({'E (keV)': energy, 'Fluence ()': fluence})
    # Calculate mean energy
    # Calculate mean air kerma
    # Take logarithms for the spectrum for later interpolation?
    # Interpolate mass transmission coefficients to the energy bins of the spectrum using Akima method
    # Interpolate conversion coefficients to the energy bins of the spectrum using Akima method

    # Append iteration results to storage dataframe
    row = [f'Iteration {i + 1}'] + [parameter[1] for parameter in random_beam_parameters] + [hvl1_al, hvl2_al, hvl1_cu, hvl2_cu]
    df = pd.concat(objs=[df, pd.DataFrame(data=[row], columns=columns)], ignore_index=True)

# Compute the mean of each column except the first one
means = df.iloc[:, 1:].mean()

# Compute the standard deviation of each column except the first one with ddof=0
standard_deviations = df.iloc[:, 1:].std(ddof=0)

# Compute the relative uncertainty by dividing standard deviations by means
relative_uncertainties = standard_deviations / means

# Convert means, standard deviations, and relative uncertainties to lists and add labels
means = ['Mean'] + list(means)
standard_deviations = ['Standard deviation'] + list(standard_deviations)
relative_uncertainties = ['Relative uncertainty'] + list(relative_uncertainties)

# Combine the means, standard deviations, and relative uncertainties into a list of lists
data = [means, standard_deviations, relative_uncertainties]

# Concatenate the list of lists with the original DataFrame to append the computed statistics as new rows
df = pd.concat(objs=[df, pd.DataFrame(data=data, columns=columns)], ignore_index=True)

# Output digest
# ----------------------------------------------------------------------------------------------------------------------
print('Output digest')
# Plot each column in a separate subplot
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(30, 30))
df.iloc[:-3, :].plot(subplots=True, ax=axes)

# Save all subplots in a single figure
fig.savefig('subplots.png')

plt.show()  # Optionally display the plot

# Calculate execution time
# ----------------------------------------------------------------------------------------------------------------------
# Measure initial time
final_time = time()

# Calculate execution time
execution_time = final_time - initial_time

print(f'Execution time: {execution_time} s')
