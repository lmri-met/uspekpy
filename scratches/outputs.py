import os
from random import random

import pandas as pd


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


# Assign random values to variables representing different magnitudes
media_energia = random()
desviacion_energia = random()
coef_variacion_energia = random()
media_kerma = random()
desviacion_kerma = random()
coef_variacion_kerma = random()
mean_hvl = random()
sd_hvl = random()
v_hvl = random()
mean_hvl2 = random()
sd_hvl2 = random()
v_hvl2 = random()
mean_hvlCu = random()
sd_hvlCu = random()
v_hvlCu = random()
mean_hvl2Cu = random()
sd_hvl2Cu = random()
v_hvl2Cu = random()

# Assign random values to variables representing hpk at different angles
hpk, hpk15, hpk30, hpk45, hpk60, hpk75, hpk90, hpk180 = (random(),) * 8
sd_hpk, sd_hpk_15, sd_hpk_30, sd_hpk_45, sd_hpk_60, sd_hpk_75, sd_hpk_90, sd_hpk_180 = (random(),) * 8
v_hpk, v_hpk_15, v_hpk_30, v_hpk_45, v_hpk_60, v_hpk_75, v_hpk_90, v_hpk_180 = (random(),) * 8

# List of magnitudes independent of the angle
magnitudes = [(media_energia, desviacion_energia, coef_variacion_energia),
              (media_kerma, desviacion_kerma, coef_variacion_kerma),
              (mean_hvl, sd_hvl, v_hvl),
              (mean_hvl2, sd_hvl2, v_hvl2),
              (mean_hvlCu, sd_hvlCu, v_hvlCu),
              (mean_hvl2Cu, sd_hvl2Cu, v_hvl2Cu)]
print(f'magnitudes: {magnitudes}')

# List of angles and their respective values and uncertainties of hpk
angles = [(0, hpk, sd_hpk, v_hpk),
          (15, hpk15, sd_hpk_15, v_hpk_15),
          (30, hpk30, sd_hpk_30, v_hpk_30),
          (45, hpk45, sd_hpk_45, v_hpk_45),
          (60, hpk60, sd_hpk_60, v_hpk_60),
          (75, hpk75, sd_hpk_75, v_hpk_75),
          (90, hpk90, sd_hpk_90, v_hpk_90),
          (180, hpk180, sd_hpk_180, v_hpk_180)]
print(f'angles: {angles}')

example_df = output_digest(operational_magnitude='h_amb_10', quality='N-60', mean_magnitudes=magnitudes,
                           conversion_coefficients=angles, execution_time=120, output_folder='./')
