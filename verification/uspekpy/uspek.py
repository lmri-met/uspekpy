# Scripto to demonstrate the usage of batch_simulation() function with Excel file input
from pathlib import Path

from uspeckpy.simulation import batch_simulation

# Define the path to the input Excel file
my_excel = 'input/input_h60_hp_07_slab_45.xlsx'

# Define the name of the sheet in the input Excel file (it is a must for Excel input files)
my_sheet = 'input'

# Define the output folder path where the simulation results will be saved
my_folder = '.'

# Define the name of the output file where the simulation results will be saved
my_output_csv = 'output.csv'

# Call the batch_simulation function with the provided input arguments
# and store the resulting DataFrame in df
df = batch_simulation(input_file_path=my_excel, sheet_name=my_sheet)

# Save results to a CSV file
df.to_csv(Path(my_folder) / my_output_csv, index=True)
