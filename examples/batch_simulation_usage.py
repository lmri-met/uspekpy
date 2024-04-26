# Scripto to demonstrate the functionality of batch_simulation() function
from pathlib import Path

from uspeckpy.simulation import batch_simulation

# Using batch_simulation() with excel file input
# ----------------------------------------------------------------------------------------------------------------------

# Define the path to the input Excel file
my_excel = 'data/input/input.xlsx'

# Define the name of the sheet in the input Excel file (it is a must for Excel input files)
my_sheet = 'input'

# Define the output folder path where the simulation results will be saved
my_folder = 'output'

# Call the batch_simulation function with the provided input arguments
# and store the resulting DataFrame in df1
df1 = batch_simulation(input_file_path=my_excel, output_folder=my_folder, sheet_name=my_sheet)

# Save results to a CSV file
df1.to_csv(Path(my_folder) / 'output.csv', index=True)

# Using batch_simulation() with CSV file input
# ----------------------------------------------------------------------------------------------------------------------

# Define the path to the input CSV file
my_csv = 'data/input/input.csv'

# Define the output folder path where the simulation results will be saved
my_folder = 'output'

# Call the batch_simulation function with the provided input arguments
# and store the resulting DataFrame in df1
df2 = batch_simulation(input_file_path=my_csv, output_folder=my_folder)

# Save results to a CSV file
df2.to_csv(Path(my_folder) / 'output.csv', index=True)
