# Scripto to demonstrate the usage of batch_simulation() function with Excel file input
from pathlib import Path

from uspeckpy.simulation import batch_simulation

# Define the path to the input Excel file
my_excel = 'data/input/input.xlsx'

# Define the name of the sheet in the input Excel file (it is a must for Excel input files)
my_sheet = 'input'

# Define the output folder path where the simulation results will be saved
my_folder = 'output'

# Call the batch_simulation function with the provided input arguments
# and store the resulting DataFrame in df1
df = batch_simulation(input_file_path=my_excel, output_folder=my_folder, sheet_name=my_sheet)

# Save results to a CSV file
df.to_csv(Path(my_folder) / 'output.csv', index=True)