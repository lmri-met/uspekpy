# Scripto to demonstrate the functionality of batch_simulation() function

from pathlib import Path

from uspeckpy.simulation import batch_simulation

# TODO: units of mu_tr_rho
# TODO: batch simulation could be a method of USpek?
# TODO: add functionalities to SpekWrapper and USpec to use with input files?

# Using batch_simulation() with excel file input
# ----------------------------------------------------------------------------------------------------------------------

# Define the path to the input Excel file
my_excel = 'data/input/input.xlsx'

# Specify the name of the sheet within the Excel file to be used
my_sheet = 'input'

# Define the output folder path where the simulation results will be saved
my_folder = 'scratches'

# Call the batch_simulation function with the provided input arguments
# and store the resulting DataFrame in df1
# df1 = batch_simulation(input_file_path=my_excel, output_folder=my_folder, sheet_name=my_sheet)

# Using batch_simulation() with CSV file input
# ----------------------------------------------------------------------------------------------------------------------

# Define the path to the input CSV file
my_csv = 'data/input/input.csv'

# Call the batch_simulation function with the provided input arguments
# and store the resulting DataFrame in df1
df2 = batch_simulation(input_file_path=my_csv, output_folder=my_folder)


