# Scripto to demonstrate the usage of batch_simulation() function with CSV file input
from pathlib import Path

from uspeckpy.simulation import batch_simulation

# Define the path to the input CSV file
my_csv = 'data/input.csv'

# Define the output folder path where the simulation results will be saved
my_folder = 'output'

# Define the name of the output file where the simulation results will be saved
my_output_csv = 'output.csv'

# Call the batch_simulation function with the provided input arguments
# and store the resulting DataFrame in df
df = batch_simulation(input_file_path=my_csv)

# Save results to a CSV file
df.to_csv(Path(my_folder) / my_output_csv, index=True)
