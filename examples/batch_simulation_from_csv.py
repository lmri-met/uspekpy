# Scripto to demonstrate the usage of batch_simulation() function with CSV file input
from pathlib import Path

from uspeckpy.simulation import batch_simulation

# Define the path to the input CSV file
my_csv = 'data/input/input.csv'

# Define the output folder path where the simulation results will be saved
my_folder = 'output'

# Call the batch_simulation function with the provided input arguments
# and store the resulting DataFrame in df
df = batch_simulation(input_file_path=my_csv, output_folder=my_folder)

# Save results to a CSV file
df.to_csv(Path(my_folder) / 'output.csv', index=True)
