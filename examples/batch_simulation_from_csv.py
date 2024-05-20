# Scripto to demonstrate the usage of batch_simulation() function with CSV file input
from uspekpy import batch_simulation

# Define the path to the input CSV file
my_csv = 'data/input.csv'

# Call the batch_simulation function with the defined input CSV file
df = batch_simulation(input_file_path=my_csv)

# Define the path of the output file where the simulation results will be saved
my_output_csv = 'output/output.csv'

# Save results to a CSV file
df.to_csv(my_output_csv, index=True)

# Print results
print(df.to_string())
