# Scripto to demonstrate the usage of batch_simulation() function with Excel file input
from uspekpy import batch_simulation

# Define the path to the input Excel file
my_excel = 'data/input.xlsx'

# Define the name of the sheet in the input Excel file
my_sheet = 'input'

# Call the batch_simulation function with the defined input Excel file and sheet
df = batch_simulation(input_file_path=my_excel, sheet_name=my_sheet)

# Define the path of the output file where the simulation results will be saved
my_output_csv = 'output/output.csv'

# Save results to a CSV file
df.to_csv(my_output_csv, index=True)

# Print results
print(df.to_string())
