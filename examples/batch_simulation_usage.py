# Scripto to demonstrate the functionality of batch_simulation() function

from pathlib import Path

from uspeckpy.simulation import batch_simulation

# Using batch_simulation
# ----------------------------------------------------------------------------------------------------------------------
# TODO: units of mu_tr_rho
# TODO: check uncertainty type for np.random.normal()
# TODO: batch simulation could be a method of USpek?
# TODO: add functionalities to SpekWrapper and USpec to use with input files?

my_csv = Path('data/input/input.csv')

my_excel = Path('data/input/input.xlsx')
my_sheet = 'input'

my_folder = Path('/scratches')

df1 = batch_simulation(excel_file_path=my_excel, sheet_name=my_sheet, output_folder=my_folder)
df2 = batch_simulation(excel_file_path=my_csv, sheet_name=my_sheet, output_folder=my_folder)


