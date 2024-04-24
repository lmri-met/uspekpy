from pathlib import Path

from uspeckpy.simulation import batch_simulation

# Using batch_simulation
# ----------------------------------------------------------------------------------------------------------------------

my_excel = Path('data/input/input.xlsx')
my_sheet = 'input'
my_folder = Path('/scratches')
df = batch_simulation(excel_file_path=my_excel, sheet_name=my_sheet, output_folder=my_folder)
