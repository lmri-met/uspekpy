from uspeckpy.main import main


def compare_text_files(file1_path, file2_path):
    with open(file1_path, 'r') as file1:
        content1 = file1.read()

    with open(file2_path, 'r') as file2:
        content2 = file2.read()

    if content1 == content2:
        return True
    else:
        return False


# H*(10)
# ----------------------------------------------------------------------------------------------------------------------

main(beam_data_file='../examples/h_amb/inputN60_no_us_@1000mm_TEST.xlsx',
     conversion_coefficients_files=['../examples/h_amb/h_amb_10.csv'],
     transmission_coefficients_file='../examples/h_amb/mutr.txt',
     transmission_coefficients_uncertainty=0,
     simulations_number=100,
     output_folder='../examples/output')

# Hp(0.07, slab)
# ----------------------------------------------------------------------------------------------------------------------

main(beam_data_file='../examples/hp_007_slab/inputN60_no_us_@2500mm_TEST.xlsx',
     conversion_coefficients_files=['../examples/hp_007_slab/hp_0.07_slab.csv'],
     transmission_coefficients_file='../examples/hp_007_slab/mutr.txt',
     transmission_coefficients_uncertainty=0,
     simulations_number=100,
     output_folder='../examples/output')

# Hp(3, cyl)
# ----------------------------------------------------------------------------------------------------------------------

main(beam_data_file='../examples/hp_3_cyl/inputN60_no_us_@2500mm_TEST.xlsx',
     conversion_coefficients_files=['../examples/hp_3_cyl/hp_3_cyl.csv'],
     transmission_coefficients_file='../examples/hp_3_cyl/mutr.txt',
     transmission_coefficients_uncertainty=0,
     simulations_number=100,
     output_folder='../examples/output')

# H'(0.07)
# ----------------------------------------------------------------------------------------------------------------------

main(beam_data_file='../examples/h_prime_007/inputN60_no_us_@2500mm_TEST.xlsx',
     conversion_coefficients_files=['../examples/h_prime_007/h_prime_0.07.csv'],
     transmission_coefficients_file='../examples/h_prime_007/mutr.txt',
     transmission_coefficients_uncertainty=0,
     simulations_number=100,
     output_folder='../examples/output')

# print(compare_text_files('examples/ref_output/h_amb_10.csv_N-60_resultados.txt',
#                          'examples/output/h_amb_10.csv_N-60_resultados.txt'))
# # Files content differ but numbers are equal.
# print(compare_text_files('examples/output/h_amb_10.csv_N-60_resultados_.txt',
#                          'examples/output/h_amb_10.csv_N-60_resultados.txt'))
# # Identical files
# print(compare_text_files('examples/output/h_prime_0.07.csv_N-60_resultados_.txt',
#                          'examples/output/h_prime_0.07.csv_N-60_resultados.txt'))
# # Identical files
# print(compare_text_files('examples/output/hp_0.07_slab.csv_N-60_resultados_.txt',
#                          'examples/output/hp_0.07_slab.csv_N-60_resultados.txt'))
# # Identical files
# print(compare_text_files('examples/output/hp_3_cyl.csv_N-60_resultados_.txt',
#                          'examples/output/hp_3_cyl.csv_N-60_resultados.txt'))
# # Identical files
